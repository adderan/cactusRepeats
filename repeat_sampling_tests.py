import argparse
import re
import os
import shutil
import xml.etree.ElementTree

from toil.common import Toil
from toil.job import Job
from sonLib.bioio import popenCatch
from sonLib.bioio import getTempFile
from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import system


def toURL(path):
    if path.startswith("http://"):
        return path
    else:
        return "file://" + os.path.abspath(path)

def makeRawSeedCounts(job, options, sequenceID):
    job.fileStore.logToMaster("Getting seed counts from lastz...")
    sequence = job.fileStore.readGlobalFile(sequenceID)
    rawSeedCounts = job.fileStore.getLocalTempFile()
    os.system("cPecanLastz --tableonly=count %s[unmask][multiple] > %s" % (sequence, rawSeedCounts))
    countsID =  job.fileStore.writeGlobalFile(rawSeedCounts)
    return countsID

def makeSeedProbabilities(job, options, sequenceID, countsTableID, graphID):
    countsTable = job.fileStore.readGlobalFile(countsTableID)
    sequence = job.fileStore.readGlobalFile(sequenceID)
    deBruijnGraph = job.fileStore.readGlobalFile(graphID)
    if options.reuseProbabilitiesFile:
        return job.fileStore.writeGlobalFile(options.seedProbabilitiesFilename)
    seedProbabilities = job.fileStore.getLocalTempFile()

    flags = ""
    if options.clusterSeeds:
        flags += "--clusterSeeds"
    if options.seedMultiplicity:
        flags += "--seedMultiplicity"
    if options.pathMultiplicity:
        flags += "--pathMultiplicity"

    os.system("cactus_repeats_seedProbabilities %s --seedMultiplicitiy --sequenceFile %s --graphFile %s --lastzFile %s > %s"
            % (flags, sequence, deBruijnGraph, countsTable, seedProbabilities))
    seedProbabilitiesID = job.fileStore.writeGlobalFile(seedProbabilities)
    if options.seedProbabilitiesFilename and not options.reuseProbabilitiesFile:
        job.fileStore.exportFile(seedProbabilitiesID, toURL(options.seedProbabilitiesFilename))
    return seedProbabilitiesID

def makeDeBruijnGraph(job, options, sequenceID):
    job.fileStore.logToMaster("Making de bruijn graph for input genome...")
    sequence = job.fileStore.readGlobalFile(sequenceID)
    tmpDir = job.fileStore.getLocalTempDir()
    os.system("twopaco -f %i %s -k %i --tmpdir %s" % (options.filterSize, sequence, 
        options.k, tmpDir))
    return job.fileStore.writeGlobalFile(os.path.join(tmpDir, "de_bruijn.bin"))


def makeSeedProbabilitiesFromDeBruinPaths(job, options, deBruijnID, sequenceID):
    sequence = job.fileStore.readGlobalFile(sequenceID)
    deBruijn = job.fileStore.readGlobalFile(deBruijnID)
    seedProbabilitiesFile = job.fileStore.getLocalTempFile()
    os.system("cactus_repeats_seedMultiplicity --graphFilename %s --sequenceFilename %s --k %i > %s" % (deBruijn, sequence, options.k))
    return job.fileStore.writeGlobalFile(seedProbabilitiesFile)

def trimSequence(job, options, sequenceID, start, stop):
    if options.region is None or start is None or stop is None:
        return sequenceID
    sequenceFile = job.fileStore.readGlobalFile(sequenceID)
    trimmedSequenceFile = job.fileStore.getLocalTempFile()
    system("samtools faidx %s" % sequenceFile)
    system("samtools faidx %s %s:%i-%i > %s" % (sequenceFile, options.region, start, stop, trimmedSequenceFile))
    return job.fileStore.writeGlobalFile(trimmedSequenceFile)

def parseStats(statsFile):
    """Parse the lastz stats file.
    """

    stats = dict()
    for line in open(statsFile, 'r'):
        line = line.strip()
        colon = line.find(":")
        if not colon:
            continue
        statName = line[:colon]
        statValue = line[colon+1:]
        try:
            parsedStat = int(statValue.replace(",", ""))
            stats[statName] = parsedStat
        except ValueError:
            continue

    return stats

def runLastz(job, options, seqID1, seqID2):
    job.fileStore.logToMaster("Running LASTZ without seed sampling...")
    seq1 = job.fileStore.readGlobalFile(seqID1)
    seq2 = job.fileStore.readGlobalFile(seqID2)
    alignmentsFile = job.fileStore.getLocalTempFile()
    statsFile = job.fileStore.getLocalTempFile()
    os.system("cPecanLastz %s[multiple] %s[multiple] --format=maf --stats=%s > %s" % (seq1, seq2, statsFile, alignmentsFile))
    stats = parseStats(statsFile)
    return (job.fileStore.writeGlobalFile(alignmentsFile), stats)


def runSampledLastz(job, options, seqID1, seqID2, seedProbabilitiesID):
    seq1 = job.fileStore.readGlobalFile(seqID1)
    seq2 = job.fileStore.readGlobalFile(seqID2)
    seedProbabilities = job.fileStore.readGlobalFile(seedProbabilitiesID)
    alignmentsFile = job.fileStore.getLocalTempFile()
    statsFile = job.fileStore.getLocalTempFile()
    cmd = "cPecanLastz %s[multiple][unmask] %s[multiple][unmask] --notrivial --format=maf --sampleSeeds --seedSamplingProbabilities=%s --stats=%s > %s" % (seq1, seq2, seedProbabilities, statsFile, alignmentsFile)
    job.fileStore.logToMaster("Running lastz command: %s" % cmd)
    messages = system(cmd)
    stats = parseStats(statsFile)
    job.fileStore.logToMaster("Number of HSPs = %i" % stats["HSPs"])
    job.fileStore.logToMaster("Seeds skipped = %i" % stats["Seeds skipped by sampling"])
    return (job.fileStore.writeGlobalFile(alignmentsFile), stats)

def alignmentPrecisionRecall(job, options, maskedAlignmentsID, sampledAlignmentsID):
    maskedAlignments = job.fileStore.readGlobalFile(maskedAlignmentsID)
    sampledAlignments = job.fileStore.readGlobalFile(sampledAlignmentsID)
    precisionRecall = job.fileStore.getLocalTempFile()
    system("mafComparator --maf1=%s --maf2=%s --out=%s" % (maskedAlignments, sampledAlignments, precisionRecall))
    precisionRecallID = job.fileStore.writeGlobalFile(precisionRecall)
    precisionRecallXML = xml.etree.ElementTree.parse(precisionRecall).getroot()
    precision = float(precisionRecallXML[0][0][0].attrib["average"])
    recall = float(precisionRecallXML[1][0][0].attrib["average"])
    job.fileStore.logToMaster("Precision = %f, recall = %f" % (precision, recall))
    if options.precisionRecallFile:
        job.fileStore.exportFile(precisionRecallID, toURL(options.precisionRecallFile))
    return (precision, recall)

def printThresholdStats(job, options, thresholdValues, statsDicts):
    job.fileStore.logToMaster("Writing output file...")
    outputFile = job.fileStore.getLocalTempFile()
    with open(outputFile, 'w') as fh:
        for threshold, stats in zip(thresholdValues, statsDicts):
            fh.write("%i\t%i\n" % (threshold, stats["HSPs"]))
    outputID = job.fileStore.writeGlobalFile(outputFile)
    job.fileStore.exportFile(outputID, toURL(options.outputFile))

def printScalabilityStats(job, options, statsList):
    outputFile = job.fileStore.getLocalTempFile()
    with open(outputFile, 'w') as fh:
        for stats in statsList:
            fh.write("%i %i\n" % (stats["target length"], stats["HSPs"]))
    outputID = job.fileStore.writeGlobalFile(outputFile)
    job.fileStore.exportFile(outputID, toURL(options.outputFile))

def printNumberOfHSPs(job, options, stats):
    job.fileStore.logToMaster("HSPs = %i" % stats["HSPs"])
def dummyJobFn(job):
    pass

def plotScalability(options):
    nPoints = 20
    with Toil(options) as toil:
        sequenceID = toil.importFile(toURL(options.sequence))
        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        seedProbabilitiesJob = Job.wrapJobFn(makeSeedProbabilities, options, sequenceID)
        dummyJob = Job.wrapJobFn(dummyJobFn)
        trimSeqJobs = []
        lastzJobs = []
        statsList = []
        seqStep = options.seqLength/nPoints
        for i in range(1, nPoints):
            trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart + i*seqStep)
            if options.sampleSeeds:
                lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedProbabilitiesJob.rv())
            else:
                lastzJob = Job.wrapJobFn(runLastz, options, trimSeqJob.rv(), trimSeqJob.rv())
            statsList.append(lastzJob.rv(1))
            lastzJobs.append(lastzJob)
            trimSeqJobs.append(trimSeqJob)

        printStatsJob = Job.wrapJobFn(printScalabilityStats, options, statsList)
        seedCountsJob.addFollowOn(seedProbabilitiesJob)
        seedProbabilitiesJob.addFollowOn(dummyJob)
        for trimSeqJob, lastzJob in zip(trimSeqJobs, lastzJobs):
            dummyJob.addChild(trimSeqJob)
            trimSeqJob.addFollowOn(lastzJob)
        dummyJob.addFollowOn(printStatsJob)

        toil.start(seedCountsJob)

def runWorkflow(options):
    with Toil(options) as toil:
        sequenceID = toil.importFile(toURL(options.sequence))
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart + options.seqLength)

        if options.sampleSeeds:
            seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
            deBruijnGraphJob = Job.wrapJobFn(makeDeBruijnGraph, options, sequenceID)
            trimSeqJob.addFollowOn(seedCountsJob)
            seedCountsJob.addFollowOn(deBruijnGraphJob)

            seedProbabilitiesJob = Job.wrapJobFn(makeSeedProbabilities, options, sequenceID, seedCountsJob.rv(), graphID=deBruijnGraphJob.rv())
            seedCountsJob.addFollowOn(seedProbabilitiesJob)
            lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedProbabilitiesJob.rv())

            seedProbabilitiesJob.addFollowOn(lastzJob)

        else:
            lastzJob = Job.wrapJobFn(runLastz, options, trimSeqJob.rv(), trimSeqJob.rv())
            trimSeqJob.addFollowOn(lastzJob)

        printOutputJob = Job.wrapJobFn(printNumberOfHSPs, options, lastzJob.rv(1))
        lastzJob.addFollowOn(printOutputJob)

        toil.start(trimSeqJob)
def comparePrecisionRecall(options):
    with Toil(options) as toil:
        sequenceID = toil.importFile(toURL(options.sequence))
        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart + options.seqLength)

        seedProbabilitiesJob = Job.wrapJobFn(makeSeedProbabilities, options, sequenceID)
        sampledLastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedProbabilitiesJob.rv())

        maskedLastzJob = Job.wrapJobFn(runLastz, options, trimSeqJob.rv(), trimSeqJob.rv())

        precisionRecallJob = Job.wrapJobFn(alignmentPrecisionRecall, options, maskedLastzJob.rv(0), sampledLastzJob.rv(0))
        trimSeqJob.addFollowOn(seedProbabilitiesJob)
        seedProbabilitiesJob.addFollowOn(sampledLastzJob)
        sampledLastzJob.addFollowOn(maskedLastzJob)
        maskedLastzJob.addFollowOn(precisionRecallJob)
        toil.start(trimSeqJob)


def plotThreshold(options):
    with Toil(options) as toil:
        thresholdValues = [options.thresholdMin + options.thresholdStep*i for i in range(options.thresholdNValues)]
        sequenceID = toil.importFile(toURL(options.sequence))
        seedProbabilitiesJob = Job.wrapJobFn(makeSeedProbabilities, options, sequenceID)
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart+options.seqLength)
        dummyJob = Job.wrapJobFn(dummyJobFn)

        lastzJobs = []
        statsDicts = []
        for threshold in thresholdValues:
            lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedProbabilitiesJob.rv())
            statsDicts.append(lastzJob.rv(1))
            lastzJobs.append(lastzJob)

        printOutputJob = Job.wrapJobFn(printThresholdStats, options, thresholdValues, statsDicts)

        seedProbabilitiesJob.addFollowOn(trimSeqJob)
        trimSeqJob.addFollowOn(dummyJob)
        for lastzJob in lastzJobs:
            dummyJob.addChild(lastzJob)
        dummyJob.addFollowOn(printOutputJob)

        toil.start(seedProbabilitiesJob)

def main():
    parser = argparse.ArgumentParser()

    #Experiment types
    parser.add_argument("--plotThreshold", action="store_true")
    parser.add_argument("--plotScalability", action="store_true")
    parser.add_argument("--comparePrecisionRecall", action="store_true")

    #parameters
    parser.add_argument("--sampleSeeds", action="store_true")
    parser.add_argument("--clusterSeeds", action="store_true")
    parser.add_argument("--seedMultiplicity", action="store_true")
    parser.add_argument("--pathMultiplicity", action="store_true")

    parser.add_argument("--scoreThreshold", type=int, default=2)
    parser.add_argument("--sampleSeedThreshold", type=int, default=10)
    parser.add_argument("--sampleSeedStrategy", type=str, default="constantSampling")
    parser.add_argument("--region", type=str, default=None)
    parser.add_argument("--seqStart", type=int, default=0)
    parser.add_argument("--seqLength", type=int, default=10000)
    parser.add_argument("--maxSeqLength", type=int, default=50000)
    parser.add_argument("--k", type=int, default=19)
    parser.add_argument("--filterSize", type=int, default=10)
    parser.add_argument("--seedPattern", type=str, default="1110100110010101111")

    #options for plotting sample threshold
    parser.add_argument("--thresholdMin", type=int, default=2)
    parser.add_argument("--thresholdStep", type=int, default=20)
    parser.add_argument("--thresholdNValues", type=int, default=1)

    #input and output
    parser.add_argument("--sequence", type=str)
    parser.add_argument("--seedProbabilitiesFilename", type=str, default=None)
    parser.add_argument("--reuseProbabilitiesFile", action="store_true")
    parser.add_argument("--precisionRecallFile", type=str, default=None)
    parser.add_argument("--outputFile", type=str)

    Job.Runner.addToilOptions(parser)

    options = parser.parse_args()

    if options.plotThreshold:
        plotThreshold(options)
    elif options.plotScalability:
        plotScalability(options)
    elif options.comparePrecisionRecall:
        comparePrecisionRecall(options)
    else:
        runWorkflow(options)


if __name__ == "__main__":
    main()

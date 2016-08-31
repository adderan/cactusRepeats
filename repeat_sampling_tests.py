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

def makeSeedScoresTable(job, options, rawSeedCountsIDs):
    if options.reuseScoresFile:
        return job.fileStore.writeGlobalFile(options.seedScoresFile)
    job.fileStore.logToMaster("Making seed scores table...")
    rawSeedCountsFiles = [job.fileStore.readGlobalFile(rawSeedCountsID) for rawSeedCountsID in rawSeedCountsIDs]
    seedScoresTable = job.fileStore.getLocalTempFile()

    clusterSeedsOptions = ""
    if options.clusterSeeds:
        clusterSeedsOptions = "--clusterSeeds"

    os.system("cactus_blast_makeSeedScoreTable %s --scoreThreshold=%i --seedScoresFile=%s %s"
            % (clusterSeedsOptions, options.scoreThreshold, seedScoresTable, " ".join(rawSeedCountsFiles)))
    seedScoresID = job.fileStore.writeGlobalFile(seedScoresTable)
    if options.seedScoresFile and not options.reuseScoresFile:
        job.fileStore.exportFile(seedScoresID, toURL(options.seedScoresFile))
    return seedScoresID

def getChunks(job, options, sequenceIDs):
    job.fileStore.logToMaster("Making chunks....")
    sequences = [job.fileStore.readGlobalFile(sequenceID) for sequenceID in sequenceIDs]
    chunksDir = os.path.join(job.fileStore.getLocalTempDir(), "chunks")
    os.makedirs(chunksDir)
    chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                                          ("DEBUG",
                                                          options.chunkSize,
                                                          options.overlapSize,
                                                          chunksDir,
                                                          " ".join(sequences))).split("\n") if chunk != "" ]
    job.fileStore.logToMaster("Made %i chunks" % len(chunks))
    return [job.fileStore.writeGlobalFile(chunk) for chunk in chunks]

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


def runSampledLastz(job, options, seqID1, seqID2, seedScoresID, sampleSeedThreshold, sampleSeedStrategy):
    job.fileStore.logToMaster("Running lastz with seed sampling, threshold = %i, sampling strategy = %s ..." % (sampleSeedThreshold, sampleSeedStrategy))
    seq1 = job.fileStore.readGlobalFile(seqID1)
    seq2 = job.fileStore.readGlobalFile(seqID2)
    seedScores = job.fileStore.readGlobalFile(seedScoresID)
    alignmentsFile = job.fileStore.getLocalTempFile()
    statsFile = job.fileStore.getLocalTempFile()
    cmd = "cPecanLastz %s[multiple][unmask] %s[multiple][unmask] --notrivial --format=maf --seedCountsTable=%s --sampleSeedThreshold=%i --sampleSeedStrategy=%s --stats=%s > %s" % (seq1, seq2, seedScores, sampleSeedThreshold, sampleSeedStrategy, statsFile, alignmentsFile)
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
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        dummyJob = Job.wrapJobFn(dummyJobFn)
        trimSeqJobs = []
        lastzJobs = []
        statsList = []
        seqStep = options.seqLength/nPoints
        for i in range(1, nPoints):
            trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart + i*seqStep)
            if options.sampleSeeds:
                lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedScoreTableJob.rv(), options.sampleSeedThreshold, options.sampleSeedStrategy)
            else:
                lastzJob = Job.wrapJobFn(runLastz, options, trimSeqJob.rv(), trimSeqJob.rv())
            statsList.append(lastzJob.rv(1))
            lastzJobs.append(lastzJob)
            trimSeqJobs.append(trimSeqJob)

        printStatsJob = Job.wrapJobFn(printScalabilityStats, options, statsList)
        seedCountsJob.addFollowOn(seedScoreTableJob)
        seedScoreTableJob.addFollowOn(dummyJob)
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
            seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
            lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedScoreTableJob.rv(), sampleSeedThreshold=options.sampleSeedThreshold, sampleSeedStrategy=options.sampleSeedStrategy)

            trimSeqJob.addFollowOn(seedCountsJob)
            seedCountsJob.addFollowOn(seedScoreTableJob)
            seedScoreTableJob.addFollowOn(lastzJob)
        else:
            lastzJob = Job.wrapJobFn(runLastz, options, trimSeqJob.rv(), trimSeqJob.rv())
            trimSeqJob.addFollowOn(lastzJob)

        printOutputJob = Job.wrapJobFn(printNumberOfHSPs, options, lastzJob.rv(1))
        lastzJob.addFollowOn(printOutputJob)

        toil.start(trimSeqJob)
def comparePrecisionRecall(options):
    with Toil(options) as toil:
        sequenceID = toil.importFile(toURL(options.sequence))
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart + options.seqLength)

        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        sampledLastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedScoreTableJob.rv(), sampleSeedThreshold=options.sampleSeedThreshold, sampleSeedStrategy = options.sampleSeedStrategy)

        maskedLastzJob = Job.wrapJobFn(runLastz, options, trimSeqJob.rv(), trimSeqJob.rv())

        precisionRecallJob = Job.wrapJobFn(alignmentPrecisionRecall, options, maskedLastzJob.rv(0), sampledLastzJob.rv(0))
        trimSeqJob.addFollowOn(seedCountsJob)
        seedCountsJob.addFollowOn(seedScoreTableJob)
        seedScoreTableJob.addFollowOn(sampledLastzJob)
        sampledLastzJob.addFollowOn(maskedLastzJob)
        maskedLastzJob.addFollowOn(precisionRecallJob)
        toil.start(trimSeqJob)


def plotThreshold(options):
    with Toil(options) as toil:
        thresholdValues = [options.thresholdMin + options.thresholdStep*i for i in range(options.thresholdNValues)]
        sequenceID = toil.importFile(toURL(options.sequence))
        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStart+options.seqLength)
        dummyJob = Job.wrapJobFn(dummyJobFn)

        lastzJobs = []
        statsDicts = []
        for threshold in thresholdValues:
            lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedScoreTableJob.rv(), sampleSeedThreshold=threshold, sampleSeedStrategy="rejectAll")
            statsDicts.append(lastzJob.rv(1))
            lastzJobs.append(lastzJob)

        printOutputJob = Job.wrapJobFn(printThresholdStats, options, thresholdValues, statsDicts)

        seedCountsJob.addFollowOn(seedScoreTableJob)
        seedScoreTableJob.addFollowOn(trimSeqJob)
        trimSeqJob.addFollowOn(dummyJob)
        for lastzJob in lastzJobs:
            dummyJob.addChild(lastzJob)
        dummyJob.addFollowOn(printOutputJob)

        toil.start(seedCountsJob)

def main():
    parser = argparse.ArgumentParser()

    #Experiment types
    parser.add_argument("--plotThreshold", action="store_true")
    parser.add_argument("--plotScalability", action="store_true")
    parser.add_argument("--comparePrecisionRecall", action="store_true")

    #parameters
    parser.add_argument("--sampleSeeds", action="store_true")
    parser.add_argument("--clusterSeeds", action="store_true");

    parser.add_argument("--scoreThreshold", type=int, default=2)
    parser.add_argument("--sampleSeedThreshold", type=int, default=10)
    parser.add_argument("--sampleSeedStrategy", type=str, default="constantSampling")
    parser.add_argument("--region", type=str, default=None)
    parser.add_argument("--seqStart", type=int, default=0)
    parser.add_argument("--seqLength", type=int, default=10000)
    parser.add_argument("--maxSeqLength", type=int, default=50000)

    #options for plotting sample threshold
    parser.add_argument("--thresholdMin", type=int, default=2)
    parser.add_argument("--thresholdStep", type=int, default=20)
    parser.add_argument("--thresholdNValues", type=int, default=1)

    #input and output
    parser.add_argument("--sequence", type=str)
    parser.add_argument("--seedScoresFile", type=str, default=None)
    parser.add_argument("--reuseScoresFile", action="store_true")
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

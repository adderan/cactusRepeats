import argparse
import re
import os
import shutil

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
    return job.fileStore.writeGlobalFile(rawSeedCounts)

def makeSeedScoresTable(job, options, rawSeedCountsIDs):
    job.fileStore.logToMaster("Making seed scores table...")
    rawSeedCountsFiles = [job.fileStore.readGlobalFile(rawSeedCountsID) for rawSeedCountsID in rawSeedCountsIDs]
    seedScoresTable = job.fileStore.getLocalTempFile()
    os.system("cactus_blast_makeSeedScoreTable --countThreshold=%i --seedScoresFile=%s %s" 
            % (options.countThreshold, seedScoresTable, " ".join(rawSeedCountsFiles)))
    seedScoresID = job.fileStore.writeGlobalFile(seedScoresTable)
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

def trimSequence(job, options, sequenceID, seqStart, seqStop):
    if options.seqStart is None or options.seqStop is None:
        return sequenceID
    sequenceFile = job.fileStore.readGlobalFile(sequenceID)
    trimmedSequenceFile = job.fileStore.getLocalTempFile()
    for name, seq in fastaRead(sequenceFile):
        name = name + "|%i_%i" % (seqStart, seqStop)
        fastaWrite(trimmedSequenceFile, name, seq[seqStart:seqStop])
    seqID = job.fileStore.writeGlobalFile(trimmedSequenceFile)
    return seqID




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
    os.system("cPecanLastz %s[multiple] %s[multiple] --format=cigar --stats=%s > %s" % (seq1, seq2, statsFile, alignmentsFile))
    return job.fileStore.writeGlobalFile(statsFile)


def runSampledLastz(job, options, seqID1, seqID2, seedScoresID, sampleSeedThreshold, sampleSeedConstant):
    job.fileStore.logToMaster("Running lastz with seed sampling, threshold = %i, sampling constant = %i ..." % (sampleSeedThreshold, sampleSeedConstant))
    seq1 = job.fileStore.readGlobalFile(seqID1)
    seq2 = job.fileStore.readGlobalFile(seqID2)
    seedScores = job.fileStore.readGlobalFile(seedScoresID)
    alignmentsFile = job.fileStore.getLocalTempFile()
    statsFile = job.fileStore.getLocalTempFile()
    cmd = "cPecanLastz %s[multiple][unmask] %s[multiple][unmask] --notrivial --format=cigar --seedCountsTable=%s --sampleSeedThreshold=%i --sampleSeedConstant=%i --stats=%s > %s" % (seq1, seq2, seedScores, sampleSeedThreshold, sampleSeedConstant, statsFile, alignmentsFile)
    job.fileStore.logToMaster("Running lastz command: %s" % cmd)
    messages = system(cmd)
    stats = parseStats(statsFile)
    job.fileStore.exportFile(seq1, toURL("/home/alden/cactusRepeats/seq.txt"))
    job.fileStore.logToMaster("Number of HSPs = %i" % stats["HSPs"])
    job.fileStore.logToMaster("Seeds skipped = %i" % stats["Seeds skipped by sampling"])

    return stats


def printThresholdVsHSPs(job, options, thresholdValues, statsDicts):
    job.fileStore.logToMaster("Writing output file...")
    outputFile = job.fileStore.getLocalTempFile()
    with open(outputFile, 'w') as fh:
        for threshold, stats in zip(thresholdValues, statsDicts):
            fh.write("%i %i\n" % (threshold, stats["HSPs"]))
    outputID = job.fileStore.writeGlobalFile(outputFile)
    job.fileStore.exportFile(outputID, toURL(options.outputFile))
    
def printNumberOfHSPs(job, options, stats):
    job.fileStore.logToMaster("HSPs = %i" % stats["HSPs"])
def dummyJobFn(job):
    pass

def runWorkflow(options):
    with Toil(options) as toil:
        sequenceID = toil.importFile(toURL(options.sequence))
        fullSequenceID = toil.importFile(toURL(options.fullSequence))

        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, fullSequenceID)
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStop)

        lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedScoreTableJob.rv(), sampleSeedThreshold=options.sampleSeedThreshold, sampleSeedConstant=options.sampleSeedConstant)

        printOutputJob = Job.wrapJobFn(printNumberOfHSPs, options, lastzJob.rv())

        seedCountsJob.addFollowOn(seedScoreTableJob)
        seedScoreTableJob.addFollowOn(trimSeqJob)
        trimSeqJob.addFollowOn(lastzJob)
        lastzJob.addFollowOn(printOutputJob)
        toil.start(seedCountsJob)

    
def plotThreshold(options):
    with Toil(options) as toil:
        thresholdValues = range(2, 20)
        sequenceID = toil.importFile(toURL(options.sequence))
        if options.fullSequence:
            fullSequenceID = toil.importFile(toURL(options.fullSequence))
            seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, fullSequenceID)
        else:
            seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        trimSeqJob = Job.wrapJobFn(trimSequence, options, sequenceID, options.seqStart, options.seqStop)
        dummyJob = Job.wrapJobFn(dummyJobFn)

        lastzJobs = []
        statsDicts = []
        for threshold in thresholdValues:
            lastzJob = Job.wrapJobFn(runSampledLastz, options, trimSeqJob.rv(), trimSeqJob.rv(), seedScoreTableJob.rv(), sampleSeedThreshold=threshold, sampleSeedConstant=0)
            statsDicts.append(lastzJob.rv())
            lastzJobs.append(lastzJob)

        printOutputJob = Job.wrapJobFn(printThresholdVsHSPs, options, thresholdValues, statsDicts)

        seedCountsJob.addFollowOn(seedScoreTableJob)
        seedScoreTableJob.addFollowOn(trimSeqJob)
        trimSeqJob.addFollowOn(dummyJob)
        for lastzJob in lastzJobs:
            dummyJob.addChild(lastzJob)
            lastzJob.addFollowOn(printOutputJob)

        toil.start(seedCountsJob)

def main():
    parser = argparse.ArgumentParser()

    #Experiment types
    parser.add_argument("--plotThreshold", action="store_true")
    parser.add_argument("--runSampled", action="store_true")
    parser.add_argument("--runUnsampled", action="store_true")

    #parameters
    parser.add_argument("--chunkSize", type=int, default=5000)
    parser.add_argument("--overlapSize", type=int, default=1000)
    parser.add_argument("--countThreshold", type=int, default=2)
    parser.add_argument("--sampleSeedThreshold", type=int, default=2)
    parser.add_argument("--sampleSeedConstant", type=int, default=2)
    parser.add_argument("--seqStart", type=int, default=None)
    parser.add_argument("--seqStop", type=int, default=None)

    #input and output
    parser.add_argument("--sequence", type=str)
    parser.add_argument("--fullSequence", type=str)
    parser.add_argument("--outputFile", type=str)

    Job.Runner.addToilOptions(parser)

    options = parser.parse_args()

    if options.plotThreshold:
        plotThreshold(options)
    elif options.runSampled:
        runWorkflow(options)
    elif options.runUnsampled:
        runWorkflow(options)


if __name__ == "__main__":
    main()

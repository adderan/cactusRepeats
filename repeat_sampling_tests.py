import argparse
import re
import os

from toil.common import Toil
from toil.job import Job
from sonLib.bioio import popenCatch

def toURL(path):
    if path.startswith("http://"):
        return path
    else:
        return "file://" + os.path.abspath(path)

def makeRawSeedCounts(job, options, sequenceID):
    job.fileStore.logToMaster("Getting seed counts from lastz...")
    sequence = job.fileStore.readGlobalFile(sequenceID)
    rawSeedCounts = job.fileStore.getLocalTempFile()
    os.system("cactus_lastz --tableonly=count %s > %s" % (sequence, rawSeedCounts))
    return job.fileStore.writeGlobalFile(rawSeedCounts)

def makeSeedScoresTable(job, options, rawSeedCountsIDs):
    job.fileStore.logToMaster("Making seed scores table...")
    rawSeedCountsFiles = [job.fileStore.readGlobalFile(rawSeedCountsID) for rawSeedCountsID in rawSeedCountsIDs]
    seedScoresTable = job.fileStore.getLocalTempFile()
    os.system("cactus_blast_makeSeedScoresTable --countThreshold=%i --seedScoresFile=%s %s" 
            % (options.countThreshold, seedScoresTable, " ".join(rawSeedCountsFiles)))
    return job.fileStore.writeGlobalFile(seedScoresTable)

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
    return [job.fileStore.writeGlobalFile(chunk) for chunk in chunks]

def parseStats(job, options, statsID):
    """Parse the lastz stats file.
    """
    statsFile = job.fileStore.readGlobalFile(statsID)
    stat = re.complie("(\S+): (\d+)")

    stats = dict()
    for line in open(stats):
        line = line.strip()
        m = stat.match(line)
        if m:
            stats[m.group(0)] = m.group(1)
    return stats


def runSampledLastz(job, options, seqID1, seqID2, seedScoresID):
    job.fileStore.logToMaster("Running lastz with seed sampling...")
    seq1 = job.fileStore.readGlobalFile(seqID1)
    seq2 = job.fileStore.readGlobalFile(seqID2)
    seedScores = job.fileStore.readGlobalFile(seedScoresID)
    alignmentsFile = job.fileStore.getLocalTempFile()
    statsFile = job.fileStore.getLocalTempFile()
    os.system("cactus_lastz %s[unmask][multiple] %s[unmask][multiple] --seedScoresFile=%s --stats=%s --sampleSeedThreshold=%i > %s" % (seq1, seq2, seedScores, options.sampleSeedThreshold, statsFile, alignmentsFile))
    return job.fileStore.writeGlobalFile(statsFile)


def printThresholdVsHSPs(options, thresholdValues, statsDicts):
    outputFile = job.fileStore.getLocalTempFile()
    with open(outputFile, 'w') as fh:
        for threshold, stats in zip(thresholdValues, statsDicts):
            fh.write("%i %i\n" % (threshold, stats["HSPs"]))
    return job.fileStore.writeGlobalFile(outputFile)
    
    
def plotThreshold(options, sequence):
    with Toil(options) as toil:
        sequenceID = toil.importFile(toURL(sequence))
        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoresTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        chunksJob = Job.wrapJobFn(getChunks, options, [sequenceID])
        lastzJob = Job.wrapJobFn(runSampledLastz, chunksJob.rv(0), chunksJob.rv(0), seedScoreTableJob.rv())
        parseStatsJob = Job.wrapJobFn(parseStats, options, statsID=lastzJob.rv())
        printOutputJob = Job.wrapJobFn(printThresholdVsHSPs, options, [2], [parseStatsJob.rv()])

        seedCountsJob.addFollowOn(seedScoreTableJob)
        seedScoreTableJob.addFollowOn(chunksJob)
        chunksJob.addFollowOn(lastzJob)
        lastzJob.addFollowOn(parseStatsJob)
        parseStatsJob.addFollowOn(printOutputJob)

        toil.start(seedCountsJob)

        toil.exportFile(printOutputJob.rv(), toURL(options.outputFile))


def main():
    parser = argparse.ArgumentParser()

    #Experiment types
    parser.add_argument("--thresholdVsHSPs", action="store_true")

    #parameters
    parser.add_argument("--chunkSize", type=int, default=5000)
    parser.add_argument("--overlapSize", type=int, default=1000)
    parser.add_argument("--countThreshold", type=int, default=2)
    parser.add_argument("--sampleSeedThreshold", type=int, default=2)

    #input and output
    parser.add_argument("--sequences", type=str)
    parser.add_argument("--outputFile", type=str)

    Job.Runner.addToilOptions(parser)

    options = parser.parse_args()
    sequences = options.sequences.split(',')

    if options.thresholdVsHSPs:
        plotThreshold(options, sequences[0])


if __name__ == "__main__":
    main()

import argparse
import system
import re
import os

from toil import Toil

def makeRawSeedCounts(job, options, sequenceID):
    sequence = job.fileStore.readGlobalFile(sequenceID)
    rawSeedCounts = job.fileStore.getLocalTempFile()
    system("cactus_lastz --tableonly=count %s > %s" % (sequence, rawSeedCounts))
    return job.fileStore.writeGlobalFile(rawSeedCounts)

def makeSeedScoresTable(job, options, rawSeedCountsIDs):
    rawSeedCountsFiles = [job.fileStore.readGlobalFile(rawSeedCountsID) for rawseedCountsID in rawSeedCountsIDs]
    seedScoresTable = job.fileStore.getLocalTempFile()
    system("cactus_blast_makeSeedScoresTable --countThreshold=%i --seedScoresFile=%s %s" 
            % (options.countThreshold, seedScoresTable, " ".join(rawSeedCountsFiles)))
    return job.fileStore.writeGlobalFile(seedScoresTable)

def getChunks(job, options, sequenceIDs):
    sequences = [job.fileStore.readGlobalFile(sequenceID) for sequenceID in sequenceIDs]
    chunksDir = os.path.join(job.fileStore.getLocalTempDir(), "chunks")
    os.makedirs(chunksDir)
    chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                                          (getLogLevelString(),
                                                          options.chunkSize,
                                                          options.overlapSize,
                                                          chunksDir,
                                                          " ".join(sequences))).split("\n") if chunk != "" ]
    return [job.fileStore.writeGlobalFile(chunk) for chunk in chunks]

def parseStats(job, statsID):
    statsFile = job.fileStore.readGlobalFile(statsID)
    stat = re.complie("(\S+): (\d+)")

    stats = dict()
    for line in open(stats):
        line = line.strip()
        m = stat.match(line)
        if m:
            stats[m.group(0)] = m.group(1)
    return stats


def runSampledLastz(job, seqID1, seqID2, seedScoresID, scoreThreshold=2):
    seq1 = job.fileStore.readGlobalFile(seqID1)
    seq2 = job.fileStore.readGlobalFile(seqID2)
    seedScores = job.fileStore.readGlobalFile(seedScoresID)
    alignmentsFile = job.fileStore.getLocalTempFile()
    statsFile = job.fileStore.getLocalTempFile()
    system("cactus_lastz %s[unmask][multiple] %s[unmask][multiple] --seedScoresFile=%s --stats=%s --sampleSeedThreshold=%i > %s" % (seq1, seq2, seedScores, scoreThreshold, statsFile, alignmentsFile))
    return job.fileStore.writeGlobalFile(statsFile)

    
def plotThreshold(options, sequence):
    with toil(options) as Toil:
        sequenceID = toil.importFile("file://" + sequence)
        seedCountsJob = Job.wrapJobFn(makeRawSeedCounts, options, sequenceID)
        seedScoreTableJob = Job.wrapJobFn(makeSeedScoreTable, options, rawSeedCountsIDs = [seedCountsJob.rv()])
        chunksJob = Job.wrapJobFn(getChunks, sequenceID)

        lastzJob = Job.wrapJobFn(runSampledLastz
        seedCountsJob.addChild(seedScoreTableJob)


def main():
    parser = argparse.ArgumentParser()

    #Experiment types
    parser.add_argument("--thresholdVsHsps", action="store_true")

    #parameters
    parser.add_argument("--chunkSize", type=int, default=5000)
    parser.add_argument("--overlapSize", type=int, default=1000)
    parser.add_argument("--sampleSeedThreshold", type=int, default=2)

    #input and output
    parser.add_argument("--sequences", type=str)
    parser.add_argument("--outputFile", type=str)

    Job.Runner.addToilOptions(parser)

    options = parser.parse_args()
    sequences = options.sequences.split(',')

    if args.thresholdVsHsps:
        plotThreshold(options, 


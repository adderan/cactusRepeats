import matplotlib
import matplotlib.pyplot as plt
import sys
import numpy
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outputGraph", type=str)
    args = parser.parse_args()
    nNodes = []
    nEdges = []
    for line in sys.stdin:
        if len(line.split()) != 2:
            continue
        nNodes.append(line.split()[0])
        nEdges.append(line.split()[1])
    fig = plt.figure()
    plt.scatter(numpy.array(nNodes), numpy.array(nEdges))
    fig.savefig(args.outputGraph)


if __name__ == "__main__":
    main()

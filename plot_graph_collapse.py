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
        nNodes.append(float(line.split()[0]))
        nEdges.append(float(line.split()[1]))
    fig = plt.figure()
    nNodes = numpy.array(nNodes)
    nEdges = numpy.array(nEdges)
    plt.scatter(nNodes, nEdges)
    plt.xlabel("Number of nodes in graph")
    plt.title("Number of random edges to reach 1 connected component")
    a, b = numpy.polyfit(nNodes, nEdges, 1)
    offset = 3000
    b = b + offset
    print "a = %f, b = %f" % (a, b)
    plt.plot(nNodes, nNodes*a + b)
    fig.savefig(args.outputGraph)


if __name__ == "__main__":
    main()

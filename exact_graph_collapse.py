import matplotlib
import matplotlib.pyplot as plt
import argparse
import scipy
import scipy.misc
import math
import numpy


cache = dict()

def P_n(n, p):
    if n == 1:
        return 1.0
    if p == 0.0:
        return 0.0

    term1 = 0
    for k in range(1, n-1):
        P_n_k_p = None
        if k in cache.keys():
            P_n_k_p = cache[k]
        else:
            P_n_k_p = P_n(k, p)
            cache[k] = P_n_k_p
            print "Caching P_n(%i, %f) = %f" % (k, p, P_n_k_p)
        if not P_n_k_p == P_n(k, p):
            print P_n_k_p
            print P_n(k, p)
            assert False
        term1 += scipy.misc.comb(n-1, k-1) * P_n_k_p * math.pow(1.0 - p, k*(n-k))

    return 1.0 - term1

def numberEdgesToCollapse(n, threshold):
    if n == 1:
        return 0
    e = 0
    P = 0
    nPossibleEdges = float(scipy.misc.comb(n,2))
    while P < threshold:
        e = e + 1
        cache = None
        cache = dict()
        P = P_n(n, e/nPossibleEdges)

    return e

def numberOfEdgesToCollapseAtInfinity(n, threshold):
    if n == 1:
        return 0
    q = math.pow((1.0-threshold)/n, 1/(n-1.0))
    p = 1.0 - q
    nEdges = p * scipy.misc.comb(n, 2)
    return nEdges
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float)
    parser.add_argument("--nMax", type=int)
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--outputFile", type=str)

    options = parser.parse_args()


    nValues = []
    eValues = []

    for n in range(1, options.nMax):
        e = numberOfEdgesToCollapseAtInfinity(n, options.threshold)
        eValues.append(e)
        nValues.append(n)
        print "n = %i, e = %i" % (n, e)
    naive = [scipy.misc.comb(i, 2) for i in nValues]


    fig = plt.figure()
    plt.plot(numpy.array(nValues), numpy.array(eValues), linewidth=4, color="blue", label="Number of edges to collapse graph")
    plt.plot(numpy.array(nValues), numpy.array(naive), linewidth=4, color="red", label="All edges")
    #plt.title("How many edges are necessary?")
    plt.xlabel("Number of nodes")
    plt.ylabel("Number of edges")
    plt.legend()
    fig.savefig(options.outputFile)


if __name__ == "__main__":
    main()

#include "sonLib.h"
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int addEdgesUntilCollapse(int n) {
    if (n == 0) return 0;
    srand(time(NULL));
    stConnectivity *graph = stConnectivity_construct();
    for (long i = 0; i < n; i++) {
        stConnectivity_addNode(graph, (void*)i);
    }
    int nEdges = 0;
    while (true) {
        long i = 0;
        long j = 0;

        while (i == j || stConnectivity_hasEdge(graph, (void*)i, (void*)j)) {
            i = rand() % n;
            j = rand() % n;
        }
        nEdges++;
        stConnectivity_addEdge(graph, (void*)i, (void*)j);
        if (stConnectivity_getNComponents(graph) == 1) {
            return nEdges;
        }
    }
    stConnectivity_destruct(graph);
}


int main(int argc, char **argv) {
    int maxNodes;
    int key;
    while (1) {
        static struct option long_options[] = { { "maxNodes", required_argument, 0, 'a' },
            { 0, 0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                sscanf(optarg, "%i", &maxNodes);
                break;
        }

        for (int i = 0; i < maxNodes; i += 10) {
            printf("%i %i\n", i, addEdgesUntilCollapse(i));
        }
            
    }
}


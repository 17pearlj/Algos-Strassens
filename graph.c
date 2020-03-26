#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// functions to construct a random graph and to print it
void graph(double p, int v_count, int **adj);
void print_graph(int v_count, int **adj);

int main(int argc, char** argv) {
    // ensure correct usage
    if (argc != 3){
        printf("Usage: ./graph p v_count");
        return 0;
    }

    // Seed RNG with current time
    srandom(time(NULL));

    // calloc space for our graph
    int v_count = atoi(argv[2]);
    int *adj[v_count];
    for (int i = 0; i < v_count; i++){
        adj[i] = (int *)calloc(v_count, sizeof(int));
    }

    // add random edges to our adjacency matrix and print the graph
    graph(atof(argv[1]), v_count, adj);
    print_graph(v_count, adj);

    // free the space taken up by our graph
    for (int i = 0; i < v_count; i++){
        free(adj[i]);
    }

    return 0;

}

// 2D graph, 1024 vertices, each edge has a probability p of being generated
void graph(double p, int v_count, int **adj) {
    for (int i = 0; i < v_count; i++) {
        for (int j = i; j < v_count; j++) {
            // ignore reflexive edges
            if (j == i) {
                continue;
            }
            // add each edge with probability p
            else if (((double)random())/((double)(RAND_MAX)) <= p) {
                adj[i][j] = 1;
                adj[j][i] = 1;
            }
            else {
                continue;
            }
        }
    }
}

void print_graph(int v_count, int **adj) {
    for (int i = 0; i < v_count; i++) {
        for (int j = 0; j < v_count; j++) {
            printf("%i ", adj[i][j]);
        }
        printf("\n");
    }
}

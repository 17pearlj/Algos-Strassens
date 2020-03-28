#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "matrix.c"

int** big_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);
void diagonals(int **three, int dim);
void graph(double p, int v_count, int **adj);
double graph_triangles(int v_count, int **adj);
int **standard_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim);
int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);
void program(int mode, int pad_dim, int dim, char* fname);

int main(int argc, char** argv) {
   // Ensure correct usage
   if (argc != 4) {
         printf("Usage: ./strassen 0 dimension inputfile \n");
         return 1;
   }
   int dim = atoi(argv[2]);
   int pad_dim = pow(2, ceil(log(dim)/ log(2)));
   program(atoi(argv[1]), pad_dim, dim, argv[3]);
   return 0;

}

void program(int mode, int pad_dim, int dim, char* fname) {
   FILE *f = fopen(fname, "r");
   if (f == 0)
   {
      //fopen returns 0, the NULL pointer, on failure
      perror("Cannot open input file\n");
      exit(-1);
	}
   int **one = (int **)calloc(pad_dim, sizeof(int *));
   int **two = (int **)calloc(pad_dim, sizeof(int *));
   int **three = (int **)calloc(pad_dim, sizeof(int *));
   for (int i = 0; i < pad_dim; i++) {
      one[i] = (int *)calloc(pad_dim, sizeof(int));
      two[i] = (int *)calloc(pad_dim, sizeof(int));
      three[i] = (int *)calloc(pad_dim, sizeof(int));
   }

   makeTwoMatrices(f, dim, one, two);
   clock_t start, end;

   start = clock();
   if (mode == 0) {
      big_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim, pad_dim);
   }
   else if (mode == 1) {
      printf("standard mult ");
      standard_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim);
   } else {
       // Seed RNG with current time
       srandom(time(NULL));

       // calloc space for our graph
       int v_count = 1024;
       int *adj[v_count];
       int *adj2[v_count];
       int *adj3[v_count];
       for (int i = 0; i < v_count; i++) {
           adj[i] = (int *)calloc(v_count, sizeof(int));
           adj2[i] = (int *)calloc(v_count, sizeof(int));
           adj3[i] = (int *)calloc(v_count, sizeof(int));
       }

       // add random edges to our adjacency matrix and print the graph
       graph(0.01, v_count, adj);

       big_mult(adj, adj, adj2, 0, 0, 0, 0, 0, 0, v_count, v_count);
       big_mult(adj2, adj, adj3, 0, 0, 0, 0, 0, 0, v_count, v_count);

       printf("num triangles: %f\n", graph_triangles(v_count, adj3));

       // free the space taken up by our graph
       for (int i = 0; i < v_count; i++) {
           free(adj[i]);
           free(adj2[i]);
           free(adj3[i]);
       }
   }
   end = clock();
   double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
   if (mode == 0) {
      diagonals(three, dim);
   }
   for (int i = 0; i < pad_dim; i ++){
      free(one[i]);
      free(two[i]);
      free(three[i]);
   }

   free(one);
   free(two);
   free(three);
}

void diagonals(int **three, int dim) {
   for (int i = 0; i < dim; i++){
      printf("%d \n", three[i][i]);
   }
   printf("\n");
}

/*big multiplication algorithm
RETURNS: a list of diagonal entries
tasks:
1) turn file into matrix form (if necessary? is this necessary?)
2) based on n, either call strassens or call standard mult
3) collect diagonals
*/
int** big_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim) {
   // check dimensions, if under crossover point, call standard, else call strassens
   if (r1 >= real_dim || c1 >= real_dim) {
      return one;
   } else if (r2 >= real_dim || c2 >= real_dim) {
      return two;
   } else if (dim < 84){
      return standard_mult(one, two, three, r1, c1, r2, c2, r3, c3, dim);
   } else {
      return strassen(one, two, three, r1, c1, r2, c2, r3, c3, real_dim, dim);
   }
}

/*standard matrix multiplication
RETURNS: int**
tasks:
1) perform standard mult
*/
int** standard_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim) {

    // we need to make copies of the matrices so we can line up the indices
    int **temp1 = (int **)calloc(dim, sizeof(int *));
    int **temp2 = (int **)calloc(dim, sizeof(int *));
    int **temp3 = (int **)calloc(dim, sizeof(int *));
    for (int i = 0; i < dim; i++) {
       temp1[i] = (int *)calloc(dim, sizeof(int));
       temp2[i] = (int *)calloc(dim, sizeof(int));
       temp3[i] = (int *)calloc(dim, sizeof(int));
    }

    // make 0-indexed copies of matrices one and two
    for (int i = r1; i < r1 + dim; i++) {
        for (int j = c1; j < c1 + dim; j++) {
            temp1[i - r1][j - c1] = one[i][j];
        }
    }
    for (int i = r2; i < r2 + dim; i++) {
        for (int j = c2; j < c2 + dim; j++) {
            temp2[i - r2][j - c2] = two[i][j];
        }
    }

    // multiply using the 0-indexed copies
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                temp3[i][j] += (temp1[i][k] * temp2[k][j]);
            }
        }
    }

    // copy the product into matrix 3 at the appropriate indices
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            three[i + r3][j + c3] = temp3[i][j];
        }
    }

    // free our temporary 0-indexed copies
    for (int i = 0; i < dim; i ++) {
       free(temp1[i]);
       free(temp2[i]);
       free(temp3[i]);
    }
    free(temp1);
    free(temp2);
    free(temp3);
    return three;
}




/*strassens algorithm
RETURNS: int**
tasks:
1) add necessary padding so that it's of a dimension that's a power of two
2) perform strassens -> recursive call to big_mult in it
3) in strassen, add each term to the correct part of the third and final matrix
*/

int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim) {
   int new_d = dim/2;
   if (new_d == 1) {
      three[r3][c3 + 1] +=  (one[r1][c1])*(two[r2][c2 + 1] - two[r2 + 1][c2 + 1]);
      three[r3 + 1][c3 + 1] += (one[r1][c1])*(two[r2][c2 + 1] - two[r2 + 1][c2 + 1]);
      three[r3][c3 + 1] += (one[r1][c1] + one[r1][c1 + 1])*(two[r2 + 1][c2 + 1]);
      three[r3][c3] -= (one[r1][c1] + one[r1][c1 + 1])*(two[r2 + 1][c2 + 1]);
      three[r3 + 1][c3] += (one[r1 + 1][c1] + one[r1 + 1][c1 + 1])*(two[r2][c2]);
      three[r3 + 1][c3 + 1] -= (one[r1 + 1][c1] + one[r1 + 1][c1 + 1])*(two[r2][c2]);
      three[r3][c3] += (one[r1 + 1][c1 + 1])*(two[r2 + 1][c2] - two[r2][c2]);
      three[r3 + 1][c3] +=(one[r1 + 1][c1 + 1])*(two[r2 + 1][c2] - two[r2][c2]);
      three[r3][c3] += (one[r1][c1] + one[r1 + 1][c1 + 1])*(two[r2][c2] * two[r2 + 1][c2 + 1]);
      three[r3 + 1][c3 + 1] += (one[r1][c1] + one[r1 + 1][c1 + 1])*(two[r2][c2] * two[r2 + 1][c2 + 1]);
      three[r3][c3] += (one[r1][c1 + 1] - one[r1 + 1][c1 + 1])*(two[r2 + 1][c2] + two[r2 + 1][c2 + 1]);
      three[r3 + 1][c3 + 1] -= (one[r1][c1] - one[r1 + 1][c1])*(two[r2][c2] + two[r2][c2 + 1]);
      return three;
   } else {
      int ***sums = (int***)calloc(2, sizeof(int**));
      for (int i = 0; i < 2; i++) {
         sums[i] = (int **)calloc(new_d, sizeof(int*));
         for (int j = 0; j < new_d; j++) {
            sums[i][j] = (int *)calloc(new_d, sizeof(int));
         }
      }


      int **p = (int**)calloc(new_d, sizeof(int*));
      for (int j = 0; j < new_d; j++) {
         p[j] = (int *)calloc(new_d, sizeof(int));
      }
      for (int i = 0; i < 10; i++) {
         switch(i) {
            case 0: // f- h
               matrixDestinationAddition(two, r2, c2 + new_d, two, r2 + new_d, c2 + new_d, sums[0], new_d, -1);
               big_mult(one, sums[0], p, r1, c1, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 , c3 + new_d, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, 1);
               break;

            case 1: // a + b
               for (int j = 0; j < new_d; j++) {
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(one, r1, c1, one, r1, c1 + new_d, sums[0], new_d, 1);
               big_mult(sums[0], two, p, 0, 0, r2 + new_d, c2 + new_d, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 , c3 + new_d, new_d, 1);
               matrixCorAddition(p, three, r3, c3, new_d, -1);
               break;

            case 2: //c + d
               matrixDestinationAddition(one, r1 + new_d, c1, one, r1 + new_d, c1 + new_d, sums[0], new_d, 1);
               for (int j = 0; j < new_d; j++) {
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               big_mult(sums[0], two, p, 0, 0, r2, c2, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 + new_d, c3, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, -1);
               break;

            case 3: //g - e
               for (int j = 0; j < new_d; j++) {
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(two, r2 + new_d, c2, two, r2, c2, sums[0], new_d, -1);
               big_mult(one, sums[0], p, r1 + new_d, c1 + new_d, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3, c3, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3, new_d, 1);
               break;

            case 4: // a + d
               matrixDestinationAddition(one, r1, c1, one, r1 + new_d, c1 + new_d, sums[0], new_d, 1);
               break;

            case 5: //e + h
               for (int j = 0; j < new_d; j++) {
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(two, r2, c2, two, r2 + new_d, c2 + new_d, sums[1], new_d, 1);
               big_mult(sums[0], sums[1], p, 0, 0, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3, c3, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, 1);
               break;

            case 6: //b - d
               matrixDestinationAddition(one, r1, c1 + new_d, one, r1 + new_d, c1 + new_d, sums[0], new_d, -1);
               break;

            case 7: // g + h
               for (int j = 0; j < new_d; j++) {
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(two, r2 + new_d, c2, two, r2 + new_d, c2 + new_d, sums[1], new_d, 1);
               big_mult(sums[0], sums[1], p, 0, 0, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3, c3, new_d, 1);
               break;

            case 8: // a - c
               matrixDestinationAddition(one, r1, c1, one, r1 + new_d, c1, sums[0], new_d, -1);
               break;

            default: //e + f
               for (int j = 0; j < new_d; j++) {
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(two, r2, c2, two, r2, c2 + new_d, sums[1], new_d, 1);
               big_mult(sums[0], sums[1], p, 0, 0, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, -1);
               break;
         }
      }

      for (int i = 0; i < 2; i++) {
         for (int j = 0; j < new_d; j++) {
               free(sums[i][j]);
         }
         free(sums[i]);
      }
      free(sums);

      for (int j = 0; j < new_d; j++) {
         free(p[j]);
      }
      free(p);

      return three;
   }
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

double graph_triangles(int v_count, int **adj) {
    double sum = 0;
    for (int i = 0; i < v_count; i++) {
        sum += adj[i][i];
    }
    return (sum / 6.);
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int** big_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);
void diagonals(int **three, int dim);
void exec_fun(int mode, int pad_dim, int dim, char* fname);
void graph(double p, int v_count, int **adj);
double graph_triangles(int v_count, int **adj);
void matrixCorAddition(int **one, int **two, int r2, int c2, int dim, int positive);
int makeTwoMatrices(FILE *f, int dim, int** one, int** two);
void matrixDestinationAddition(int **one, int r1, int c1, int **two, int r2, int c2, int**three, int dim, int positive);
int **standard_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim);
int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);

int crossover = 73;

int main(int argc, char** argv) {
   // ensure correct usage
   if (argc != 4) {
         printf("Usage: ./strassen 0 dimension inputfile \n");
         return 1;
   }
   // execute the program, padding the matrix dimension as needed
   int dim = atoi(argv[2]);
   int pad_dim = pow(2, ceil(log(dim)/ log(2)));
   exec_fun(atoi(argv[1]), pad_dim, dim, argv[3]);
   return 0;
}

void exec_fun(int mode, int pad_dim, int dim, char* fname) {
   // open file containing matrix
   FILE *f = fopen(fname, "r");
   if (f == 0)
   {
      // fopen returns 0, the NULL pointer, on failure
      perror("Cannot open input file\n");
      fclose(f);
      exit(-1);
 	 }
   // allocate space for the two matrices we're multiplying and a third for the product
   int **one = (int **) calloc(pad_dim, sizeof(int *));
   int **two = (int **) calloc(pad_dim, sizeof(int *));
   int **three = (int **) calloc(pad_dim, sizeof(int *));
   for (int i = 0; i < pad_dim; i++) {
      one[i] = (int *) calloc(pad_dim, sizeof(int));
      two[i] = (int *) calloc(pad_dim, sizeof(int));
      three[i] = (int *) calloc(pad_dim, sizeof(int));
   }

   makeTwoMatrices(f, dim, one, two);

   // start timing
   clock_t start, end;
   start = clock();

   // mode 0 is default, others were used for development
   if (mode == 0) {
      big_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim, pad_dim);
   }
   else if (mode == 1) {
      standard_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim);
   }
   else {
       // seed RNG with current time
       srandom(time(NULL));
       double p_val[5] = {0.01, 0.02, 0.03, 0.04, 0.05};

       // run five trials for each p_val
       for (int a = 0; a < 5; a++) {
           printf("5 trials, p_val: %f \n", p_val[a]);
           for (int trial = 0; trial < 5; trial++) {
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

               // add random edges to adjacency matrix
               graph(p_val[a], v_count, adj);

               // cube the adjacency matrix to calculate the number of triangles in our graph
               big_mult(adj, adj, adj2, 0, 0, 0, 0, 0, 0, v_count, v_count);
               big_mult(adj2, adj, adj3, 0, 0, 0, 0, 0, 0, v_count, v_count);
               printf("num triangles: %f\n", graph_triangles(v_count, adj3));

               // free the space taken up by our matrices
               for (int i = 0; i < v_count; i++) {
                   free(adj[i]);
                   free(adj2[i]);
                   free(adj3[i]);
               }
           }
           printf("\n");
       }
   }
   end = clock();
   double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

   // print diagonals of our product matrix
   if (mode != 2) {
      diagonals(three, dim);
   }

   // free used memory
   for (int i = 0; i < pad_dim; i ++) {
      free(one[i]);
      free(two[i]);
      free(three[i]);
   }

   free(one);
   free(two);
   free(three);

   fclose(f);
}

// print the diagonals in matrix three
void diagonals(int **three, int dim) {
   for (int i = 0; i < dim; i++) {
      printf("%d \n", three[i][i]);
   }
}

/*big multiplication algorithm
RETURNS: a list of diagonal entries
tasks:
1) based on n, either call strassens or call standard mult
2) collect diagonals
*/
int** big_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim) {
   // check dimensions, if under crossover point, call standard, else call strassen
   if (r1 >= real_dim || c1 >= real_dim) {
      return one;
   } else if (r2 >= real_dim || c2 >= real_dim) {
      return two;
   } else if (dim < crossover) {
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

    // allocate space for 0-indexed copies of the matrices
    int **temp1 = (int **) calloc(dim, sizeof(int *));
    int **temp2 = (int **) calloc(dim, sizeof(int *));
    int **temp3 = (int **) calloc(dim, sizeof(int *));
    for (int i = 0; i < dim; i++) {
       temp1[i] = (int *) calloc(dim, sizeof(int));
       temp2[i] = (int *) calloc(dim, sizeof(int));
       temp3[i] = (int *) calloc(dim, sizeof(int));
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

    // multiply the 0-indexed copies of one and two
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                // store the product in 0-indexed matrix three
                temp3[i][j] += (temp1[i][k] * temp2[k][j]);
            }
        }
    }

    // copy the 0-indexed product into matrix three at the appropriate indices
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            three[i + r3][j + c3] = temp3[i][j];
        }
    }

    // free used memory
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


// generate 2D graph, 1024 vertices, each edge has a probability p of being generated
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

// calculate the number of triangles in graph represented by adj matrix
double graph_triangles(int v_count, int **adj) {
    // return the sum of the diagonal divided by 6
    double sum = 0;
    for (int i = 0; i < v_count; i++) {
        sum += adj[i][i];
    }
    return (sum / 6.);
}

int makeTwoMatrices(FILE *f, int dim, int** one, int** two) {
		char ch, buffer[11];
    ch = '\0';
    memset(buffer, '\0', sizeof(buffer));
		int i = 0, ar = 0, row = 0, col = 0, counter = 0;
		// while both arrays not read
		while (counter < (2*(pow(dim, 2.)))) {
				// reads char
				ch = fgetc(f);
			  // if EOF is encountered, error
				if(ch == EOF) {
						perror("File is too short");
						exit(-1);
						break;
				}
				else if(ch == '\n') {
						switch (ar) {
								case 0:
					      		one[row][col] = atoi(buffer);
										break;
								default: two[row][col] = atoi(buffer);
						}
						col++;
						if (col == dim) {
								col = 0;
								row++;
								if (row == dim) {
										ar++;
										row = 0;
								}
						}
			      bzero(buffer, 10);
						i = 0;
			      counter++;
						continue;
				}
				else {
						buffer[i] = ch;
						i++;
				}
		}
		return 0;
}


void matrixCorAddition(int **one, int **two, int r2, int c2, int dim, int positive) {
    for(int r = 0; r < dim; r++) {
				for (int c = 0; c < dim; c++) {
		            two[r2 + r][c2 + c] += (one[r][c]*positive);
				}
		}
}

void matrixDestinationAddition(int **one, int r1, int c1, int **two, int r2, int c2, int**three, int dim, int positive) {
    for(int r = 0; r < dim; r++) {
				for (int c = 0; c < dim; c++) {
		    		three[r][c] = one[r + r1][c + c1] + (two[r + r2][c + c2]*positive);
				}
		}
}

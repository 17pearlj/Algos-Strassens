#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.c"
#include <time.h>

int **standard_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim);
int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);
void diagonals(int **three, int dim);
int **big_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);
int **strassenOpt(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim);

int main(int argc, char** argv) {
   // Ensure correct usage
   if (argc != 4) {
         printf("Usage: ./strassen 0 dimension inputfile \n");
         return 1;
   }
   FILE *f = fopen(argv[3], "r");
   if (f == 0)
   {
      //fopen returns 0, the NULL pointer, on failure
      perror("Cannot open input file\n");
      exit(-1);
	}
   int dim = atoi(argv[2]);
   int pad_dim = pow(2, ceil(log(dim)/ log(2)));
   char ch, buffer[10];
   int **one = (int **)calloc(pad_dim, sizeof(int *));
   int **two = (int **)calloc(pad_dim, sizeof(int *));
   int **three = (int **)calloc(pad_dim, sizeof(int *));
   for (int i = 0; i < pad_dim; i++){
      one[i] = (int *)calloc(pad_dim, sizeof(int));
      two[i] = (int *)calloc(pad_dim, sizeof(int));
      three[i] = (int *)calloc(pad_dim, sizeof(int));
   }
   
   makeTwoMatrices(f, dim, one, two);
   clock_t start, end;
   
   start = clock();
   if (atoi(argv[1]) == 0){
      printf("big mult");
      big_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim, pad_dim);
   } else {
      standard_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim);
   }
   // big_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim, pad_dim);
   // strassen(one, two, three, 0, 0, 0, 0, 0, 0, dim, pad_dim);
   standard_mult(one, two, three, 0, 0, 0, 0, 0, 0, dim);
   end = clock();
   double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

   
   // printMatrices(three, pad_dim);
   // diagonals(three, dim);
   printf("TIME :%f \n", cpu_time_used);
   for (int i = 0; i < pad_dim; i ++){
      free(one[i]);
      free(two[i]);
      free(three[i]);
   }

   free(one);
   free(two);
   free(three);
   return 0;

}


void diagonals(int **three, int dim){
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

int** big_mult_onlyStrassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim){
   // check dimensions, if under crossover point, call standard, else call strassens
   return strassenOpt(one, two, three, r1, c1, r2, c2, r3, c3, real_dim, dim);
}

int** big_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim){
   // check dimensions, if under crossover point, call standard, else call strassens
   if (r1 >= real_dim || c1 >= real_dim){
      return one;
   } else if (r2 >= real_dim || c2 >= real_dim){
      return two;
   } else if (dim < 10){
      return standard_mult(one, two, three, r1, c1, r2, c2, r3, c3, dim);
   } else {
      return strassenOpt(one, two, three, r1, c1, r2, c2, r3, c3, real_dim, dim);
   }   
}

/*standard matrix multiplication
RETURNS: int**
tasks:
1) perform standard mult
*/
int** standard_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim){

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
    for (int i = 0; i < dim; i ++){
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
int **strassenOpt1(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim){
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
      for (int i = 0; i < 2; i++){
         sums[i] = (int **)calloc(new_d, sizeof(int*));
         for (int j = 0; j < new_d; j++){
            sums[i][j] = (int *)calloc(new_d, sizeof(int));
         }
      }


      int **p = (int**)calloc(new_d, sizeof(int*));
      for (int j = 0; j < new_d; j++){
         p[j] = (int *)calloc(new_d, sizeof(int));   
      }

      for (int i = 0; i < 10; i++){
         switch(i){
            case 0: // f- h
               matrixDestinationAddition(two, r2, c2 + new_d, two, r2 + new_d, c2 + new_d, sums[0], new_d, -1);
               big_mult(one, sums[0], p, r1, c1, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 , c3 + new_d, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, 1);
               break;

            case 1: // a + b
               for (int j = 0; j < new_d; j++){
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(one, r1, c1, one, r1, c1 + new_d, sums[0], new_d, 1);
               big_mult(sums[0], two, p, 0, 0, r2 + new_d, c2 + new_d, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 , c3 + new_d, new_d, 1);
               matrixCorAddition(p, three, r3, c3, new_d, -1);
               break;

            case 2: //c + d
               matrixDestinationAddition(one, r1 + new_d, c1, one, r1 + new_d, c1 + new_d, sums[0], new_d, 1);
               for (int j = 0; j < new_d; j++){
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               big_mult(sums[0], two, p, 0, 0, r2, c2, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 + new_d, c3, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, -1);
               break;

            case 3: //g - e
               for (int j = 0; j < new_d; j++){
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
               for (int j = 0; j < new_d; j++){
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
               for (int j = 0; j < new_d; j++){
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
               for (int j = 0; j < new_d; j++){
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(two, r2, c2, two, r2, c2 + new_d, sums[1], new_d, 1);
               big_mult(sums[0], sums[1], p, 0, 0, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, -1);
               break;
         }
      }
      
      for (int i = 0; i < 2; i++){
         for (int j = 0; j < new_d; j++){
               free(sums[i][j]);
         }
         free(sums[i]);
      }
      free(sums);

      for (int j = 0; j < new_d; j++){
         free(p[j]);
      }
      free(p);
      
      return three;
   } 
}   

int **strassenOpt(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim){
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
      for (int i = 0; i < 2; i++){
         sums[i] = (int **)calloc(new_d, sizeof(int*));
         for (int j = 0; j < new_d; j++){
            sums[i][j] = (int *)calloc(new_d, sizeof(int));
         }
      }


      int **p = (int**)calloc(new_d, sizeof(int*));
      for (int j = 0; j < new_d; j++){
         p[j] = (int *)calloc(new_d, sizeof(int));   
      }
      for (int i = 0; i < 10; i++){
         switch(i){
            case 0: // f- h
               matrixDestinationAddition(two, r2, c2 + new_d, two, r2 + new_d, c2 + new_d, sums[0], new_d, -1);
               big_mult(one, sums[0], p, r1, c1, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 , c3 + new_d, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, 1);
               break;

            case 1: // a + b
               for (int j = 0; j < new_d; j++){
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(one, r1, c1, one, r1, c1 + new_d, sums[0], new_d, 1);
               big_mult(sums[0], two, p, 0, 0, r2 + new_d, c2 + new_d, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 , c3 + new_d, new_d, 1);
               matrixCorAddition(p, three, r3, c3, new_d, -1);
               break;

            case 2: //c + d
               matrixDestinationAddition(one, r1 + new_d, c1, one, r1 + new_d, c1 + new_d, sums[0], new_d, 1);
               for (int j = 0; j < new_d; j++){
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               big_mult(sums[0], two, p, 0, 0, r2, c2, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 + new_d, c3, new_d, 1);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, -1);
               break;

            case 3: //g - e
               for (int j = 0; j < new_d; j++){
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
               for (int j = 0; j < new_d; j++){
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
               for (int j = 0; j < new_d; j++){
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
               for (int j = 0; j < new_d; j++){
                  memset(p[j], 0, sizeof(int) * new_d);
               }
               matrixDestinationAddition(two, r2, c2, two, r2, c2 + new_d, sums[1], new_d, 1);
               big_mult(sums[0], sums[1], p, 0, 0, 0, 0, 0, 0, real_dim, new_d);
               matrixCorAddition(p, three, r3 + new_d, c3 + new_d, new_d, -1);
               break;
         }
      }
      
      for (int i = 0; i < 2; i++){
         for (int j = 0; j < new_d; j++){
               free(sums[i][j]);
         }
         free(sums[i]);
      }
      free(sums);

      for (int j = 0; j < new_d; j++){
         free(p[j]);
      }
      free(p);
      
      return three;
   } 
}  
/*
int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int real_dim, int dim){
   int new_d = dim/2;
   if (new_d == 1){
      // actually multiply & add to three
      if (r1 >= real_dim || c1 >= real_dim || r2 >= real_dim || c2 >= real_dim){
         three[r3][c3] = 0;
         three[r3][c3 + 1]= 0;
         three[r3 + 1][c3] = 0;
         three[r3 + 1][c3 + 1] = 0;
      } else {
         int a = one[r1][c1], b = one[r1][c1 + 1], c = one[r1 + 1][c1], d = one[r1 + 1][c1 + 1];
         int e = two[r2][c2], f = two[r2][c2 + 1], g = two[r2 + 1][c2], h = two[r2 + 1][c2 + 1];
         int p1 = a*(f - h);
         int p2 = (a + b)*h;
         int p3 = (c + d)*e;
         int p4 = d*(g - e);
         int p5 = (a + d)*(e + h);
         int p6 = (b - d)*(g + h);
         int p7 = (a - c)*(e + f);
         three[r3][c3] = p5 + p4 - p2 + p6;
         three[r3][c3 + 1]= p1 + p2;
         three[r3 + 1][c3] = p3 + p4;
         three[r3 + 1][c3 + 1] = p1 + p5 - p3 - p7;
      }
      return three;
   } else {
      bool all_zeroes = false, one_zeroes = false, two_zeroes = false;
      
      if (r1 >= real_dim || c1 >= real_dim){
         one_zeroes = true;
      }
      if (r2 >= real_dim || c2 >= real_dim){
         two_zeroes = true;
      }
      if (one_zeroes & two_zeroes){
         all_zeroes = true;
      }

      if (!all_zeroes){
         int ***sums = (int***)calloc(10, sizeof(int**));
         for (int i = 0; i < 10; i++){
            sums[i] = (int **)calloc(new_d, sizeof(int*));
            for (int j = 0; j < new_d; j++){
               sums[i][j] = (int *)calloc(new_d, sizeof(int));
            }
         }
         for (int i = 0; i < 10; i++){
            switch(i){
               case 0: // f- h
                  if (!two_zeroes){
                     matrixDestinationAddition(two, r2, c2 + new_d, two, r2 + new_d, c2 + new_d, sums[i], new_d, -1);
                     
                  }
                  break;

               case 1: // a + b
                  if (!one_zeroes){
                     matrixDestinationAddition(one, r1, c1, one, r1, c1 + new_d, sums[i], new_d, 1);

                  }
                  break;

               case 2: //c + d
                  if (!one_zeroes){
                     matrixDestinationAddition(one, r1 + new_d, c1, one, r1 + new_d, c1 + new_d, sums[i], new_d, 1);
                     
                  }
                  break;

               case 3: //g - e
                  if (!two_zeroes){ 
                     matrixDestinationAddition(two, r2 + new_d, c2, two, r2, c2, sums[i], new_d, -1);
                  }
                  break;

               case 4: // a + d
                  if (!one_zeroes){ 
                     matrixDestinationAddition(one, r1, c1, one, r1 + new_d, c1 + new_d, sums[i], new_d, 1);
                  }
                  break;

               case 5: //e + h
                  if (!two_zeroes){ 
                     matrixDestinationAddition(two, r2, c2, two, r2 + new_d, c2 + new_d, sums[i], new_d, 1);
                  }
                  break;

               case 6: //b - d
                  if (!one_zeroes){
                     matrixDestinationAddition(one, r1, c1 + new_d, one, r1 + new_d, c1 + new_d, sums[i], new_d, -1);
                  }
                  break;

               case 7: // g + h
                  if (!two_zeroes){ 
                     matrixDestinationAddition(two, r2 + new_d, c2, two, r2 + new_d, c2 + new_d, sums[i], new_d, 1);
                  }
                  break;

               case 8: // a - c
                  if (!one_zeroes){
                     matrixDestinationAddition(one, r1, c1, one, r1 + new_d, c1, sums[i], new_d, -1);
                  }
                  break;

               default: //e + f
                  if (!two_zeroes){
                     matrixDestinationAddition(two, r2, c2, two, r2, c2 + new_d, sums[i], new_d, 1);
                  }
                  break;
            }
         }
         
         int ***p = (int ***)calloc(7, sizeof(int**));
         for (int i = 0; i < 7; i++){
            p[i] = (int**)calloc(new_d, sizeof(int*));
            for (int j = 0; j < new_d; j++){
               p[i][j] = (int *)calloc(new_d, sizeof(int));
            }
         }
 
         p[0] = big_mult(one, sums[0], p[0], r1, c1, 0, 0, 0, 0, real_dim, new_d);  
         p[1] = big_mult(sums[1], two, p[1], 0, 0, r2 + new_d, c2 + new_d, 0, 0, real_dim, new_d);
         p[2] = big_mult(sums[2], two, p[2], 0, 0, r2, c2, 0, 0, real_dim, new_d);
         p[3] = big_mult(one, sums[3], p[3], r1 + new_d, c1 + new_d, 0, 0, 0, 0, real_dim, new_d);
         p[4] = big_mult(sums[4], sums[5], p[4], 0, 0, 0, 0, 0, 0, real_dim, new_d);
         p[5] = big_mult(sums[6], sums[7], p[5], 0, 0, 0, 0, 0, 0, real_dim, new_d);
         p[6] = big_mult(sums[8], sums[9], p[6], 0, 0, 0, 0, 0, 0, real_dim, new_d);
         matrixThirdsAddition(4, 3, 1, 5, 1, 1, -1, 1, p, 7, three, r3, c3, new_d);
         matrixThirdsAddition(0, 1, -1, -1, 1, 1, 1, 1, p, 7, three, r3, c3 + new_d, new_d);
         matrixThirdsAddition(2, 3, -1, -1, 1, 1, 1, 1, p, 7, three, r3 + new_d, c3, new_d);
         matrixThirdsAddition(0, 4, 2, 6, 1, 1, -1, -1, p, 7, three, r3 + new_d, c3 + new_d, new_d);
      
         for (int i = 0; i < 10; i++){
            for (int j = 0; j < new_d; j++){
                  free(sums[i][j]);
            }
            free(sums[i]);
         }
         free(sums);
         for (int i = 0; i < 7; i++){
            for (int j = 0; j < new_d; j++){
               free(p[i][j]);
            }
            free(p[i]);
         }
         free(p);
      }
      return three;

   }
}
*/




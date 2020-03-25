#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.c"

int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim);

int main(int argc, char** argv) {
   // Ensure correct usage
   if (argc != 4) {
         printf("Usage: ./strassen 0 dimension inputfile \n");
         return 1;
   }
   printf("%s \n", argv[3]);
   FILE *f = fopen(argv[3], "r");
   if (f == 0)
   {
      //fopen returns 0, the NULL pointer, on failure
      perror("Cannot open input file\n");
      exit(-1);
	}
   printf("opened the file \n");
   int dim = atoi(argv[2]);
   printf("%i \n", dim);
   char ch, buffer[10];
   int **one = (int **)calloc(dim, sizeof(int *)); 
   int **two = (int **)calloc(dim, sizeof(int *));
   int **three = (int **)calloc(dim, sizeof(int *));
   for (int i = 0; i < dim; i++){
      one[i] = (int *)calloc(dim, sizeof(int));
      two[i] = (int *)calloc(dim, sizeof(int));
      three[i] = (int *)calloc(dim, sizeof(int));
   }
   makeTwoMatrices(f, dim, one, two);
   printf("\n \n");

   printf("\n \n");
   strassen(one, two, three, 0, 0, 0, 0, 0, 0, dim);
   printMatrices(three, dim);
   for (int i = 0; i < dim; i ++){
      free(one[i]);
      free(two[i]);
      free(three[i]);
   }

   free(one); 
   free(two);
   free(three);
   return 0;

}



/*big multiplication algorithm
RETURNS: a list of diagonal entries 
tasks:
1) turn file into matrix form (if necessary? is this necessary?)
2) based on n, either call strassens or call standard mult
3) collect diagonals
*/
int big_mult(){
		return 0;
}

/*standard matrix multiplication 
RETURNS: ?
tasks:
1) perform mult -> recursively call itself? or refer to big_mult
*/
int** standard_mult(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim){

}




/*strassens algorithm
RETURNS: ?
tasks:
1) add necessary padding so that it's of a dimension that's a power of two
2) perform strassens -> recursive call to big_mult in it
3) in strassen, add each term to the correct part of the third and final matrix
*/
int **strassen(int **one, int **two, int **three, int r1, int c1, int r2, int c2, int r3, int c3, int dim){
   int new_d = dim/2;
   if (new_d == 1){
      // actually multiply & add to three
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
      return three;
   } else {
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
               matrixDestinationAddition(two, r2, c2 + new_d, two, r2 + new_d, c2 + new_d, sums[i], new_d, -1);
               break;

            case 1: // a + b
               matrixDestinationAddition(one, r1, c1, one, r1, c1 + new_d, sums[i], new_d, 1);
               break;
            
            case 2: //c + d
               matrixDestinationAddition(one, r1 + new_d, c1, one, r1 + new_d, c1 + new_d, sums[i], new_d, 1);
               break;
            
            case 3: //g - e
               matrixDestinationAddition(two, r2 + new_d, c2, two, r2, c2, sums[i], new_d, -1);
               break;

            case 4: // a + d
               matrixDestinationAddition(one, r1, c1, one, r1 + new_d, c1 + new_d, sums[i], new_d, 1);
               break;

            case 5: //e + h
               matrixDestinationAddition(two, r2, c2, two, r2 + new_d, c2 + new_d, sums[i], new_d, 1);
               break;

            case 6: //b - d
               matrixDestinationAddition(one, r1, c1 + new_d, one, r1 + new_d, c1 + new_d, sums[i], new_d, -1);
               break;

            case 7: // g + h
               matrixDestinationAddition(two, r2 + new_d, c2, two, r2 + new_d, c2 + new_d, sums[i], new_d, 1);
               break;

            case 8: // a - c
               matrixDestinationAddition(one, r1, c1, one, r1 + new_d, c1, sums[i], new_d, -1);
               break;

            default: //e + f
               matrixDestinationAddition(two, r2, c2, two, r2, c2 + new_d, sums[i], new_d, 1);
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
      p[0] = strassen(one, sums[0], p[0], r1, c1, 0, 0, 0, 0, new_d); 
      p[1] = strassen(sums[1], two, p[1], 0, 0, r2 + new_d, c2 + new_d, 0, 0, new_d); 
      p[2] = strassen(sums[2], two, p[2], 0, 0, r2, c2, 0, 0, new_d); 
      p[3] = strassen(one, sums[3], p[3], r1 + new_d, c1 + new_d, 0, 0, 0, 0, new_d); 
      p[4] = strassen(sums[4], sums[5], p[4], 0, 0, 0, 0, 0, 0, new_d); 
      p[5] = strassen(sums[6], sums[7], p[5], 0, 0, 0, 0, 0, 0, new_d); 
      p[6] = strassen(sums[8], sums[9], p[6], 0, 0, 0, 0, 0, 0, new_d); 
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
      return three;

   }
}





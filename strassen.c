#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.c"


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
   int **one = (int **)malloc(dim * sizeof(int *)); 
   int **two = (int **)malloc(dim * sizeof(int *));
   int **three = (int **)malloc(dim * sizeof(int *));
   for (int i = 0; i < dim; i++){
      one[i] = (int *)malloc(dim * sizeof(int));
      two[i] = (int *)malloc(dim * sizeof(int));
      three[i] = (int *)malloc(dim * sizeof(int));
   }
   makeTwoMatrices(f, dim, one, two);
   printf("\n \n");
   matrixAddition(one, two, dim, 1);
   printMatrices(one, dim);
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
int* big_mult(){
		return 0;
}

/*standard matrix multiplication 
RETURNS: ?
tasks:
1) perform mult -> recursively call itself? or refer to big_mult
*/
void standard_mult(){

}
/*strassens algorithm
RETURNS: ?
tasks:
1) add necessary padding so that it's of a dimension that's a power of two
2) perform strassens -> recursive call to big_mult in it
3) in strassen, add each term to the correct part of the third and final matrix
*/
void strassen(){

}



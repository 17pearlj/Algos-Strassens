// #include <float.h>
// #include <math.h>
// #include <stdbool.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <time.h>


int main(int argc, char** argv) {

  // Ensure correct usage
  if (argc != 4) {
    printf("Usage: ./strassen 0 dimension inputfile \n");
    return 1;
  }
  
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
*/
void strassen(){

}

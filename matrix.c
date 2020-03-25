
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

void printMatrices(int **matrix, int dim);
void matrixAddition(int **one, int **two, int dim, int positive);
void matrixDestinationAddition(int **one, int r1, int c1, int **two, int r2, int c2, int**three, int r3, int c3, int dim, int positive);
int makeTwoMatrices(FILE *f, int dim, int** one, int** two);

int makeTwoMatrices(FILE *f, int dim, int** one, int** two){
	char ch, buffer[10];
	int i = 0, ar = 0, row = 0, col = 0, counter = 0;
		// while both arrays not read
	while (counter < (2*(pow(dim, 2.)))){
			// Reads char
		ch = fgetc(f);
	   // If EOF is encountered, error 
		if(ch == EOF){
			perror("File is too short");
			exit(-1);
			break;
		} else if(ch == '\n'){
			switch (ar){
				case 0: 
               one[row][col] = atoi(buffer);
					break;
				default: two[row][col] = atoi(buffer);
			}
			col++;
			if (col == dim){
				col = 0;
				row++;
				if (row == dim){
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
	
	printMatrices(one, dim);
    printf("\n \n");
	printMatrices(two, dim);
	return 0;

}

void printMatrices(int **matrix, int dim){
    for(int r = 0; r < dim; r++){
      printf("[ ");
		for (int c = 0; c < dim; c++){
			printf("%d, ", matrix[r][c]);
		}
      printf("] \n");
	}
}

void matrixAddition(int **one, int **two, int dim, int positive){
    for(int r = 0; r < dim; r++){
		for (int c = 0; c < dim; c++){
            one[r][c] += (two[r][c]*positive);
		}
	}
    
}
void matrixDestinationAddition(int **one, int r1, int c1, int **two, int r2, int c2, int**three, int r3, int c3, int dim, int positive){
    for(int r = 0; r < dim; r++){
		for (int c = 0; c < dim; c++){
            three[r + r1][c + c1] = one[r + r1][c + c1] + (two[r + r1][c + c1]*positive);
		}
	}
    
}
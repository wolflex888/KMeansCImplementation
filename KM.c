/*
 * main.c
 *
 *  Created on: Sep 15, 2017
 *      Author: DavidLu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>




// define functions //
double**      	malloc_matrix   (int n1, int n2);
void 			free_matrix(int n1, int n2, double **a);
double 			cal_dis(double *x, double *y, int dim);
void 			cluster(double **data, double **cent, int numData, int numDim, int k, int *cent_count, double **cent_sum);
void 			update_centroid(double **cent, int k, double **cent_sum, int numData, int numDim, int *cent_count);
double 			sse_culc(double *a, double *b, int numDim);

int main(int argc, const char * argv[]){
	if (argc != 4 || argv[1] == NULL || argv[2] == NULL || argv[3] == NULL) { //check for invalid input
	        printf("please give me three parameters\n1.input file path\n2.number of cluster\n3.output file path\n"); // usage
	        printf("path to executable: %s\n", argv[0]);
	        printf("%d\n", argc);
	        exit(1);
	    }
	FILE* fp;
	FILE* ofp;
	// put argument into right variable //
	const char* inputFile = argv[1];
	int k = atoi(argv[2]);
	const char* outputFile = argv[3];

	// create file handle //
	fp = fopen(inputFile, "r");
// check if the file is empty //
	if(fp == NULL){
		printf("inputfile io error");
		return 1;
	}
	int numData; // load the first line of the data //
	int numDim;
	if (fscanf(fp, "%d %d", &numData, &numDim) !=2){
		printf("FILE FORMAT ERROR");
	}
	// declair necessary arrays for necessary data //
	double **data = malloc_matrix(numData, numDim+1); // initiate data matrix
	double **cent = malloc_matrix(k, numDim); // initiate centroid matrix
	double **prev_cent = malloc_matrix(k, numDim);
	double **cent_sum = malloc_matrix(k, numDim);
	int *cent_count = (int *)malloc(k * sizeof(int));

	// load data to the matrix //
	for(int i = 0; i<numData; i++){ // load data into matrix
			for(int j = 0; j<numDim; j++){
				fscanf(fp, "%le", &data[i][j]);
				if(i < k){
					cent[i][j] = data[i][j]; // initialize the centroid as the first k points
				}
			}
		}
	double sse = 0; // initialize sse
	int iteration = 0; // count iteration
	int flag = 0; // make sure the algorithm won't jump out at the first iteration
	do{
		sse = 0;
		flag = 1;
		for(int i= 0; i < k; i++){
						cent_count[i] = 0;
						for(int j = 0; j < numDim; j++){
							cent_sum[i][j] = 0; // reset the matrix
						}
					}

			cluster(data, cent, numData, numDim, k, cent_count, cent_sum); // step 1 -> assign cluster
			for(int i = 0; i < k; i++){
				for(int j = 0; j < numDim; j++){
					prev_cent[i][j] = cent[i][j]; // store current centroid for future calculation
					//printf("%f\t", cent[i][j]);
				}
			}
			// reset cent_sum, cent_count

			update_centroid(cent, k, cent_sum, numData, numDim, cent_count); // according to the cluster, update centroid location
			for(int i = 0; i < k; i++){
				double count = (double)cent_count[i]; //get count in double for future calculation
				//printf("count is %f\n", count);
				double temp = sse_culc(prev_cent[i], cent[i], numDim);
				//printf("before dividing is %f\n", temp);
				sse += sse_culc(prev_cent[i], cent[i], numDim)/count;
			}
			printf("%f\n", sse);
			iteration++;

			printf("=============== iteration %d ================\n", iteration);
	}while(sse > 0.0001 && iteration < 1000 && flag == 1); // threshold set as 0.0001 and max iteration set at 1000

	ofp = fopen(outputFile, "w"); // get output handler
	fprintf(ofp, "%d %d %d\n", numData, numDim, k); // write the first line number of data, number of dimension and number of cluster
	for(int i = 0; i < numData; i++){ //
		for (int j = 0; j < numDim +1; j++){
			fprintf(ofp, "%lf ", data[i][j]);// write everything into file
		}
		fprintf(ofp, "\n");
	}
////////////// evaluation "the elbow test" //////////////////////
	double eval = 0;
	for(int i = 0; i < numData; i++){
		printf("initial eval.. %f\n", eval);
		int cluster = (int)data[i][numDim];
		printf("cluster: %f\n", data[i][numDim]);
		printf("cluster_act: %d\n", cluster);
		eval += sse_culc(data[i], cent[cluster], numDim);
		for(int j = 0; j < numDim; j++){
			printf("%f\t%f\n", cent[cluster][j], data[i][j]);
		}
		printf("eval now.. %f\n", eval);
		}
	printf("Sum of Sqaure Error for k value %d is %f\n", k, eval);
///////////////////////////////////////////////////////////////////

	// free memory //
	free_matrix(k, numDim, cent_sum);
	free_matrix(numData, numDim, data);
	free_matrix(k, numDim, cent);
	free(cent_count);
///////////////////
}
// calculate sse // given rows and dimension of the row
double sse_culc(double *a, double *b, int numDim){
	double sse = 0;
	for(int i = 0; i < numDim; i++){
		double r = (a[i] - b[i]);
		sse += r*r;
	}


	return sse;
}

// calculate center according to the current cluster points//
void update_centroid(double **cent, int k, double **cent_sum, int numData, int numDim, int *cent_count){
	for(int i = 0; i < k; i++){
		for(int j = 0; j < numDim; j++){
			//printf("%f\n", cent[i][j]);
			cent[i][j] = cent_sum[i][j] / cent_count[i];
			//printf("%f\n", cent[i][j]);
		}
	}
}

// assign cluster using euclidean distance //
void cluster(double **data, double **cent, int numData, int numDim, int k, int *cent_count, double **cent_sum){
	for(int i = 0; i < numData; i++){
			double min = INFINITY; // set the minimum to infinity//
			int cluster;
			for(int l = 0; l<k; l++){
				double temp = cal_dis(data[i], cent[l], numDim);
				printf("distance is: %f\n", temp);
				if(min > temp){
					//printf("yoyoyo, %f\n", temp);
					min = temp;
					cluster = l;
				}
				//printf("%d %f\n", l, min);

			}
			//printf("%d\n", cluster);
			data[i][numDim] = cluster; // write the cluster to the last column of the data //
			cent_count[cluster]++;
			for(int j = 0; j < numDim; j++){
				cent_sum[cluster][j] += data[i][j];
			}
		}
}
// calculate distance // .... i realize this is the same as sse calculation... please don't laugh at me.. it's 3 in the morning
double cal_dis(double *x, double *y, int dim){
	int i;
	double t = 0;
	double sum=0;
	for(i = 0; i<dim; i++){
		t=x[i]-y[i], sum+=t*t;
	}
	return sum;
}

// allocate memory //
double **malloc_matrix(int n1, int n2) {

    double **mat = NULL;       // pointer to the matrix
    if(n1 <= 0 || n2 <= 0 || n1 == NULL || n2 == NULL){ // check for invalid input
    			 printf("invalid input\n");
    			 exit(0);
    		}
    mat = (double **)malloc(n1 * sizeof(double *)); // allocate memory for the matrix
    for (int i = 0; i < n1; i++){
    		mat[i] = (double *)malloc(n2*sizeof(double)); //allocate second dimension
    }
    return mat;
}

void free_matrix(int n1, int n2, double **a) {

	if(n1 <= 0 || n2 <= 0 || n1 == NULL || n2 == NULL){ // check for invalid input
		 printf("invalid input\n");
		 exit(0);
	}
	for (int i = 0; i < n1; i++){
		free(a[i]); //free second demension
	}
	free(a); // free first dimension
}
















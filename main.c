#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "math.h"
#include "complex.h"

#define MANDELBROT 0 
#define JULIA 1 
#define NUM_COLORS 256

//structure that is sent from master to slaves
typedef struct{
	int type;
	int width;
	int height;
	double minOx, maxOx, minOy, maxOy;
	double resolution;
	double reJulia, imJulia;
	int steps;
}Buffer;


int main(int argc, char **argv){
	if(argc != 3) {
		printf("USAGE: ./tema3 input_file output_image");
		return 0;
	}
	int width;
	int i,j,k;
	int height;
	int type;
	double minOx, maxOx, minOy, maxOy;
	double resolution;
	double reJulia, imJulia;
	int steps;
	int step;
    
	int num_tasks, rank, len, rc = 0;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("ERROR: could not start MPI program. QUIT.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return 0;
	}
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank == 0){
		//Master Process
		//Reading from input file
		FILE *in = fopen(argv[1], "r");
		if(!in) {
			printf("ERROR: while opening %s \n", argv[1]);
		}	
		fscanf(in, "%d", &type);
		fscanf(in, "%lf %lf %lf %lf", &minOx, &maxOx, &minOy, &maxOy);
		fscanf(in, "%lf", &resolution);
		fscanf(in, "%d", &steps);
		if(type == JULIA) {
			fscanf(in, "%lf %lf", &reJulia, &imJulia);
		} else if(type == MANDELBROT){
			reJulia = 0.0;
			imJulia = 0.0;
		}
		fclose(in); 
		width = (maxOx - minOx)/resolution;
		height = (maxOy - minOy)/resolution;
		int matrix[height][width];
		//Parsing the first part of the input(the master part)
		for(i = 0; i < height/num_tasks;i++){
			for(j = 0; j < width;j++){
				step = 0;
				double complex z;
				double complex c;
				if(type == JULIA){
					z = (minOx + j * resolution) + (minOy + i * resolution)* I;
					c = reJulia + imJulia * I;
				}else if(type == MANDELBROT){
					c = (minOx + j * resolution) + (minOy + i * resolution)* I;
					z = 0 + 0 * I;
				}
				while((sqrt(creal(z)*creal(z) + cimag(z)*cimag(z)) < 2) && (step < steps)){
					z = z * z + c;
					step++;
				}
				int color = step % NUM_COLORS;
				matrix[i][j] = color;
			}
		}
		//Sending Buffer to slaves
		for (i = 1; i < num_tasks; ++i) {
			
			Buffer b;
			b.type = type;
			b.width = width;
			b.height = height;
			b.minOx = minOx;
			b.minOy = minOy;
			b.maxOx = maxOx;
			b.maxOy = maxOy;
			b.resolution = resolution;
			b.reJulia = reJulia;
			b.imJulia = imJulia;
			b.steps = steps;
			MPI_Send(&b, sizeof(Buffer), MPI_BYTE, i, 1, MPI_COMM_WORLD);
		}
		//Receiving info from slaves and updating the main matrix
		for (k = 1; k < num_tasks; ++k) {
			int rmatrix[height/num_tasks][width];
			MPI_Recv(&rmatrix, (height/num_tasks)*width, MPI_INT, k, 1, MPI_COMM_WORLD,&status);
			for(i = 0; i < height/num_tasks;i++){
				for(j = 0; j < width;j++){
					matrix[i+k*height/num_tasks][j] = rmatrix[i][j];
				}
			}
		}
		//Writing results in output file
		FILE *f;
		f = fopen(argv[2], "wt");
		int i, j;
		fprintf(f, "%s\n","P2");
		fprintf(f, "%d %d\n", width, height);
		fprintf(f, "%d\n", NUM_COLORS - 1);
		for(i = 0; i < height; ++i) {
			for(j = 0; j < width; ++j){
				fprintf(f, "%d ", matrix[height - i - 1][j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}else{
		//Slave processes
		//Receiving the buffer from master
		Buffer r;
		MPI_Recv(&r, sizeof(Buffer), MPI_BYTE, 0, 1, MPI_COMM_WORLD, &status);
		int matrix[r.height/num_tasks][r.width];
		//Parsing a part of the matrix(the  ones remaining for the slaves)
		for(i = 0; i < r.height/num_tasks;i++){
			for(j = 0; j < r.width;j++){
				step = 0;
				double complex c;
				double complex z;
				if(r.type == JULIA){
					c = r.reJulia + r.imJulia * I;
					z = (r.minOx + j * r.resolution) + (r.minOy + (i + ((rank * r.height )/ num_tasks))* r.resolution)* I;
				}else if(r.type == MANDELBROT){
					z = 0 + 0 * I;
					c = (r.minOx + j * r.resolution) + (r.minOy + (i + ((rank * r.height )/ num_tasks))* r.resolution)* I;
				}
				while((sqrt(creal(z)*creal(z) + cimag(z)*cimag(z)) < 2) && (step < r.steps)){
					z = z * z + c;
					step++;
				}
				int color = step % NUM_COLORS;
				matrix[i][j] = color;
			}
		}
		MPI_Send(&matrix, (r.height/num_tasks)*r.width, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}
	MPI_Get_processor_name(hostname, &len);
	MPI_Finalize();
	return 0;
}

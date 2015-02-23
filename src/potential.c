#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int count_lines(char *filename);

// Gives an array with the gravitational potential for each particle
// Requires:
//    argv[1]: name of the file with the positions of each particle in x,y and z
//    argv[2]: name of the file to export the potential of each particle
// Returns:
//    a file with the potential for each particle

int main(int argc, char** argv){

 int i,j,n_points;
 double count;
 double* x;
 double* y;
 double* z;
 double* k;
 FILE *rfile;
 FILE *cfile;

 n_points = count_lines(argv[1]);

 rfile = fopen(argv[1], "r");
 cfile = fopen(argv[2], "w");

 x = malloc(n_points*sizeof(double));
 y = malloc(n_points*sizeof(double));
 z = malloc(n_points*sizeof(double));
 k = malloc(n_points*sizeof(double));

 for (i = 0; i < n_points; i ++){
   fscanf(rfile,"%lf,%lf,%lf\n",&(x[i]),&(y[i]),&(z[i]));
   k[i] = 0;
 }

for (i = 0; i < n_points; i ++){
   for (j = i+1; j < n_points; j ++){
     count = 1.0/sqrt(pow(x[i]-x[j],2.0) + pow(y[i]-y[j],2.0) + pow(z[i]-z[j],2.0));
     k[i] += count;
     k[j] += count;    
   }
 }
 for (i = 0; i < n_points; i ++){
   fprintf(cfile,"%lf\n",k[i]);
 }
}

// Counts lines of a file, copied from github.com/forero/ComputationalPhysicsUniandes repository
// Requires
//    filename: a string with the name of the file
// Returns:
//    number of lines in the file
int count_lines(char *filename){
  FILE *in;
  int n_lines;
  int c;
  if(!(in=fopen(filename, "r"))){
    printf("problem opening file %s\n", filename);
    exit(1);
  }

  n_lines = 0;
  do{
    c = fgetc(in);
    if(c=='\n'){
      n_lines++;
    }
  }while(c!=EOF);

  rewind(in);
  fclose(in);
  return n_lines;
}

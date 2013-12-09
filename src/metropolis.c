#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int count_lines(char *filename);
double rand_generator(double min, double max);
double random_normal(double sigma, double mu);
double chisq(double* obs, double* mod, int n_points);
double loglogmass(double* logr, double logrho0, double logrs, int n_points);

int main()
{
  srand48(getpid());
  return 0;
}

int count_lines(char *filename)
{
  FILE *in;
  int n_lines;
  int c;
  if(!(in=fopen(filename, "r")))
  {
    printf("problem opening file %s\n", filename);
    exit(1);
  }

  n_lines = 0;
  do
  {
    c = fgetc(in);
    if(c=='\n')
    {
      n_lines++;
    }
  }
  while(c!=EOF);
  
  rewind(in);
  fclose(in);
  return n_lines;
}

double rand_generator(double min, double max)
{
  double ans;
  ans = min + drand48()*(max - min) 
  return ans;
}

double random_normal(double sigma, double mu)
{
  double ans;
  ans = rand_generator(mu-sigma,mu+sigma)
  return ans;
}

double chisq(double* obs, double* mod, int n_points)
{
  int i;
  double ans;
  ans = 0;
  for(i=0;i<n_points,i++)
  {
    ans -= pow((mod[i]-obs[i]),2)/2;
  }
  return ans;
}

double loglogmass(double* logr, double logrho0, double logrs, int n_points)
{
  int i;
  double* ans;
  ans = malloc(n_points*sizeof(double));
  for(i=0;i<n_points,i++)
  {
    ans[i] = log(4*PI_M)+3*logrs+logrho0+log(log(1+exp(logr-logrs))-pow(1+exp(logrs-logr),-1));
  }
  return ans;
}

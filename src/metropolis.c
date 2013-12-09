#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int count_lines(char *filename);
double rand_generator(double min, double max);
double random_normal(double sigma, double mu);
double chisq(double* obs, double* mod, int n_points);
double loglogmass(double* logr, double logrho0, double logrs, int n_points);

int main(int argc, char** argv)
{
  srand48(getpid());
  int i;
  int n_points;
  int n_iterations;

  double guess_a;
  double guess_b;

  double* a_walk;
  double* b_walk;
  double* chi2;

  FILE *data;
  FILE *a_file;
  FILE *b_file;
  FILE *chi2_file;

  double* logR;
  double* logM;
 
  n_points = count_lines(argv[1]);
  n_iterations = atoi(argv[4]);

  guess_a = atof(argv[2]);
  guess_b = atof(argv[3]);

  a_walk = malloc(n_iterations*sizeof(double));
  b_walk = malloc(n_iterations*sizeof(double));
  chi2 = malloc(n_iterations*sizeof(double));

  data = fopen(argv[1], "r");
  a_file = fopen("a_walk.dat", "w");
  b_file = fopen("b_walk.dat", "w");
  chi2_file = fopen("chi2.dat", "w");
  
  for(i=0;i<n_points,i++)
  {
    fscanf(data,"%lf %lf\n",&(logR[i]),&(logM[i]));
  }

  a_walk[0] = guess_a;
  b_walk[0] = guess_b;
  chi2[0] = chisq(logM,loglogmass(logR,guess_a,guess_b,n_points),n_points);

  for(i=0;i<n_iterations,i++)
  {
    double a_prime,b_prime,chi2_init,chi2_prime,alpha,ratio,beta;
    
    a_prime = random_normal(a_walk[i],0.001);
    b_prime = random_normal(b_walk[i],0.001);
    
    chi2_init = chisq(logM,loglogmass(logR,a_walk[i],b_walk[i],n_points),n_points);
    chi2_prime = chisq(logM,loglogmass(logR,a_prime,b_prime,n_points),n_points);
    
    alpha =  exp(chi2_prime-chi2_init);
    ratio = chi2_init/chi2_prime;
    
    if(ratio >= 1)
    {
      a_walk[i+1] = a_prime;
      b_walk[i+1] = b_prime;
      chi2[i+1] = chi2_prime;
    }
    else
    {
      beta = random(0,1);
      if(alpha >= beta)
      {
	a_walk[i+1] = a_prime;
	b_walk[i+1] = b_prime;
	chi2[i+1] = chi2_prime;
      }
      else
      {
	a_walk[i+1] = a_walk[i];
	b_walk[i+1] = b_walk[i];
	chi2[i+1] = chi2[i];
      }
    }
  }
  
  for (i = 0; i < n_iterations; i ++)
  {
    fprintf(a_file,"%lf\n",a_walk[i]);
    fprintf(b_file,"%lf\n",b_walk[i]);
    fprintf(chi2_file,"%lf\n",chi2[i]);
  } 
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

/* 2D Jacobi smoothing to solve -u'' = f */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

      
int main (int argc, char **argv)
{
  long  i,j,n;
  double *u, *unew, *f;
  double res = 1, h2, tmp=0.0, stopping_criteria;
  long p=0, max_iters = 1000;
  int nthreads = 8;
  
  if (argc != 3){
    n = 100;
  }
  else{
    n = atol(argv[1]);
    nthreads  = atol(argv[2]);}
  
#ifdef _OPENMP
omp_set_num_threads(nthreads);
#endif
 
#pragma omp parallel
  {
#ifdef _OPENMP
    int my_threadnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
#else
    int my_threadnum = 0;
    int numthreads = 1;
#endif
    printf("Thread %d of %d starting...\n", my_threadnum, numthreads);
  }
 
  h2 =  1 / (double) pow((n+1),2);
  stopping_criteria = sqrt((double) n*n) / (double) pow(10, 4);

  /* Allocation of vectors, including boundary ghost points*/
  u = (double *) malloc(sizeof(double) * pow((n+2),2));
  unew = (double *) malloc(sizeof(double) * pow((n+2),2));
  f = (double *) malloc(sizeof(double) * pow((n+2),2));

  /* Initialize all the vectors */
#pragma omp parallel for default(none) shared(n,unew,u,f) private(i)
  for (i = 0; i < (n+2)*(n+2)-1; ++i) {
    u[i] = 0.0;
    unew[i] = 0.0;
    f[i] = 1.0 ;
  }

  timestamp_type time1, time2;
  get_timestamp(&time1);


  
  while( p < max_iters && res > stopping_criteria){
    
/*  Update unew */
    
#pragma omp parallel for default(none) shared(n,unew,u,f, h2) private(i,j)
    for (i = 1; i < n+1; ++i){
      for (j = 1; j < n+1; ++j){
	unew[j+(n+2)*i] = (h2 * f[j+(n+2)*i] + u[j+(n+2)*i+1] + u[j+(n+2)*i-1]+u[j+(n+2)*(i-1)]+u[j+(n+2)*(i+1)]) / 4.0f;
      }
    }		 

/* Compute residual */
    res = 0.0;
#pragma omp parallel for default(none) shared(unew,f, n, h2) private(i,j,tmp) reduction(+:res)
      for (i = 1; i < n+1; ++i){
	for (j = 1; j < n+1; ++j){
	  tmp = (double)pow((f[j+(n+2)*i]- (4.0f* unew[j+(n+2)*i] - unew[j+(n+2)*i+1] - unew[j+(n+2)*i-1]-unew[j+(n+2)*(i-1)]-unew[j+(n+2)*(i+1)])/h2),2);
	  res += tmp;
	}
      }

    
/* Pass the value in unew onto u */
      double *utemp;
      utemp = u;
      u = unew;
      unew = utemp;
      
    p = p + 1;
  }

  printf("The number of grid points: N = %ld \n", n);
  printf("The total number of iteration is: %ld \n", p);
  printf("The norm of the final iteration is: %e \n", res);
  
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Time elapsed is %f seconds.\n", elapsed);

  free(u);
  free(unew);
  free(f);
  
  return 0;
}

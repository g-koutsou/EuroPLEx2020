#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

/* Allows modifying at compile time, e.g. -DNREP=20 */
#ifndef NREP
#define NREP 20
#endif

/* Speed of light */
#define C (3e8)

/* 
   Structure which holds space-time coords 
*/
typedef struct {
  float x;
  float y;
  float z;
  float t;  
  float s;
} st_coords;

/*
 * Returns seconds elapsed since t0
 */
double
stop_watch(double t0)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  double t1 = tp.tv_sec + tp.tv_usec*1e-6;  
  return t1-t0;
}

/*
 * Allocate memory with minimal error detection
 */
void *
alloc(size_t size)
{
  void *ptr;
  posix_memalign(&ptr, 64, size);
  if(ptr == NULL) {
    fprintf(stderr, "malloc() returned NULL, quitting\n");
    exit(3);
  }
  return ptr;
}

/*
 * Print usage (to stderr)
 */
void
usage(char *argv[])
{
  fprintf(stderr,
	  " Usage: %s <L>\n", argv[0]);
  return;
}

/*
 * compute: s = t**2 - x**2 - y**2 - z**2, with s, t, x, y, z
 * member variables of the st_coords structure arr
 */
void
comp_s(st_coords *arr, int L)
{
#pragma omp parallel for
  for(int i=0; i<L; i++) {
    float x = arr[i].x;
    float y = arr[i].y;
    float z = arr[i].z;
    float t = arr[i].t;
    float s = t*t - (x*x + y*y + z*z);
    arr[i].s = s;
  }
  return;
}

int
main(int argc, char *argv[])
{
  if(argc != 2) {
    usage(argv);
    return 1;
  }  
  int L = atoi(argv[1]);
  st_coords *arr = alloc(sizeof(st_coords)*L);
  for(int i=0; i<L; i++) {
    arr[i].x = drand48();
    arr[i].y = drand48();
    arr[i].z = drand48();
    arr[i].t = drand48()*C;
  }

  {
    /* Warm-up */
    comp_s(arr, L);
    double t0acc = 0;
    double t1acc = 0;
    int n = 1;
    /* 
       Loop accumulating run-time. Stop when the average time has less
       than a 10% error
    */
    while(1) {
      double t0 = stop_watch(0);
      for(int i=0; i<NREP; i++)
	comp_s(arr, L);
      t0 = stop_watch(t0)/(double)NREP;
      t0acc += t0;
      t1acc += t0*t0;
      if(n > 2) {
	double ave = t0acc/n;
	double err = sqrt(t1acc/n - ave*ave)/sqrt(n);
	if(err/ave < 0.1) {
	  t0acc = ave;
	  t1acc = err;
	  break;
	}
      }
      n++;
    }
    /*
     * t0acc is the average time, t1acc is its error
     */
    double beta_fp = (7.0*L)/t0acc/1e9               /* _TODO_A_: insert number of Gflop/sec based on timing */;
    double beta_io = (5.0*sizeof(float)*L)/t0acc/1e9 /* _TODO_A_: insert number of Gbyte/sec based on timing */;
    int nthr = 0;
#pragma omp parallel
    {
      nthr = omp_get_num_threads();
    }
    printf(" Nthr = %d L = %d in %3.1e +/- %3.1e secs, %g Gflop/s, %g Gbytes/s\n",
	   nthr, L, t0acc, t1acc, beta_fp, beta_io);
  }

  free(arr);
  return 0;
}

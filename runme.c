#include <omp.h>
#include <stdio.h>

int
main(int argc, char *argv[])
{
#pragma omp parallel
  {
    int n_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();
    printf(" n_threads = %d thread_id = %d\n", n_threads, thread_id);
  }
  return 0;
}

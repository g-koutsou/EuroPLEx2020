#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <complex.h>

#define ND 2
#define ALIGNMENT 32
#ifndef VL
#define VL 4
#endif

/*
 * Structure for scalar field
 */
typedef struct {
  float phi_re[VL];
  float phi_im[VL];
} field;

/*
 * Structure for gauge links
 */
typedef struct {
  float u_re[ND][VL];
  float u_im[ND][VL];
} link;

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
 * malloc with minimal error detection
 */
void *
alloc(size_t size)
{
  void *ptr;
  posix_memalign(&ptr, ALIGNMENT, size);
  if(ptr == NULL) {
    fprintf(stderr, " malloc() returned NULL. Out of memory?\n");
    exit(-1);
  }
  return ptr;
}

/*
 * allocates a new field and returns its starting address
 */
field *
new_field(size_t L)
{
  field *ptr = alloc(sizeof(field)*L*L/VL);
  return ptr;
}

/*
 * allocates a new U(1) gauge field and returns its starting address
 */
link *
new_links(size_t L)
{
  link *ptr = alloc(sizeof(link)*L*L/VL);
  return ptr;
}

/*
 * Fills u with random entries on the unit circle, gaussianly
 * distributed around 1 + 0*i, using Box-Mueller
 *
 * Sites are traversed in the same order as the non-vectorized
 * version, so that an identical random gauge field is produced
 * between the codes
 */
void
rand_links(size_t L, link *g)
{
  for(int y=0; y<L; y++) {
    int y0 = y % (L/VL);
    int y1 = y / (L/VL);
    for(int x=0; x<L; x++) {
      int v = y0*L + x;
      for(int d=0; d<ND; d++) {
	double u0 = drand48();
	double u1 = drand48();
	double phi = sqrt(-2*log(u0))*sin(2*M_PI*u1)*M_PI;
	g[v].u_re[d][y1] = cos(phi);
	g[v].u_im[d][y1] = sin(phi);
      }
    }
  }
    
  return;
}

/*
 * Fills x with random entries on the unit circle, gaussianly
 * distributed around 1 + 0*i, using Box-Mueller
 *
 * Sites are traversed in the same order as the non-vectorized
 * version, so that an identical field is produced between the codes
 */
void
rand_field(size_t L, field *z)
{
  for(int y=0; y<L; y++) {
    int y0 = y % (L/VL);
    int y1 = y / (L/VL);
    for(int x=0; x<L; x++) {
      int v = y0*L + x;
      double u0 = drand48();
      double u1 = drand48();
      double phi = sqrt(-2*log(u0))*sin(2*M_PI*u1)*M_PI;
      z[v].phi_re[y1] = cos(phi);
      z[v].phi_im[y1] = sin(phi);
    }
  }
  return;
}

/*
 * Fills x with zeros
 */
void
zero_field(size_t L, field *x)
{
  for(int i=0; i<L*(L/VL); i++) 
    for(int j=0; j<VL; j++) {
      x[i].phi_re[j] = 0.0;
      x[i].phi_im[j] = 0.0;
    }
  return;
}

/*
 * Shift the array of length VL by one element up
 */
static void
shift_up(float x[VL])
{
  float swap = x[0];
  for(int j=0; j<VL-1; j++)
    x[j] = x[j+1];
  x[VL-1] = swap;
  return;
}

/*
 * Shift the array of length VL by one element down
 */
static void
shift_dn(float x[VL])
{
  float swap = x[VL-1];
  for(int j=0; j<VL-1; j++)
    x[VL-j-1] = x[VL-j-2];
  x[0] = swap;
  return;
}

/*
 * Applies U(1) gauge laplacian to phi_in, with background field u,
 * and returns in phi_out
 */
void
lapl(size_t L, field *out, field *in, link *g)
{
#pragma omp parallel
  {
    float NDx2 = 2*ND;
    int ly = L/VL;
#pragma omp for    
    for(int y=0; y<ly; y++) {
      int y0 = y*L;
      int yp = ((y + 1)%ly)*L;
      int ym = ((y + ly - 1)%ly)*L;
      for(int x=0; x<L; x++) {
	int v00 = x + y0;
	int vp0 = x + yp;
	int vm0 = x + ym;
	int v0p = (x+1)%L + y0;
	int v0m = (L+x-1)%L + y0;

	float *in_v0p_re = &in[v0p].phi_re[0];
	float *in_v0p_im = &in[v0p].phi_im[0];
	float *in_v0m_re = &in[v0m].phi_re[0];
	float *in_v0m_im = &in[v0m].phi_im[0];

	float in_vp0_re[VL];
	float in_vp0_im[VL];
	float in_vm0_re[VL];
	float in_vm0_im[VL];
	for(int j=0; j<VL; j++) {
	  in_vp0_re[j] = in[vp0].phi_re[j];
	  in_vp0_im[j] = in[vp0].phi_im[j];
	  in_vm0_re[j] = in[vm0].phi_re[j];
	  in_vm0_im[j] = in[vm0].phi_im[j];
	}

	float *g_v00_0_re = &g[v00].u_re[0][0];
	float *g_v00_0_im = &g[v00].u_im[0][0];
	float *g_v00_1_re = &g[v00].u_re[1][0];
	float *g_v00_1_im = &g[v00].u_im[1][0];

	float *g_v0m_1_re = &g[v0m].u_re[1][0];
	float *g_v0m_1_im = &g[v0m].u_im[1][0];
	float g_vm0_0_re[VL];
	float g_vm0_0_im[VL];
	for(int j=0; j<VL; j++) {
	  g_vm0_0_re[j] = g[vm0].u_re[0][j];
	  g_vm0_0_im[j] = g[vm0].u_im[0][j];
	}

	/*
	 * _TODO_B_
	 *
	 * If y == 0, we need to use one of the shift functions on the
	 * affected fields
	 *
	 */
	if(y == 0) {
	  /*
	   * use shift_dn() or shift_up() here accordingly
	   */
	  shift_dn(g_vm0_0_re);
	  shift_dn(g_vm0_0_im);
	  shift_dn(in_vm0_re);
	  shift_dn(in_vm0_im);
	}
	
	/*
	 * _TODO_B_
	 *
	 * If y == ly-1, we need to use one of the shift functions on
	 * the affected fields
	 *
	 */
	if(y == ly-1) {
	  /*
	   * use shift_dn() or shift_up() here accordingly
	   */
	  shift_up(in_vp0_re);
	  shift_up(in_vp0_im);
	}
	
	for(int j=0; j<VL; j++) {
	  float p_re = 0;
	  float p_im = 0;
	  /*
	   * _TODO_B_
	   *
	   * Complete the missing parts using the structures in the
	   * new data layout. 
	   *
	   */
	  p_re  =   in_v0p_re[j]*g_v00_1_re[j] - in_v0p_im[j]*g_v00_1_im[j]; /* Direction x+, real part */
	  p_im  =   in_v0p_re[j]*g_v00_1_im[j] + in_v0p_im[j]*g_v00_1_re[j]; /* Direction x+, imag part */
		                                                           
	  p_re +=   in_v0m_re[j]*g_v0m_1_re[j] + in_v0m_im[j]*g_v0m_1_im[j]; /* Direction x-, real part */
	  p_im +=  -in_v0m_re[j]*g_v0m_1_im[j] + in_v0m_im[j]*g_v0m_1_re[j]; /* Direction x-, imag part */
	  	                                                           
	  p_re +=   in_vp0_re[j]*g_v00_0_re[j] - in_vp0_im[j]*g_v00_0_im[j]; /* Direction y+, real part */
	  p_im +=   in_vp0_re[j]*g_v00_0_im[j] + in_vp0_im[j]*g_v00_0_re[j]; /* Direction y+, imag part */
		                                                           
	  p_re +=   in_vm0_re[j]*g_vm0_0_re[j] + in_vm0_im[j]*g_vm0_0_im[j]; /* Direction y-, real part */
	  p_im +=  -in_vm0_re[j]*g_vm0_0_im[j] + in_vm0_im[j]*g_vm0_0_re[j]; /* Direction y-, imag part */
	  
	  out[v00].phi_re[j] = NDx2*in[v00].phi_re[j] - p_re;
	  out[v00].phi_im[j] = NDx2*in[v00].phi_im[j] - p_im;
	}
      }
    }
  }
  return;
}

/*
 * returns y = x^H x for vector x of length L*L
 */
double
xdotx(size_t L, field *x)
{
  double y = 0; 
  for(int i=0; i<L*(L/VL); i++)
    for(int j=0; j<VL; j++) {
      float phi_re = x[i].phi_re[j];
      float phi_im = x[i].phi_im[j];
      y += phi_re*phi_re + phi_im*phi_im;
    }
  return y;
}

/*
 * returns z = x^H y for vectors x, y of length L*L
 */
_Complex double
xdoty(size_t L, field *x, field *y)
{
  _Complex double z = 0; 
  for(int i=0; i<L*(L/VL); i++)
    for(int j=0; j<VL; j++) {
      float x_re = x[i].phi_re[j];
      float x_im = x[i].phi_im[j];
      float y_re = y[i].phi_re[j];
      float y_im = y[i].phi_im[j];
      float z_re = x_re*y_re + x_im*y_im;
      float z_im = x_re*y_im - x_im*y_re;
      z += z_re + _Complex_I*z_im;
    }
  return z;
}

/*
 * returns y = x - y for vectors y, x, of length L*L
 */
void
xmy(size_t L, field *x, field *y)
{
  for(int i=0; i<L*(L/VL); i++)
    for(int j=0; j<VL; j++) {
      y[i].phi_re[j] = x[i].phi_re[j] - y[i].phi_re[j];
      y[i].phi_im[j] = x[i].phi_im[j] - y[i].phi_im[j];
    }
  return;
}

/*
 * returns y = x for vectors y, x, of length L*L
 */
void
xeqy(size_t L, field *x, field *y)
{
  for(int i=0; i<L*(L/VL); i++)
    for(int j=0; j<VL; j++) {
      x[i].phi_re[j] = y[i].phi_re[j];
      x[i].phi_im[j] = y[i].phi_im[j];
    }
  return;
}

/*
 * returns y = a*x+y for vectors y, x, of length L*L and scalar a
 */
void
axpy(size_t L, _Complex float a, field *x, field *y)
{
  for(int i=0; i<L*(L/VL); i++)
    for(int j=0; j<VL; j++) {
      y[i].phi_re[j] = creal(a)*x[i].phi_re[j] - cimag(a)*x[i].phi_im[j] + y[i].phi_re[j];
      y[i].phi_im[j] = creal(a)*x[i].phi_im[j] + cimag(a)*x[i].phi_re[j] + y[i].phi_im[j];
    }
  return;
}

/*
 * returns y = x+a*y for vectors y, x, of length L*L and scalar a
 */
void
xpay(size_t L, field *x, _Complex float a, field *y)
{
  for(int i=0; i<L*(L/VL); i++)
    for(int j=0; j<VL; j++) {
      y[i].phi_re[j] = x[i].phi_re[j] + creal(a)*y[i].phi_re[j] - cimag(a)*y[i].phi_im[j];
      y[i].phi_im[j] = x[i].phi_im[j] + creal(a)*y[i].phi_im[j] + cimag(a)*y[i].phi_re[j];
    }
  return;
}

/*
 * Solves lapl(u) x = b, for x, given b, using Conjugate Gradient
 */
void
cg(size_t L, field *x, field *b, link *g)
{
  int max_iter = 100;
  float tol = 1e-9;

  /* Temporary fields needed for CG */
  field *r = new_field(L);
  field *p = new_field(L);
  field *Ap = new_field(L);

  /* Initial residual and p-vector */
  lapl(L, r, x, g);
  xmy(L, b, r);
  xeqy(L, p, r);

  /* Initial r-norm and b-norm */
  float rr = xdotx(L, r);  
  float bb = xdotx(L, b);
  double t_lapl = 0;
  int iter = 0;
  for(iter=0; iter<max_iter; iter++) {
    printf(" %6d, res = %+e\n", iter, rr/bb);
    if(sqrt(rr/bb) < tol)
      break;
    double t = stop_watch(0);
    lapl(L, Ap, p, g);
    t_lapl += stop_watch(t);
    float pAp = xdoty(L, p, Ap);
    float alpha = rr/pAp;
    axpy(L, alpha, p, x);
    axpy(L, -alpha, Ap, r);
    float r1r1 = xdotx(L, r);
    float beta = r1r1/rr;
    xpay(L, r, beta, p);
    rr = r1r1;
  }

  /* Recompute residual after convergence */
  lapl(L, r, x, g);
  xmy(L, b, r);
  rr = xdotx(L, r);

  /* 
   * _TODO_A_
   *
   * Compute the susstained flops and bandwidth. Note that t_lapl has
   * accumulated the time spent in the laplacian operator application.
   *
   */
  int N_fp_per_site = 34;
  int N_io_per_site = sizeof(float)*(2 + 4 + 2);
  double beta_fp = N_fp_per_site*((double)L*L)/(t_lapl/(double)iter)*1e-9;
  double beta_io = N_io_per_site*((double)L*L)/(t_lapl/(double)iter)*1e-9;
  printf(" Converged after %6d iterations, res = %+e\n", iter, rr/bb);  
  printf(" Time in lapl(): %+6.3e sec/call   %4.2e Gflop/s   %4.2e GB/s\n",
	 t_lapl/(double)iter, beta_fp, beta_io);  

  free(r);
  free(p);
  free(Ap);
  return;
}

/*
 * Usage info
 */
void
usage(char *argv[])
{
  fprintf(stderr, " Usage: %s L\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{
  if(argc != 2) {
    usage(argv);
    exit(1);
  }

  char *e;
  size_t L = (int)strtoul(argv[1], &e, 10);
  if(*e != '\0') {
    usage(argv);
    exit(2);
  }

  srand48(7);
  
  field *b = new_field(L);
  field *x = new_field(L);
  link *g = new_links(L);

  rand_links(L, g);
  rand_field(L, b);
  zero_field(L, x);

  cg(L, x, b, g);
  
  free(b);
  free(x);
  free(g);
  return 0;
}

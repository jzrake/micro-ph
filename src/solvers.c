


#include <stdio.h>
#include <math.h>

static int verbose = 0;

static const int    MAX_INTEGRAL_ITER    = 1000000;
static const double INTEGRATE_ACC_GOAL   = 1e-10;
static const double ZERO_SLOPE_REACHED   = 1e-12;
static const int    ZERO_SLOPE_REPEATED  = 10;
static const double ZERO_SECANT_REACHED  = 1e-6;
static const int    MAX_SECANT_ITER      = 500;



typedef double (*Dfunc)(double t, void *p);
double rootfind_secant(Dfunc f, Dfunc g, void *p);
double rootfind_newton(Dfunc f, Dfunc g, void *p);
double integrate_to_infinite(Dfunc f, void *p);
void plot_function(Dfunc f, void *p, double a, double b, const char *fname);


static double step_rk4 (Dfunc f, double t, void *p, double dt);
static double step_rk4a(Dfunc f, double t, void *p, double *h);


void solvers_set_verbose(int v)
{
  verbose = v;
}

void plot_function(Dfunc f, void *p, double a, double b, const char *fname)
{
  FILE *out = fopen(fname, "w");
  int i;

  for (i=0; i<1000; ++i) {
    double x = a + i*(b-a)/1000.0;
    fprintf(out, "%e %e\n", x, f(x,p));
  }

  fclose(out);
}

double rootfind_secant(Dfunc f, Dfunc g, void *p)
{
  int niter = 0;
  double x0 = -0.01; // starting guess values sandwich the root if possible
  double x1 = +0.01;
  double x2;

  while (f(x0, p) * f(x1, p) > 0.0) {
    x0 *= 2.0;
    x1 *= 2.0;
    if (verbose >= 1) {
      printf("[%s]: widening the initial bracket, [%f %f] -> [%f %f]\n",
	     __FUNCTION__, x0, x1, f(x0, p), f(x1, p));
    }
  }

  while (1) {

    double f0 = f(x0, p);
    double f1 = f(x1, p);

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x0 = x1;
    x1 = x2;

    if (fabs(f1) < ZERO_SECANT_REACHED) {
      break;
    }
    else if (++niter == MAX_SECANT_ITER) {
      printf("[%s]: warning! convergence took too many iterations\n",
	     __FUNCTION__);
      break;
    }
  }
  if (verbose >= 1) {
    printf("[%s]: converged to %e after %d iterations\n", __FUNCTION__,
	   x2, niter);
  }
  return x2;
}

double rootfind_newton(Dfunc f, Dfunc g, void *p)
{
  int niter = 0;
  double dx = 1e-8;
  double x = 0.0;

  while (1) {

    double f0 = f(x, p);
    double g0 = g ? g(x, p) : (f(x+0.5*dx, p) - f(x-0.5*dx, p))/dx;

    if (verbose >= 1) {
      printf("[%s]: f(%f) = %e\n", __FUNCTION__, x, f0);
    }

    x -= f0/g0;

    if (fabs(f0) < ZERO_SECANT_REACHED) {
      break;
    }
    else if (++niter == MAX_SECANT_ITER) {
      printf("[%s]: warning! convergence took too many iterations\n",
	     __FUNCTION__);
      break;
    }
  }
  if (verbose >= 1) {
    printf("[%s]: converged to %e after %d iterations\n", __FUNCTION__,
	   x, niter);
  }
  return x;
}



double integrate_to_infinite(Dfunc f, void *p)
// -----------------------------------------------------------------------------
// Calls the RK4 routine to integrate the function 'f' until its value is no
// longer changing. Typical algorithm parameters are ZERO_SLOPE_REACHED = 1e-10,
// ZERO_SLOPE_REPEATED = 10.
// -----------------------------------------------------------------------------
{
  int n_zero_slope = 0, niter = 0;
  double dx = 1e-3; // initial choice does not matter if adaptive step size
  double dy = dx;
  double x = 0.0;
  double y = 0.0;

  while (1) {
    
    if (verbose >= 2) {
      printf("[%s]: x=%f y(x)=%e\n", __FUNCTION__, x, y);
    }

    dy = step_rk4a(f, x, p, &dx);
    y += dy;
    x += dx;

    if (fabs(dy/dx) < ZERO_SLOPE_REACHED) {
      n_zero_slope += 1;
    }
    else {
      n_zero_slope = 0;
    }
    if (n_zero_slope == ZERO_SLOPE_REPEATED) {
      break;
    }
    else if (++niter == MAX_INTEGRAL_ITER) {
      printf("[%s]: warning! convergence took too many iterations\n",
	     __FUNCTION__);
      break;
    }
  }
  if (verbose >= 1) {
    printf("[%s]: converged to %e after %d iterations\n", __FUNCTION__,
	   y, niter);
  }
  return y;
}


double step_rk4(Dfunc f, double t, void *p, double dt)
// -----------------------------------------------------------------------------
// Take a 4th-order Runge-Kutta step for ODE integration, returns dy.
//
// http://en.wikipedia.org/wiki/Runge-Kutta_methods
//
// @inputs: (1) f(t,y) := y'(t)  function pointer
//          (2) t                time (or independent variable)
//          (3) y                solution value at time 't'
//          (4) dt               step size
//
// @return: dt * y'(t)
// -----------------------------------------------------------------------------
{
  const double k1 = dt * f(t, p);
  const double k2 = dt * f(t + 0.5*dt, p);
  const double k3 = dt * f(t + 0.5*dt, p);
  const double k4 = dt * f(t + 1.0*dt, p);

  return (1./6.)*(k1 + 2*k2 + 2*k3 + k4);
}

double esterr(Dfunc f, double t, void *p, double dt, double *dy)
{
  const double dy0     = step_rk4(f, t+0.0*dt, p, dt*0.5); // y(t + dt/2)
  const double dy_half = step_rk4(f, t+0.5*dt, p, dt*0.5) + dy0;
  const double dy_full = step_rk4(f, t, p, dt);

  if (dy != NULL) *dy = dy_half;

  if (fabs(dy_half) < 1e-16) {
    return fabs(dy_full - dy_half);
  }
  else {
    return fabs((dy_full - dy_half) / dy_half);
  }
}

double step_rk4a(Dfunc f, double t, void *p, double *h)
// -----------------------------------------------------------------------------
// Same as step_rk4, but uses step doubling as a means of adaptive step size
// control.
// -----------------------------------------------------------------------------
{
  static double dt = 1.0; // won't matter what the initial value is

  while (esterr(f,t,p,dt,NULL) < INTEGRATE_ACC_GOAL) {
    dt *= 2.0;
  }

  while (esterr(f,t,p,dt,NULL) > INTEGRATE_ACC_GOAL) {
    dt *= 0.5;
  }

  *h = dt;
  return step_rk4(f,t,p,dt);
}




#include <stdio.h>
#include <math.h>


static int    verbose = 0;
static int    MAX_INTEGRAL_ITER    = 1000000;
static double INTEGRATE_STEP_SIZE  = 1e-2;
static double ZERO_SLOPE_REACHED   = 1e-12;
static int    ZERO_SLOPE_REPEATED  = 10;
static double ZERO_SECANT_REACHED  = 1e-8;
static int    MAX_SECANT_ITER      = 500;

static double step_rk4(double (*f)(double t, double y), double t, double y, double dt);

void solvers_set_verbose(int v)
{
  verbose = v;
}

double rootfind_secant(double (*f)(double, void*), void *p)
{
  int niter = 0;
  double x0 = -0.01; // starting guess values sandwich the root if possible
  double x1 = +0.01;
  double x2;

  while (f(x0, p) * f(x1, p) > 0.0) {
    x0 *= 2.0;
    x1 *= 2.0;
    if (verbose > 0) {
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
  if (verbose > 0) {
    printf("[%s]: converged to %e after %d iterations\n", __FUNCTION__,
	   x2, niter);
  }
  return x2;
}

double rootfind_newton(double (*f)(double, void*), void *p)
{
  int niter = 0;
  double dx = 1e-10;
  double x = 0.0;

  while (1) {

    double f0 =  f(x, p);
    double g0 = (f(x+0.5*dx, p) - f(x-0.5*dx, p))/dx;

    if (verbose > 0) {
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
  if (verbose > 0) {
    printf("[%s]: converged to %e after %d iterations\n", __FUNCTION__,
	   x, niter);
  }
  return x;
}



double integrate_to_infinite(double (*f)(double t, double y))
// -----------------------------------------------------------------------------
// Calls the RK4 routine to integrate the function 'f' until its value is no
// longer changing. Typical algorithm parameters are ZERO_SLOPE_REACHED = 1e-10,
// ZERO_SLOPE_REPEATED = 10, and INTEGRATE_STEP_SIZE = 1e-2.
// -----------------------------------------------------------------------------
{
  int n_zero_slope = 0, niter = 0;
  double dx = INTEGRATE_STEP_SIZE;
  double dy = dx;
  double x = 0.0;
  double y = 0.0;

  while (1) {

    dy = step_rk4(f, x, y, dx);
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

  return y;
}




double step_rk4(double (*f)(double t, double y), double t, double y, double dt)
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
  const double k1 = dt * f(t, y);
  const double k2 = dt * f(t + 0.5*dt, y + 0.5*k1);
  const double k3 = dt * f(t + 0.5*dt, y + 0.5*k2);
  const double k4 = dt * f(t + 1.0*dt, y + 1.0*k3);

  return (1./6.)*(k1 + 2*k2 + 2*k3 + k4);
}




#include <stdio.h>
#include <math.h>

static int verbose = 0;

static const int    MAX_INTEGRAL_ITER    = 1000000;
static const double INTEGRATE_STEP_SIZE  = 1.0;
static const double INTEGRATE_ACC_GOAL   = 1e-10;
static const double ZERO_SLOPE_REACHED   = 1e-12;
static const int    ZERO_SLOPE_REPEATED  = 10;
static const double ZERO_SECANT_REACHED  = 1e-8;
static const int    MAX_SECANT_ITER      = 500;


typedef double (*Dfunc)(double t, double y); // derivative function

static double step_rk4 (Dfunc f, double t, double y, double dt);
static double step_rk4a(Dfunc f, double t, double y, double *h);


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



double integrate_to_infinite(Dfunc f)
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

    dy = step_rk4a(f, x, y, &dx);
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
  if (verbose > 0) {
    printf("[%s]: converged to %e after %d iterations\n", __FUNCTION__,
	   y, niter);
  }
  return y;
}




double step_rk4(Dfunc f, double t, double y, double dt)
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

double esterr(Dfunc f, double t, double y, double dt, double *dy)
{
  const double dy0 = step_rk4(f, t, y, dt*0.5); // y(t + dt/2)
  const double dy_half = step_rk4(f, t+0.5*dt, y+dy0, dt*0.5) + dy0;
  const double dy_full = step_rk4(f, t, y, dt);

  if (dy != NULL) *dy = dy_half;

  if (fabs(dy_half) < 1e-16) {
    return fabs(dy_full - dy_half);
  }
  else {
    return fabs(dy_full - dy_half) / dy_half;
  }
}

double step_rk4a(Dfunc f, double t, double y, double *h)
// -----------------------------------------------------------------------------
// Same as step_rk4, but uses step doubling as a means of adaptive step size
// control.
// -----------------------------------------------------------------------------
{
  static double dt = INTEGRATE_STEP_SIZE;

  while (esterr(f,t,y,dt,NULL) < INTEGRATE_ACC_GOAL) {
    dt *= 2.0;
  }

  while (esterr(f,t,y,dt,NULL) > INTEGRATE_ACC_GOAL) {
    dt *= 0.5;
  }

  *h = dt;
  return step_rk4(f, t, y, dt);
}

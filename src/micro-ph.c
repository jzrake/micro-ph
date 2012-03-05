


#include <stdio.h>
#include <math.h>


static double INTEGRATE_STEP_SIZE = 1e-2;
static double ZERO_SLOPE_REACHED  = 1e-10;
static int    ZERO_SLOPE_REPEATED = 10;
static int    MAXIMUM_ITERATION   = 1000000;


static double  EtaValue = 1.0; // ... mu/kT (degeneracy parameter)
static double BetaValue = 1.0; // ... kT


double step_rk4(double (*f)(double t, double y), double t, double y, double dt)
// -----------------------------------------------------------------------------
// Take a 4th-order Runge-Kutta step for ODE integration.
//
// http://en.wikipedia.org/wiki/Runge-Kutta_methods
//
// @inputs: (1) f(t,y) := y'(t)  function pointer
//          (2) t                time (or independent variable)
//          (3) y                solution value at time 't'
//          (4) dt               step size
//
// @return: y'(t)
// -----------------------------------------------------------------------------
{
  const double k1 = dt * f(t + 0.0*dt, y + 0.0);
  const double k2 = dt * f(t + 0.5*dt, y + 0.5*k1);
  const double k3 = dt * f(t + 0.5*dt, y + 0.5*k2);
  const double k4 = dt * f(t + 1.0*dt, y + 1.0*k3);

  return (1./6.)*(k1 + 2*k2 + 2*k3 + k4);
}


// -----------------------------------------------------------------------------
// These functions generate the integrand for the dimensionless number density
// of electrons (np) and positrons (ne) according to Sekiguchi (2010). The
// independent variable, x := pc. The answer is the number of particles in
// the volume pi^2 (hbar / me c)^3
// -----------------------------------------------------------------------------
double integrand_ne_42(double x, double t)
{
  return (x*x) / (exp(BetaValue * (sqrt(1 + x*x) - 1) - EtaValue) + 1);
}
double integrand_np_45(double x, double t)
{
  return (x*x) / (exp(BetaValue * (sqrt(1 + x*x) + 1) + EtaValue) + 1);
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
    else if (++niter == MAXIMUM_ITERATION) {
      printf("[%s]: warning! convergence took too many iterations\n", __FUNCTION__);
      break;
    }
  }

  return y;
}



double evaluate_ne_42(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_ne_42);
}
double evaluate_np_45(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_np_45);
}


void test_number_densities()
// -----------------------------------------------------------------------------
// Compare results of numerical integration with those giben by Mathematica's
// NIntegrate.
// -----------------------------------------------------------------------------
{
  printf("ne(1.0, 1.0) = %18.15e (8.518609867257567e+00)\n", evaluate_ne_42(1.0, 1.0));
  printf("np(1.0, 1.0) = %18.15e (2.176235567609514e-01)\n", evaluate_np_45(1.0, 1.0));
  printf("ne(1.0, 0.1) = %18.15e (4.696042780715451e+03)\n", evaluate_ne_42(1.0, 0.1));
  printf("np(1.0, 0.1) = %18.15e (6.390164655753489e+02)\n", evaluate_np_45(1.0, 0.1));
}


int main()
{
  test_number_densities();
  return 0;
}

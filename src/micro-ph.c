


#include <stdio.h>
#include <math.h>


static int    MAX_INTEGRAL_ITER    = 1000000;
static double INTEGRATE_STEP_SIZE  = 1e-2;
static double ZERO_SLOPE_REACHED   = 1e-10;
static int    ZERO_SLOPE_REPEATED  = 10;
static double ZERO_SECANT_REACHED  = 1e-12;
static int    MAX_SECANT_ITER      = 50;


static double  EtaValue = 1.0; // mu/kT     ... degeneracy parameter
static double BetaValue = 1.0; // me c^2/kT ... unitless inverse temperature


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


// -----------------------------------------------------------------------------
// These functions generate the integrand for the dimensionless number density,
// pressure, and energy density of electrons and positrons according to
// Sekiguchi (2010). The independent variable, x := pc.
// 
//                               *** Units ***
//
// Volume ... pi^2 (hbar / me c)^3
// Energy ... me c^2
// -----------------------------------------------------------------------------


// Number densities
// -----------------------------------------------------------------------------
double integrand_ne_42(double x, double t)
{
  return (x*x) / (exp(BetaValue * (sqrt(1 + x*x) - 1) - EtaValue) + 1);
}
double integrand_np_45(double x, double t)
{
  return (x*x) / (exp(BetaValue * (sqrt(1 + x*x) + 1) + EtaValue) + 1);
}

// Pressure
// -----------------------------------------------------------------------------
double integrand_Pe_43(double x, double t)
{
  const double numer = pow(x,4) / sqrt(1 + x*x);
  return numer / (exp(BetaValue * (sqrt(1 + x*x) - 1) - EtaValue) + 1);
}
double integrand_Pp_46(double x, double t)
{
  const double numer = pow(x,4) / sqrt(1 + x*x);
  return numer / (exp(BetaValue * (sqrt(1 + x*x) + 1) + EtaValue) + 1);
}

// Internal energy density
// -----------------------------------------------------------------------------
double integrand_ue_44(double x, double t)
{
  const double numer = pow(x,2) * (sqrt(1 + x*x) - 1);
  return numer / (exp(BetaValue * (sqrt(1 + x*x) - 1) - EtaValue) + 1);
}
double integrand_up_47(double x, double t)
{
  const double numer = pow(x,2) * (sqrt(1 + x*x) + 1);
  return numer / (exp(BetaValue * (sqrt(1 + x*x) + 1) + EtaValue) + 1);
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

double evaluate_Pe_43(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_Pe_43);
}
double evaluate_Pp_46(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_Pp_46);
}

double evaluate_ue_44(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_ue_44);
}
double evaluate_up_47(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_up_47);
}



double solve_for_eta(double beta, double C)
// -----------------------------------------------------------------------------
// Solves the implicit equation ne(e,b) - np(e,b) = C for e := eta, where C is a
// constant, typically the total number of positive charges in the
// characteristic volume V0 := pi^2 (hbar / me c)^3.
// -----------------------------------------------------------------------------
{
  int niter = 0;
  double eta0 = -1.0; // starting guess values sandwich the root if possible
  double eta1 = +1.0;
  double eta2;

  while (1) {

    double f0 = evaluate_ne_42(eta0, beta) - evaluate_np_45(eta0, beta) - C;
    double f1 = evaluate_ne_42(eta1, beta) - evaluate_np_45(eta1, beta) - C;

    eta2 = eta1 - f1 * (eta1 - eta0) / (f1 - f0);
    eta0 = eta1;
    eta1 = eta2;

    if (fabs(f1) < ZERO_SECANT_REACHED) {
      break;
    }
    else if (++niter == MAX_SECANT_ITER) {
      printf("[%s]: warning! convergence took too many iterations\n", __FUNCTION__);
      break;
    }
  }

  return eta2;
}


void microph_test_npu()
// -----------------------------------------------------------------------------
// Compare results of numerical integration with those given by Mathematica's
// NIntegrate. Integrates the equations for number density, pressure, and
// internal energy density.
// -----------------------------------------------------------------------------
{
  const char *sep = "**************************************************\n";

  printf("\ntesting number densities...\n");
  printf(sep);
  printf("ne(1.0, 1.0) = %18.15e (8.518609867257567e+00)\n", evaluate_ne_42(1.0, 1.0));
  printf("np(1.0, 1.0) = %18.15e (2.176235567609514e-01)\n", evaluate_np_45(1.0, 1.0));
  printf("ne(1.0, 0.1) = %18.15e (4.696042780715451e+03)\n", evaluate_ne_42(1.0, 0.1));
  printf("np(1.0, 0.1) = %18.15e (6.390164655753489e+02)\n", evaluate_np_45(1.0, 0.1));

  printf("\ntesting pressures...\n");
  printf(sep);
  printf("Pe(1.0, 1.0) = %18.15e (2.970525729398766e+01)\n", evaluate_Pe_43(1.0, 1.0));
  printf("Pp(1.0, 1.0) = %18.15e (6.562580252817735e-01)\n", evaluate_Pp_46(1.0, 1.0));
  printf("Pe(1.0, 0.1) = %18.15e (1.571436949436383e+05)\n", evaluate_Pe_43(1.0, 0.1));
  printf("Pp(1.0, 0.1) = %18.15e (1.953555533242760e+04)\n", evaluate_Pp_46(1.0, 0.1));

  printf("\ntesting internal energy densities...\n");
  printf(sep);
  printf("Pe(1.0, 1.0) = %18.15e (2.390933965258743e+01)\n", evaluate_ue_44(1.0, 1.0));
  printf("Pp(1.0, 1.0) = %18.15e (9.540921687756961e-01)\n", evaluate_up_47(1.0, 1.0));
  printf("Pe(1.0, 0.1) = %18.15e (1.526403704611330e+05)\n", evaluate_ue_44(1.0, 0.1));
  printf("Pp(1.0, 0.1) = %18.15e (2.020504623818944e+04)\n", evaluate_up_47(1.0, 0.1));

}

void microph_test_eta()
{
  const char *sep = "**************************************************\n";

  printf("\ntesting solution to chemical potential...\n");
  printf(sep);
  printf("%+18.15e (%+18.15e)\n", solve_for_eta(1e0, 1.0), -0.6525037133686798);
  printf("%+18.15e (%+18.15e)\n", solve_for_eta(1e1, 1.0),  7.316055681629137);
  printf("%+18.15e (%+18.15e)\n", solve_for_eta(1e2, 1.0), 75.47842301516384);
}



int main()
{
  microph_test_npu();
  microph_test_eta();

  return 0;
}

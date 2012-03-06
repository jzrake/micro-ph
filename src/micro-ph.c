

/*------------------------------------------------------------------------------
 * FILE: micro-ph.c
 *
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 *
 * DESCRIPTION: Calculates contributions to the EOS from photons, pairs, and
 * neutrinos.
 *
 *
 * REFERENCES:
 * 
 * (+) Sekiguchi (2010) http://arxiv.org/abs/1009.3320
 *
 *
 * NOTES:
 * 
 * (+) Throughout the code comments, the letter h means 'h bar'
 *
 *------------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>



const int mphElectrons   =  (1 << 0);
const int mphPositrons   =  (1 << 1);
const int mphNeutrinos   =  (1 << 2);
const int mphPhotons     =  (1 << 3);


static const double HBAR_C           = 1.973269718e+02; // MeV-fm
static const double ELECTRON_MASS    = 5.110998928e-01; // MeV
static const double ATOMIC_MASS_UNIT = 9.314940612e+02; // MeV

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
// Sekiguchi (2010). The independent variable, x := p/mc^2.
// 
//                               *** Units ***
//
// Volume ... pi^2 (h/mc)^3
// Energy ... m c^2
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
double integrand_pe_43(double x, double t)
{
  const double numer = pow(x,4) / sqrt(1 + x*x);
  return numer / (exp(BetaValue * (sqrt(1 + x*x) - 1) - EtaValue) + 1);
}
double integrand_pp_46(double x, double t)
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


// -----------------------------------------------------------------------------
// These functions generate the integrand for the dimensionless number density,
// pressure, and energy density of electrons and positrons assuming neutrinos
// have zero mass. The independent variable, x := pc/kT.
// 
//                               *** Units ***
//
// Volume ... pi^2 (hc/kT)^3
// Energy ... kT
// -----------------------------------------------------------------------------
double integrand_neutrino_n(double x, double t)
{
  return (x*x) / (exp(x - EtaValue) + 1);
}
double integrand_neutrino_p(double x, double t)
{
  return (x*x*x) / (exp(x - EtaValue) + 1);
}
double integrand_neutrino_u(double x, double t)
{
  return (x*x*x) / (exp(x - EtaValue) + 1);
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



double evaluate_ne(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_ne_42);
}
double evaluate_np(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_np_45);
}

double evaluate_pe(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_pe_43);
}
double evaluate_pp(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_pp_46);
}
double evaluate_ue(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_ue_44);
}
double evaluate_up(double eta, double beta)
{
  EtaValue = eta;
  BetaValue = beta;
  return integrate_to_infinite(integrand_up_47);
}
double evaluate_neutrino_n(double eta)
{
  EtaValue = eta;
  return integrate_to_infinite(integrand_neutrino_n);
}
double evaluate_neutrino_p(double eta)
{
  EtaValue = eta;
  return integrate_to_infinite(integrand_neutrino_p);
}
double evaluate_neutrino_u(double eta)
{
  EtaValue = eta;
  return integrate_to_infinite(integrand_neutrino_u);
}



double solve_for_eta_pairs(double beta, double C)
// -----------------------------------------------------------------------------
// Solves the implicit equation ne(e,b) - np(e,b) = C for e := eta, where C is a
// constant, typically the total number of positive charges in the
// characteristic volume V0 := pi^2 (h/mc)^3.
// -----------------------------------------------------------------------------
{
  int niter = 0;
  double eta0 = -1.0; // starting guess values sandwich the root if possible
  double eta1 = +1.0;
  double eta2;

  while (1) {

    double f0 = evaluate_ne(eta0, beta) - evaluate_np(eta0, beta) - C;
    double f1 = evaluate_ne(eta1, beta) - evaluate_np(eta1, beta) - C;

    eta2 = eta1 - f1 * (eta1 - eta0) / (f1 - f0);
    eta0 = eta1;
    eta1 = eta2;

    if (fabs(f1) < ZERO_SECANT_REACHED) {
      break;
    }
    else if (++niter == MAX_SECANT_ITER) {
      printf("[%s]: warning! convergence took too many iterations\n",
	     __FUNCTION__);
      break;
    }
  }

  return eta2;
}

double solve_for_eta_neutrino(double beta, double C)
// -----------------------------------------------------------------------------
// Solves the implicit equation n(e,b) = C for e := eta, where C is a constant,
// typically the total neutrino number of a particular species in the
// characteristic volume V0 := pi^2 (h/mc)^3.
// -----------------------------------------------------------------------------
{
  int niter = 0;
  double eta0 = -1.0; // starting guess values sandwich the root if possible
  double eta1 = +1.0;
  double eta2;

  while (1) {

    double f0 = evaluate_ne(eta0, beta) - C; // equation S2010-4.9
    double f1 = evaluate_ne(eta1, beta) - C;

    eta2 = eta1 - f1 * (eta1 - eta0) / (f1 - f0);
    eta0 = eta1;
    eta1 = eta2;

    if (fabs(f1) < ZERO_SECANT_REACHED) {
      break;
    }
    else if (++niter == MAX_SECANT_ITER) {
      printf("[%s]: warning! convergence took too many iterations\n",
	     __FUNCTION__);
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
  const double T0[] = { 8.518609867257567, 0.21762355676095138, // true
			4696.042780715451, 639.0164655753489 };
  printf("ne(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_ne(1.0, 1.0), T0[0]);
  printf("np(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_np(1.0, 1.0), T0[1]);
  printf("ne(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_ne(1.0, 0.1), T0[2]);
  printf("np(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_np(1.0, 0.1), T0[3]);

  printf("\ntesting pressures...\n");
  printf(sep);
  const double T1[] = { 29.705257293987664, 0.6562580252817735,
			157143.69494363826, 19535.555332427597 };
  printf("pe(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_pe(1.0, 1.0), T1[0]);
  printf("pp(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_pp(1.0, 1.0), T1[1]);
  printf("pe(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_pe(1.0, 0.1), T1[2]);
  printf("pp(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_pp(1.0, 0.1), T1[3]);

  printf("\ntesting internal energy densities...\n");
  printf(sep);
  const double T2[] = { 23.909339652587427, 0.9540921687756961,
			152640.37046113302, 20205.04623818944 };
  printf("pe(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_ue(1.0, 1.0), T2[0]);
  printf("pp(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_up(1.0, 1.0), T2[1]);
  printf("pe(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_ue(1.0, 0.1), T2[2]);
  printf("pp(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_up(1.0, 0.1), T2[3]);
}


void microph_test_eta()
{
  const char *sep = "**************************************************\n";
  double (*f)(double, double);

  f = solve_for_eta_pairs;
  printf("\ntesting solution to chemical potential of pairs...\n");
  printf(sep);
  printf("%+18.15e (%+18.15e)\n", f(1e0, 1.0), -0.6525037133686798);
  printf("%+18.15e (%+18.15e)\n", f(1e1, 1.0),  7.316055681629137);
  printf("%+18.15e (%+18.15e)\n", f(1e2, 1.0), 75.47842301516384);

  f = solve_for_eta_neutrino;
  printf("\ntesting solution to chemical potential of neutrinos...\n");
  printf(sep);
  printf("%+18.15e\n", f(1e-1, 1.0));
  printf("%+18.15e\n", f(1e+0, 1.0));
  printf("%+18.15e\n", f(1e+1, 1.0));
}


double microph_get(int flags)
// -----------------------------------------------------------------------------
// Evaluate the pressure (in MeV/fm^3), including terms specified by 'flags'.
//
// Implemented terms are:
//
// (+) mphElectrons
// (+) mphPositrons
// (+) mphNeutrinos
// (+) mphPhotons
//
// -----------------------------------------------------------------------------
{
  const double kT = 1.0; // MeV
  const double eta = 1.0;

  double p = 0.0;

  // Electrons
  // ---------------------------------------------------------------------------
  if (flags & mphElectrons) {
    const double Volume = pow(M_PI, 2) * pow(HBAR_C/ELECTRON_MASS, 3);
    const double Energy = ELECTRON_MASS;
    const double beta = ELECTRON_MASS / kT;

    p += (Energy / Volume) * evaluate_pe(eta, beta);
  }

  // Positrons
  // ---------------------------------------------------------------------------
  if (flags & mphPositrons) {

    const double Volume = pow(M_PI, 2) * pow(HBAR_C/ELECTRON_MASS, 3);
    const double Energy = ELECTRON_MASS;
    const double beta = ELECTRON_MASS / kT;

    p += (Energy / Volume) * evaluate_pp(eta, beta);
  }

  // Neutrinos
  // ---------------------------------------------------------------------------
  if (flags & mphNeutrinos) {
    const double Volume = pow(M_PI, 2) * pow(HBAR_C/kT, 3);
    const double Energy = kT;

    p += (Energy / Volume) * evaluate_neutrino_p(eta);
  }

  // Photons
  // ---------------------------------------------------------------------------
  if (flags & mphPhotons) {

  }

  return p;
}


int main()
{
  microph_test_npu();
  microph_test_eta();

  microph_get(mphPositrons | mphPhotons | mphNeutrinos);

  return 0;
}

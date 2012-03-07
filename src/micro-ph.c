

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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>


// External function declarations
// -----------------------------------------------------------------------------
double rootfind_secant(double (*f)(double, void*), void *p);
double rootfind_newton(double (*f)(double, void*), void *p);

double integrate_to_infinite(double (*f)(double t, double y));
void solvers_set_verbose(int v);

const int mphElectrons   =  (1 << 0); // option flags
const int mphPositrons   =  (1 << 1);
const int mphNeutrinos   =  (1 << 2);
const int mphPhotons     =  (1 << 3);

struct ThermalState
{
  double density;  // baryon mass density (g/cm^3)
  double kT;       // temperature (MeV)
  double mu_ep;    // chemical potential for e+/e- pairs (MeV)
  double mu_nu;    // chemical potential for neutrino (MeV)
  double Ye;       // electron/proton fraction
  double n;        // number density (1/fm^3)
  double p;        // pressure (MeV/fm^3)
  double u;        // internal energy (MeV/fm^3)
};


static double (*rootfind)(double (*f)(double, void*), void *p);

static const double LIGHT_SPEED      = 2.997924580e+10; // cm/s
static const double HBAR_C           = 1.973269718e+02; // MeV-fm
static const double ELECTRON_MASS    = 5.110998928e-01; // MeV
static const double ATOMIC_MASS_UNIT = 9.314940612e+02; // MeV
static const double MEV_TO_ERG       = 1.602176487e-06;
static const double FM3_TO_CM3       = 1.000000000e-39;

static double  EtaValue = 1.0; // mu/kT     ... degeneracy parameter
static double BetaValue = 1.0; // me c^2/kT ... unitless inverse temperature




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




double evaluate_ne(double eta, double beta)
{
  printf("[%s]: eta = %f, beta = %f\n", __FUNCTION__, eta, beta);
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


double relation_eta_pairs(double eta, void *p)
// -----------------------------------------------------------------------------
// Defines the implicit equation ne(e,b) - np(e,b) = C for e := eta, where C is
// a dimensionless constant, typically the total number of positive charges in
// the characteristic volume V0 := pi^2 (h/mc)^3.
// -----------------------------------------------------------------------------
{
  const double beta = ((double*)p)[0];
  const double C    = ((double*)p)[1];
  return evaluate_ne(eta, beta) - evaluate_np(eta, beta) - C;
}
double solve_for_eta_pairs(double beta, double C)
{
  double p[2] = { beta, C };
  return rootfind(relation_eta_pairs, p);
}

double relation_eta_neutrino(double eta, void *p)
// -----------------------------------------------------------------------------
// Defines the implicit equation n(e,b) = C for e := eta, where C is a constant,
// typically the total neutrino number of a particular species in the
// characteristic volume V0 := pi^2 (h/mc)^3.
// -----------------------------------------------------------------------------
{
  const double beta = ((double*)p)[0];
  const double C    = ((double*)p)[1];
  return evaluate_ne(eta, beta) - C; // equation S2010-4.9
}
double solve_for_eta_neutrino(double beta, double C)
{
  double p[2] = { beta, C };
  return rootfind(relation_eta_neutrino, p);
}


void microph_get_chemical_potential_pairs(struct ThermalState *S)
{
  const double c2 = LIGHT_SPEED*LIGHT_SPEED;
  const double Vol = pow(M_PI, 2) * pow(HBAR_C/ELECTRON_MASS, 3);
  const double Erest = S->density * c2 * FM3_TO_CM3 / MEV_TO_ERG;
  const double N = Vol * S->Ye * Erest / ATOMIC_MASS_UNIT;
  const double beta = ELECTRON_MASS / S->kT;

  printf("[%s]: beta = %f N = %f\n", __FUNCTION__, beta, N);

  double param[2] = { beta, N };
  const double eta = rootfind(relation_eta_pairs, param);
  S->mu_ep = eta * S->kT;
}


void microph_get(struct ThermalState *S, int flags)
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
  S->n = 0.0;
  S->p = 0.0;
  S->u = 0.0;

  // Electrons
  // ---------------------------------------------------------------------------
  if (flags & mphElectrons) {
    const double Volume = pow(M_PI, 2) * pow(HBAR_C/ELECTRON_MASS, 3);
    const double Energy = ELECTRON_MASS;
    const double eta = S->mu_ep / S->kT;
    const double beta = ELECTRON_MASS / S->kT;

    S->n += (1.0    / Volume) * evaluate_ne(eta, beta);
    S->p += (Energy / Volume) * evaluate_pe(eta, beta);
    S->u += (Energy / Volume) * evaluate_ue(eta, beta);
  }

  // Positrons
  // ---------------------------------------------------------------------------
  if (flags & mphPositrons) {
    const double Volume = pow(M_PI, 2) * pow(HBAR_C/ELECTRON_MASS, 3);
    const double Energy = ELECTRON_MASS;
    const double eta = S->mu_ep / S->kT;
    const double beta = ELECTRON_MASS / S->kT;

    S->n += (1.0    / Volume) * evaluate_np(eta, beta);
    S->p += (Energy / Volume) * evaluate_pp(eta, beta);
    S->u += (Energy / Volume) * evaluate_up(eta, beta);
  }

  // Neutrinos
  // ---------------------------------------------------------------------------
  if (flags & mphNeutrinos) {
    const double Volume = pow(M_PI, 2) * pow(HBAR_C/S->kT, 3);
    const double Energy = S->kT;
    const double eta = S->mu_nu / S->kT;

    S->n += (Energy / Volume) * evaluate_neutrino_n(eta);
    S->p += (Energy / Volume) * evaluate_neutrino_p(eta);
    S->u += (Energy / Volume) * evaluate_neutrino_u(eta);
  }

  // Photons
  // ---------------------------------------------------------------------------
  if (flags & mphPhotons) {

  }
}



void microph_test_npu()
// -----------------------------------------------------------------------------
// Compare results of numerical integration with those given by Mathematica's
// NIntegrate. Integrates the equations for number density, pressure, and
// internal energy density.
// -----------------------------------------------------------------------------
{
  const char *sep = "**************************************************\n";

  printf("\ntesting number densities\n");
  printf(sep);
  const double T0[] = { 8.518609867257567, 0.21762355676095138, // true
                        4696.042780715451, 639.0164655753489 };
  printf("ne(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_ne(1.0, 1.0), T0[0]);
  printf("np(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_np(1.0, 1.0), T0[1]);
  printf("ne(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_ne(1.0, 0.1), T0[2]);
  printf("np(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_np(1.0, 0.1), T0[3]);

  printf("\ntesting pressures\n");
  printf(sep);
  const double T1[] = { 29.705257293987664, 0.6562580252817735,
                        157143.69494363826, 19535.555332427597 };
  printf("pe(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_pe(1.0, 1.0), T1[0]);
  printf("pp(1.0, 1.0) = %18.15e (%18.15e)\n", evaluate_pp(1.0, 1.0), T1[1]);
  printf("pe(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_pe(1.0, 0.1), T1[2]);
  printf("pp(1.0, 0.1) = %18.15e (%18.15e)\n", evaluate_pp(1.0, 0.1), T1[3]);

  printf("\ntesting internal energy densities\n");
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
  printf("\ntesting solution to chemical potential of pairs\n");
  printf(sep);
  printf("%+18.15e (%+18.15e)\n", f(1e0, 1.0), -0.6525037133686798);
  printf("%+18.15e (%+18.15e)\n", f(1e1, 1.0),  7.316055681629137);
  printf("%+18.15e (%+18.15e)\n", f(1e2, 1.0), 75.47842301516384);

  printf("\ntesting the rootfinder a bit harder\n");
  printf(sep);
  {
    double beta = 0.010;
    double C0 = 10.0;
    double eta = solve_for_eta_pairs(beta, C0);
    double C1 = evaluate_ne(eta, beta) - evaluate_np(eta, beta);
    printf("error = %e\n", (C0 - C1)/C0);
  }
  {
    double beta = 0.10;
    double C0 = 0.01;
    double eta = solve_for_eta_pairs(beta, C0);
    double C1 = evaluate_ne(eta, beta) - evaluate_np(eta, beta);
    printf("error = %e\n", (C0 - C1)/C0);
  }

  f = solve_for_eta_neutrino;
  printf("\ntesting solution to chemical potential of neutrinos\n");
  printf(sep);
  printf("%+18.15e\n", f(1e-1, 1.0));
  printf("%+18.15e\n", f(1e+0, 1.0));
  printf("%+18.15e\n", f(1e+1, 1.0));
}


void microph_test_eos()
{
  const char *sep = "**************************************************\n";
  printf("\ntesting pressure contributions from various terms\n");
  printf(sep);

  struct ThermalState *S = (struct ThermalState*)
    malloc(sizeof(struct ThermalState));

  S->kT = 40.0;
  S->Ye = 0.08;
  S->density = 1e13;

  microph_get_chemical_potential_pairs(S);
  printf("Chemical potential of pairs = %e MeV (eta = %f)\n",
         S->mu_ep, S->mu_ep/S->kT);

  microph_get(S, mphElectrons);
  printf("Electron pressure = %e MeV/fm^3\n", S->p);

  microph_get(S, mphPositrons);
  printf("Positron pressure = %e MeV/fm^3\n", S->p);

  free(S);
}




#include <ctype.h>
#include <getopt.h>

int main(int argc, char **argv)
{
  rootfind = rootfind_newton;

  int c;
  while ((c = getopt(argc, argv, "vh")) != -1) {
    switch (c) {
    case 'v':
      solvers_set_verbose(1);
      break;
    case 'h':
      printf("usage: micro-ph [-v]\n");
      return 1;
    case '?':
      printf("usage: micro-ph [-v]\n");
      return 1;
    default:
      abort();
    }
  }

  microph_test_npu();
  //  microph_test_eta();
  //  microph_test_eos();

  return 0;
}

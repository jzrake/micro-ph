

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
typedef double (*Dfunc)(double t, void *p);
double rootfind_secant(Dfunc f, Dfunc g, void *p);
double rootfind_newton(Dfunc f, Dfunc g, void *p);
double integrate_to_infinite(Dfunc f, void *p);
void plot_function(Dfunc f, void *p, double a, double b, const char *fname);

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


static double (*rootfind)(Dfunc f, Dfunc g, void *p);

static const double LIGHT_SPEED      = 2.997924580e+10; // cm/s
static const double HBAR_C           = 1.973269718e+02; // MeV-fm
static const double ELECTRON_MASS    = 5.110998928e-01; // MeV
static const double ATOMIC_MASS_UNIT = 9.314940612e+02; // MeV
static const double MEV_TO_ERG       = 1.602176487e-06;
static const double FM3_TO_CM3       = 1.000000000e-39;


double pdf_fermion(double x, void *p)
// -----------------------------------------------------------------------------
// x     := pc/mc^2   ... independent variable
// eta   := mu/kT     ... degeneracy parameter
// beta  := mc^2/kT   ... unitless inverse temperature
//
//                    *** Units ***
//
// Volume             ... pi^2 (h/mc)^3
// Energy             ... m c^2
// -----------------------------------------------------------------------------
{
  const double sgn  = ((double*)p)[0]; // (-) for anti-fermions
  const double eta  = ((double*)p)[1];
  const double beta = ((double*)p)[2];
  return (x*x) / (exp(beta * (sqrt(1 + x*x) - sgn) - sgn*eta) + 1);
}
double pdf_fermion_deriv(double x, void *p)
// -----------------------------------------------------------------------------
// The derivative with respect to chemical potential.
// -----------------------------------------------------------------------------
{
  const double sgn  = ((double*)p)[0]; // (-) for anti-fermions
  const double eta  = ((double*)p)[1];
  const double beta = ((double*)p)[2];
  return 0.5 * sgn * (x*x) / (1 + cosh(beta*sqrt(1 + x*x) - sgn*(beta + eta)));
}
double ferm_number_density(double x, void *p)  { return 1.0; }
double ferm_pressure(double x, void *p) { return pow(x,2) / (sqrt(1 + x*x)); }
double ferm_internal_energy(double x, void *p)
{
  const double sgn = ((double*)p)[0];
  return sqrt(1 + x*x) - sgn*1;
}


double pdf_fermion_massless(double x, void *p)
// -----------------------------------------------------------------------------
// x     := pc/kT     ... independent variable
// eta   := mu/kT     ... degeneracy parameter
//
//                    *** Units ***
//
// Volume             ... pi^2 (hc/kT)^3
// Energy             ... kT
// -----------------------------------------------------------------------------
{
  const double sgn  = ((double*)p)[0]; // (-) for anti-fermions
  const double eta  = ((double*)p)[1];
  return (x*x) / (exp(x - sgn*eta) + 1);
}
double pdf_fermion_massless_deriv(double x, void *p)
// -----------------------------------------------------------------------------
// The derivative with respect to chemical potential.
// -----------------------------------------------------------------------------
{
  const double sgn  = ((double*)p)[0]; // (-) for anti-fermions
  const double eta  = ((double*)p)[1];
  return 0.5 * sgn * (x*x) / (1 + cosh(x - sgn*eta));
}
double fmml_number_density(double x, void *p)  { return 1.0; }
double fmml_pressure(double x, void *p)        { return x; }
double fmml_internal_energy(double x, void *p) { return x; }



double integrand_n(double x, void *p)
{
  return pdf_fermion(x,p) * ferm_number_density(x,p);
}
double integrand_p(double x, void *p)
{
  return pdf_fermion(x,p) * ferm_pressure(x,p);
}
double integrand_u(double x, void *p)
{
  return pdf_fermion(x,p) * ferm_internal_energy(x,p);
}

double integrand_n_massless(double x, void *p)
{
  return pdf_fermion_massless(x,p) * ferm_number_density(x,p);
}
double integrand_p_massless(double x, void *p)
{
  return pdf_fermion_massless(x,p) * ferm_pressure(x,p);
}
double integrand_u_massless(double x, void *p)
{
  return pdf_fermion_massless(x,p) * ferm_internal_energy(x,p);
}


double dntegrand_n(double x, void *p)
{
  return pdf_fermion_deriv(x,p) * ferm_number_density(x,p);
}
double dvaluate_ne(double eta, double beta)
{
  double p[3] = { +1, eta, beta };
  return integrate_to_infinite(dntegrand_n, p);
}
double dvaluate_np(double eta, double beta)
{
  double p[3] = { -1, eta, beta };
  return integrate_to_infinite(dntegrand_n, p);
}


double evaluate_ne(double eta, double beta)
{
  double p[3] = { +1, eta, beta };
  return integrate_to_infinite(integrand_n, p);
}
double evaluate_np(double eta, double beta)
{
  double p[3] = { -1, eta, beta };
  return integrate_to_infinite(integrand_n, p);
}

double evaluate_pe(double eta, double beta)
{
  double p[3] = { +1, eta, beta };
  return integrate_to_infinite(integrand_p, p);
}
double evaluate_pp(double eta, double beta)
{
  double p[3] = { -1, eta, beta };
  return integrate_to_infinite(integrand_p, p);
}

double evaluate_ue(double eta, double beta)
{
  double p[3] = { +1, eta, beta };
  return integrate_to_infinite(integrand_u, p);
}
double evaluate_up(double eta, double beta)
{
  double p[3] = { -1, eta, beta };
  return integrate_to_infinite(integrand_u, p);
}


double evaluate_neutrino_n(double eta)
{
  double p[2] = { +1, eta };
  return integrate_to_infinite(integrand_n_massless, p);
}
double evaluate_neutrino_p(double eta)
{
  double p[2] = { +1, eta };
  return integrate_to_infinite(integrand_p_massless, p);
}
double evaluate_neutrino_u(double eta)
{
  double p[2] = { +1, eta };
  return integrate_to_infinite(integrand_u_massless, p);
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
double delation_eta_pairs(double eta, void *p) // its derivative
{
  const double beta = ((double*)p)[0];
  return dvaluate_ne(eta, beta) - dvaluate_np(eta, beta);
}

double solve_for_eta_pairs(double beta, double C)
{
  double p[2] = { beta, C };
  //  return rootfind(relation_eta_pairs, delation_eta_pairs, p);
  return rootfind(relation_eta_pairs, NULL, p);
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
  return rootfind(relation_eta_neutrino, NULL, p);
}


void microph_get_chemical_potential_pairs(struct ThermalState *S)
// -----------------------------------------------------------------------------
// This is just the dimensional version of solve_for_eta_pairs, it returns the
// chemical potential in MeV.
// -----------------------------------------------------------------------------
{
  const double c2 = LIGHT_SPEED*LIGHT_SPEED;
  const double Vol = pow(M_PI, 2) * pow(HBAR_C/ELECTRON_MASS, 3);
  const double Erest = S->density * c2 * FM3_TO_CM3 / MEV_TO_ERG;
  const double C = Vol * S->Ye * Erest / ATOMIC_MASS_UNIT;
  const double beta = ELECTRON_MASS / S->kT;
  const double eta = solve_for_eta_pairs(beta, C);
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

  /*
  double param_e[] = { +1,1,1 };
  double param_p[] = { -1,1,1 };
  plot_function(pdf_fermion_deriv, param_e, 0.0, 10.0, "electron.dat");
  plot_function(pdf_fermion_deriv, param_p, 0.0, 10.0, "positron.dat");
  */

  printf("\ntesting derivative of fermi integral\n");
  printf(sep);
  const double T3[] = { 6.574627294589577, -0.21539712325750585 };
  printf("dne(1.0, 1.0) = %18.15e (%18.15e)\n", dvaluate_ne(1.0, 1.0), T3[0]);
  printf("dnp(1.0, 1.0) = %18.15e (%18.15e)\n", dvaluate_np(1.0, 1.0), T3[1]);
}


void microph_test_eta()
{
  const char *sep = "**************************************************\n";
  double (*f)(double, double);

  f = solve_for_eta_pairs;
  printf("\ntesting solution to chemical potential of pairs\n");
  printf(sep);
  printf("%+18.15e (%+18.15e)\n", f(1e0, 1.0), -0.6525037133686798);
  //  printf("%+18.15e (%+18.15e)\n", f(1e1, 1.0),  7.316055681629137);
  //  printf("%+18.15e (%+18.15e)\n", f(1e2, 1.0), 75.47842301516384);

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
  while ((c = getopt(argc, argv, "v:h")) != -1) {
    switch (c) {
    case 'v':
      printf("using verbose level %s\n", optarg);
      solvers_set_verbose(atoi(optarg));
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

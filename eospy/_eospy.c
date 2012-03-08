
/*------------------------------------------------------------------------------
 * FILE: _eospy.c
 *
 *
 * AUTHORS: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *          Dan Foreman-Mackey, NYU CCPP: danfm@nyu.edu
 *
 *
 * DESCRIPTION: Calculates contributions to the EOS from photons, pairs, and
 * neutrinos.
 *
 *
 *------------------------------------------------------------------------------
 */


// External function declarations
// -----------------------------------------------------------------------------
typedef double (*Dfunc)(double t, void *p);
double rootfind_secant(Dfunc f, Dfunc g, void *p);
double rootfind_newton(Dfunc f, Dfunc g, void *p);
double integrate_to_infinite(Dfunc f, void *p);
void plot_function(Dfunc f, void *p, double a, double b, const char *fname);
void microph_set_verbose(int v);

#define mphElectrons     (1 << 0) // option flags
#define mphPositrons     (1 << 1)
#define mphNeutrinos     (1 << 2)
#define mphPhotons       (1 << 3)

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



#include <Python.h>


static char doc[] = "";
static char eospy_doc[] = "";

PyMODINIT_FUNC init_eospy(void);
static PyObject *eospy_set_verbose(PyObject *self, PyObject *args);
static PyObject *eospy_pressure(PyObject *self, PyObject *args, PyObject *kwargs);

static PyMethodDef module_methods[] = {
  {"set_verbose", (PyCFunction)eospy_set_verbose, METH_VARARGS, NULL},
  {"pressure", (PyCFunction)eospy_pressure, METH_VARARGS|METH_KEYWORDS, NULL},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_eospy(void)
{
  PyObject *m = Py_InitModule3("_eospy", module_methods, doc);
  if (m == NULL) return;
}

PyObject *eospy_set_verbose(PyObject *self, PyObject *args)
{
  int v=0;

  if (!PyArg_ParseTuple(args, "i", &v)) {
    return NULL;
  }
  microph_set_verbose(v);

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *eospy_pressure(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int photons=0, electrons=0, positrons=0;
  static char *kwlist[] = {"density", "temperature", "proton_fraction",
			   "photons", "electrons", "positrons", NULL};

  struct ThermalState S;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddd|iii", kwlist,
				   &S.density, &S.kT, &S.Ye,
				   &photons, &electrons, &positrons)) {
    return NULL;
  }

  int flags = 0;

  if (photons  ) flags |= mphPhotons;
  if (electrons) flags |= mphElectrons;
  if (positrons) flags |= mphPositrons;

  microph_get_chemical_potential_pairs(&S);
  microph_get(&S, flags);

  PyObject *ret;

  ret = Py_BuildValue("d", S.p);

  if (ret == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Couldn't build output tuple.");
    return NULL;
  }

  return ret;
}

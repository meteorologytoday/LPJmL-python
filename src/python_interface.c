/**************************************************************************************/
/**                                                                                \n**/
/**     p y t h o n _ i n t e r f a c e . c                                        \n**/
/**                                                                                \n**/
/**     NumPy C-API wrapper exposing driver.c as the Python module                 \n**/
/**     "LPJmLCoupler".                                                             \n**/
/**                                                                                \n**/
/**     Three functions:                                                            \n**/
/**       lpjml_init(config_filename, dt_fast=3600) -> dict                        \n**/
/**       lpjml_update(year, month, day, hour, minute, seconds,                    \n**/
/**                    prec_in, temp_mean_in, swdown_in, lwnet_in,                 \n**/
/**                    t_flux_in, qflux_in, dedq_in, dhdt_in, drdt_in,             \n**/
/**                    drag_q_in, p_surf_in, wind_in, co2_in,                      \n**/
/**                    land_frac_in, q_ca_old_in,                                  \n**/
/**                    carbon_flux_out, q_ca_out, runoff_out,                      \n**/
/**                    roughness_length_out, surface_temp_out,                     \n**/
/**                    albedo_out, evap1_out, dedt_out, gc_out, exlpj_out)         \n**/
/**       lpjml_end() -> None                                                      \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file   \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3              \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/              \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                   \n**/
/**                                                                                \n**/
/**************************************************************************************/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL LPJmLCoupler_ARRAY_API

#include <Python.h>
#include <numpy/arrayobject.h>

#ifdef USE_MPI
#  include <mpi.h>
#endif

/*---------------------------------------------------------------------------*/
/* Forward declarations from driver.c                                        */
/*---------------------------------------------------------------------------*/

void lpj_init_(const int  *dt_fast_in,
               const char *config_filename,
               int        *global_number_lpj_cells,
               int        *firstyear_out,
               int        *lastyear_out,
               int        *nspinup_out);

void lpj_update_(const int    *Time_year,
                 const int    *Time_month,
                 const int    *Time_day,
                 const int    *Time_hour,
                 const int    *Time_minute,
                 const int    *Time_seconds,
                 const double *prec_in,
                 const double *temp_mean_in,
                 const double *swdown_in,
                 const double *lwnet_in,
                 const double *t_flux_in,
                 const double *qflux_in,
                 const double *dedq_in,
                 const double *dhdt_in,
                 const double *drdt_in,
                 const double *drag_q_in,
                 const double *p_surf_in,
                 const double *wind_in,
                 const double *co2_in,
                 const double *land_frac_in,
                 const double *q_ca_old_in,
                 double       *carbon_flux_out,
                 double       *q_ca_out,
                 double       *runoff_out,
                 double       *roughness_length_out,
                 double       *surface_temp_out,
                 double       *albedo_out,
                 double       *evap1_out,
                 double       *dedt_out,
                 double       *gc_out,
                 double       *exlpj_out);

void lpj_end_(void);

/*---------------------------------------------------------------------------*/
/* Module-level state: set by lpjml_init, used to validate array sizes       */
/*---------------------------------------------------------------------------*/

static int g_ncells = -1;

/*---------------------------------------------------------------------------*/
/* Helper: return a C-contiguous float64 1-D numpy array of length g_ncells. */
/* For output arrays pass writeable=1 so we write into the caller's buffer   */
/* rather than a temporary copy.                                              */
/* Returns a new reference on success, NULL (with exception) on failure.     */
/*---------------------------------------------------------------------------*/

static PyArrayObject *as_double_1d(PyObject *obj, const char *name,
                                    int writeable)
{
  int flags = NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED;
  if (writeable) {
    flags |= NPY_ARRAY_WRITEABLE;
  }
  PyArrayObject *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_DOUBLE,
                                                          flags);
  if (!arr) {
    if (!PyErr_Occurred())
      PyErr_Format(PyExc_TypeError,
                   "argument '%s' must be a contiguous float64 array", name);
    return NULL;
  }
  if (PyArray_NDIM(arr) != 1 || (int)PyArray_SIZE(arr) != g_ncells) {
    PyErr_Format(PyExc_ValueError,
                 "argument '%s' must be 1-D with %d elements "
                 "(got ndim=%d, size=%zd)",
                 name, g_ncells,
                 PyArray_NDIM(arr), (Py_ssize_t)PyArray_SIZE(arr));
    Py_DECREF(arr);
    return NULL;
  }
  return arr;
}

/*---------------------------------------------------------------------------*/
/* lpjml_init(config_filename, dt_fast=3600) -> dict                         */
/*                                                                           */
/*   Initialises LPJmL by calling lpj_init_() in driver.c.                   */
/*   dt_fast is the fast-timestep length in seconds; default 3600 (one hour) */
/*   so that every lpjml_update() call performs a full day of physics.       */
/*                                                                           */
/*   Returns {"ncells": int, "firstyear": int, "lastyear": int,              */
/*             "nspinup": int}                                                */
/*---------------------------------------------------------------------------*/

static PyObject *py_lpjml_init(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"config_filename", "dt_fast", NULL};
  const char *config_filename;
  int dt_fast = 3600;

  if (g_ncells >= 0) {
    PyErr_SetString(PyExc_RuntimeError,
                    "lpjml_init() already called; call lpjml_end() first");
    return NULL;
  }

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|i", kwlist,
                                    &config_filename, &dt_fast)) {
    return NULL;
  }

#ifdef USE_MPI
  {
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) {
      MPI_Init(NULL, NULL);
    }
  }
#endif

  int ncells = 0, firstyear = 0, lastyear = 0, nspinup = 0;
  lpj_init_(&dt_fast, config_filename, &ncells, &firstyear, &lastyear, &nspinup);
  g_ncells = ncells;

  int rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  return Py_BuildValue("{s:i,s:i,s:i,s:i,s:i}",
                       "ncells",    ncells,
                       "firstyear", firstyear,
                       "lastyear",  lastyear,
                       "nspinup",   nspinup,
                       "rank",      rank);
}

/*---------------------------------------------------------------------------*/
/* lpjml_update(year, month, day, hour, minute, seconds,                     */
/*              prec_in, temp_mean_in, swdown_in, lwnet_in,                  */
/*              t_flux_in, qflux_in, dedq_in, dhdt_in, drdt_in,              */
/*              drag_q_in, p_surf_in, wind_in, co2_in,                       */
/*              land_frac_in, q_ca_old_in,                                   */
/*              carbon_flux_out, q_ca_out, runoff_out,                       */
/*              roughness_length_out, surface_temp_out,                      */
/*              albedo_out, evap1_out, dedt_out, gc_out, exlpj_out)          */
/*                                                                           */
/*   Advances LPJmL by one fast timestep.  Calls lpj_update_() in driver.c. */
/*                                                                           */
/*   Time args are plain Python ints.                                        */
/*   All array args must be 1-D float64 numpy arrays of length ncells.      */
/*   Output arrays are modified in-place; pass pre-allocated arrays.        */
/*   co2_in may be None, in which case CO2 is read from file.               */
/*   Returns None.                                                           */
/*---------------------------------------------------------------------------*/

/* Descriptor for one array argument: source PyObject, destination
   PyArrayObject, argument name, and whether it must be writeable.    */
typedef struct {
  PyObject      **obj;
  PyArrayObject **arr;
  const char     *name;
  int             writeable;
} ArrSpec;

static PyObject *py_lpjml_update(PyObject *self, PyObject *args, PyObject *kwargs)
{
  if (g_ncells <= 0) {
    PyErr_SetString(PyExc_RuntimeError,
                    "lpjml_init() must be called before lpjml_update()");
    return NULL;
  }

  static char *kwlist[] = {
    "year", "month", "day", "hour", "minute", "seconds",
    "prec_in", "temp_mean_in", "swdown_in", "lwnet_in",
    "t_flux_in", "qflux_in", "dedq_in", "dhdt_in", "drdt_in",
    "drag_q_in", "p_surf_in", "wind_in", "co2_in",
    "land_frac_in", "q_ca_old_in",
    "carbon_flux_out", "q_ca_out", "runoff_out",
    "roughness_length_out", "surface_temp_out",
    "albedo_out", "evap1_out", "dedt_out", "gc_out", "exlpj_out",
    NULL
  };

  /* time arguments */
  int year, month, day, hour, minute, seconds;

  /* raw Python objects for array arguments */
  PyObject *prec_obj, *temp_obj, *swdown_obj, *lwnet_obj;
  PyObject *t_flux_obj, *qflux_obj, *dedq_obj, *dhdt_obj, *drdt_obj;
  PyObject *drag_q_obj, *p_surf_obj, *wind_obj, *co2_obj;
  PyObject *land_frac_obj, *q_ca_old_obj;
  PyObject *carbon_flux_obj, *q_ca_obj, *runoff_obj, *roughness_obj;
  PyObject *surface_temp_obj, *albedo_obj, *evap1_obj, *dedt_obj;
  PyObject *gc_obj, *exlpj_obj;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                    "iiiiii"           /* 6 time ints */
                                    "OOOOOOOOOOOOOOO"  /* 15 input arrays */
                                    "OOOOOOOOOO",      /* 10 output arrays */
                                    kwlist,
                                    &year, &month, &day, &hour, &minute, &seconds,
                                    &prec_obj, &temp_obj, &swdown_obj, &lwnet_obj,
                                    &t_flux_obj, &qflux_obj, &dedq_obj,
                                    &dhdt_obj, &drdt_obj, &drag_q_obj,
                                    &p_surf_obj, &wind_obj, &co2_obj,
                                    &land_frac_obj, &q_ca_old_obj,
                                    &carbon_flux_obj, &q_ca_obj, &runoff_obj,
                                    &roughness_obj, &surface_temp_obj, &albedo_obj,
                                    &evap1_obj, &dedt_obj, &gc_obj, &exlpj_obj)) {
    return NULL;
  }

  /* all array pointers initialised to NULL so the cleanup loop is always safe */
  PyArrayObject *prec_in_arr          = NULL, *temp_mean_in_arr    = NULL;
  PyArrayObject *swdown_in_arr        = NULL, *lwnet_in_arr         = NULL;
  PyArrayObject *t_flux_in_arr        = NULL, *qflux_in_arr         = NULL;
  PyArrayObject *dedq_in_arr          = NULL, *dhdt_in_arr          = NULL;
  PyArrayObject *drdt_in_arr          = NULL, *drag_q_in_arr        = NULL;
  PyArrayObject *p_surf_in_arr        = NULL, *wind_in_arr          = NULL;
  PyArrayObject *land_frac_in_arr     = NULL, *q_ca_old_in_arr      = NULL;
  PyArrayObject *co2_in_arr           = NULL;
  PyArrayObject *carbon_flux_out_arr  = NULL, *q_ca_out_arr         = NULL;
  PyArrayObject *runoff_out_arr       = NULL, *roughness_out_arr    = NULL;
  PyArrayObject *surface_temp_out_arr = NULL, *albedo_out_arr       = NULL;
  PyArrayObject *evap1_out_arr        = NULL, *dedt_out_arr         = NULL;
  PyArrayObject *gc_out_arr           = NULL, *exlpj_out_arr        = NULL;

  /* Descriptor table: one entry per array argument.
     The loop below iterates until it hits the sentinel (obj == NULL).  */
  ArrSpec specs[] = {
    {&prec_obj,        &prec_in_arr,           "prec_in",             0},
    {&temp_obj,        &temp_mean_in_arr,       "temp_mean_in",        0},
    {&swdown_obj,      &swdown_in_arr,          "swdown_in",           0},
    {&lwnet_obj,       &lwnet_in_arr,           "lwnet_in",            0},
    {&t_flux_obj,      &t_flux_in_arr,          "t_flux_in",           0},
    {&qflux_obj,       &qflux_in_arr,           "qflux_in",            0},
    {&dedq_obj,        &dedq_in_arr,            "dedq_in",             0},
    {&dhdt_obj,        &dhdt_in_arr,            "dhdt_in",             0},
    {&drdt_obj,        &drdt_in_arr,            "drdt_in",             0},
    {&drag_q_obj,      &drag_q_in_arr,          "drag_q_in",           0},
    {&p_surf_obj,      &p_surf_in_arr,          "p_surf_in",           0},
    {&wind_obj,        &wind_in_arr,            "wind_in",             0},
    {&co2_obj,         &co2_in_arr,             "co2_in",              0},
    {&land_frac_obj,   &land_frac_in_arr,       "land_frac_in",        0},
    {&q_ca_old_obj,    &q_ca_old_in_arr,        "q_ca_old_in",         0},
    {&carbon_flux_obj, &carbon_flux_out_arr,    "carbon_flux_out",     1},
    {&q_ca_obj,        &q_ca_out_arr,           "q_ca_out",            1},
    {&runoff_obj,      &runoff_out_arr,         "runoff_out",          1},
    {&roughness_obj,   &roughness_out_arr,      "roughness_length_out",1},
    {&surface_temp_obj,&surface_temp_out_arr,   "surface_temp_out",    1},
    {&albedo_obj,      &albedo_out_arr,         "albedo_out",          1},
    {&evap1_obj,       &evap1_out_arr,          "evap1_out",           1},
    {&dedt_obj,        &dedt_out_arr,           "dedt_out",            1},
    {&gc_obj,          &gc_out_arr,             "gc_out",              1},
    {&exlpj_obj,       &exlpj_out_arr,          "exlpj_out",           1},
    /* Sentinel entry: obj == NULL signals end-of-table to both loops below.
       The conversion loop stops here instead of reading past the array.
       The cleanup loop uses the same check so it always frees exactly the
       entries that were converted, no more and no less.
       IMPORTANT: this entry must always remain last in the table.          */
    {NULL, NULL, NULL, 0}
  };

  /* Convert and validate all regular arrays; stop at first failure. */
  int convert_ok = 1;
  for (int i = 0; specs[i].obj && convert_ok; i++) {
    *specs[i].arr = as_double_1d(*specs[i].obj, specs[i].name,
                                  specs[i].writeable);
    if (!*specs[i].arr) convert_ok = 0;
  }

  if (convert_ok) {
    lpj_update_(&year, &month, &day, &hour, &minute, &seconds,
                (const double *)PyArray_DATA(prec_in_arr),
                (const double *)PyArray_DATA(temp_mean_in_arr),
                (const double *)PyArray_DATA(swdown_in_arr),
                (const double *)PyArray_DATA(lwnet_in_arr),
                (const double *)PyArray_DATA(t_flux_in_arr),
                (const double *)PyArray_DATA(qflux_in_arr),
                (const double *)PyArray_DATA(dedq_in_arr),
                (const double *)PyArray_DATA(dhdt_in_arr),
                (const double *)PyArray_DATA(drdt_in_arr),
                (const double *)PyArray_DATA(drag_q_in_arr),
                (const double *)PyArray_DATA(p_surf_in_arr),
                (const double *)PyArray_DATA(wind_in_arr),
                (const double *)PyArray_DATA(co2_in_arr),
                (const double *)PyArray_DATA(land_frac_in_arr),
                (const double *)PyArray_DATA(q_ca_old_in_arr),
                (double *)PyArray_DATA(carbon_flux_out_arr),
                (double *)PyArray_DATA(q_ca_out_arr),
                (double *)PyArray_DATA(runoff_out_arr),
                (double *)PyArray_DATA(roughness_out_arr),
                (double *)PyArray_DATA(surface_temp_out_arr),
                (double *)PyArray_DATA(albedo_out_arr),
                (double *)PyArray_DATA(evap1_out_arr),
                (double *)PyArray_DATA(dedt_out_arr),
                (double *)PyArray_DATA(gc_out_arr),
                (double *)PyArray_DATA(exlpj_out_arr));
  }

  /* Unified cleanup: always runs; Py_XDECREF is a no-op on NULL. */
  for (int i = 0; specs[i].obj; i++)
    Py_XDECREF(*specs[i].arr);

  return convert_ok ? (Py_INCREF(Py_None), Py_None) : NULL;
}

/*---------------------------------------------------------------------------*/
/* lpjml_end() -> None                                                       */
/*                                                                           */
/*   Flushes and closes all output, frees LPJmL memory.                     */
/*   Calls lpj_end_() in driver.c.                                           */
/*---------------------------------------------------------------------------*/

static PyObject *py_lpjml_end(PyObject *self, PyObject *args)
{
  lpj_end_();
  g_ncells = -1;

#ifdef USE_MPI
  {
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) {
      MPI_Finalize();
    }
  }
#endif

  Py_RETURN_NONE;
}

/*---------------------------------------------------------------------------*/
/* Module method table and init                                              */
/*---------------------------------------------------------------------------*/

static PyMethodDef LPJmLCoupler_methods[] = {
  {"lpjml_init",
   (PyCFunction)py_lpjml_init,
   METH_VARARGS | METH_KEYWORDS,
   "lpjml_init(config_filename, dt_fast=3600) -> dict\n"
   "\n"
   "Initialise LPJmL.  Returns a dict with keys:\n"
   "  ncells    -- number of local land cells on this MPI rank\n"
   "  firstyear -- first simulation year\n"
   "  lastyear  -- last  simulation year\n"
   "  nspinup   -- number of spinup years\n"
   "  rank      -- MPI rank of this process (0 on non-MPI builds)\n"},

  {"lpjml_update",
   (PyCFunction)py_lpjml_update,
   METH_VARARGS | METH_KEYWORDS,
   "lpjml_update(year, month, day, hour, minute, seconds,\n"
   "             prec_in, temp_mean_in, swdown_in, lwnet_in,\n"
   "             t_flux_in, qflux_in, dedq_in, dhdt_in, drdt_in,\n"
   "             drag_q_in, p_surf_in, wind_in, co2_in,\n"
   "             land_frac_in, q_ca_old_in,\n"
   "             carbon_flux_out, q_ca_out, runoff_out,\n"
   "             roughness_length_out, surface_temp_out,\n"
   "             albedo_out, evap1_out, dedt_out, gc_out, exlpj_out)\n"
   "\n"
   "Advance LPJmL by one fast timestep.  All array arguments must be\n"
   "1-D contiguous float64 numpy arrays of length ncells.  Output arrays\n"
   "are modified in-place.\n"},

  {"lpjml_end",
   py_lpjml_end,
   METH_NOARGS,
   "lpjml_end() -> None\n"
   "\n"
   "Flush output, free LPJmL memory, and (if MPI was initialised here)\n"
   "call MPI_Finalize.\n"},

  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef LPJmLCoupler_module = {
  PyModuleDef_HEAD_INIT,
  "LPJmLCoupler",
  "Python bindings for the LPJmL standalone driver (driver.c).",
  -1,
  LPJmLCoupler_methods
};

PyMODINIT_FUNC PyInit_LPJmLCoupler(void)
{
  import_array(); /* initialise NumPy C-API; returns NULL on failure */
  return PyModule_Create(&LPJmLCoupler_module);
}

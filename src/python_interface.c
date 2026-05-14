/**************************************************************************************/
/**                                                                                \n**/
/**     p y t h o n _ i n t e r f a c e . c                                        \n**/
/**                                                                                \n**/
/**     NumPy C-API wrapper exposing driver.c as the Python module                 \n**/
/**     "LPJmLCoupler".                                                             \n**/
/**                                                                                \n**/
/**     Three functions:                                                            \n**/
/**       lpjml_init(config_filename, dt_fast=86400) -> dict                       \n**/
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
/* lpjml_init(config_filename, dt_fast=86400) -> dict                        */
/*                                                                           */
/*   Initialises LPJmL by calling lpj_init_() in driver.c.                  */
/*   dt_fast is the fast-timestep length in seconds; default 86400 (one day) */
/*   so that every lpjml_update() call performs a full day of physics.       */
/*                                                                           */
/*   Returns {"ncells": int, "firstyear": int, "lastyear": int,              */
/*             "nspinup": int}                                                */
/*---------------------------------------------------------------------------*/

static PyObject *py_lpjml_init(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"config_filename", "dt_fast", NULL};
  const char *config_filename;
  int dt_fast = 86400;

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

  return Py_BuildValue("{s:i,s:i,s:i,s:i}",
                       "ncells",    ncells,
                       "firstyear", firstyear,
                       "lastyear",  lastyear,
                       "nspinup",   nspinup);
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

static PyObject *py_lpjml_update(PyObject *self, PyObject *args)
{
  if (g_ncells < 0) {
    PyErr_SetString(PyExc_RuntimeError,
                    "lpjml_init() must be called before lpjml_update()");
    return NULL;
  }

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

  if (!PyArg_ParseTuple(args,
                         "iiiiii"          /* 6 time ints */
                         "OOOOOOOOOOOOOO"  /* 14 input arrays + co2 */
                         "OOOOOOOOOO",     /* 10 output arrays */
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

  /* all array pointers initialised to NULL for safe Py_XDECREF in cleanup */
  PyArrayObject *prec_in_arr      = NULL, *temp_mean_in_arr = NULL;
  PyArrayObject *swdown_in_arr    = NULL, *lwnet_in_arr     = NULL;
  PyArrayObject *t_flux_in_arr    = NULL, *qflux_in_arr     = NULL;
  PyArrayObject *dedq_in_arr      = NULL, *dhdt_in_arr      = NULL;
  PyArrayObject *drdt_in_arr      = NULL, *drag_q_in_arr    = NULL;
  PyArrayObject *p_surf_in_arr    = NULL, *wind_in_arr      = NULL;
  PyArrayObject *land_frac_in_arr = NULL, *q_ca_old_in_arr  = NULL;
  PyArrayObject *co2_in_arr       = NULL;

  PyArrayObject *carbon_flux_out_arr      = NULL, *q_ca_out_arr      = NULL;
  PyArrayObject *runoff_out_arr           = NULL, *roughness_out_arr = NULL;
  PyArrayObject *surface_temp_out_arr     = NULL, *albedo_out_arr    = NULL;
  PyArrayObject *evap1_out_arr            = NULL, *dedt_out_arr      = NULL;
  PyArrayObject *gc_out_arr               = NULL, *exlpj_out_arr     = NULL;

  /* convert and validate all arrays */
  prec_in_arr      = as_double_1d(prec_obj,      "prec_in",      0); if (!prec_in_arr)      goto cleanup;
  temp_mean_in_arr = as_double_1d(temp_obj,      "temp_mean_in", 0); if (!temp_mean_in_arr) goto cleanup;
  swdown_in_arr    = as_double_1d(swdown_obj,    "swdown_in",    0); if (!swdown_in_arr)    goto cleanup;
  lwnet_in_arr     = as_double_1d(lwnet_obj,     "lwnet_in",     0); if (!lwnet_in_arr)     goto cleanup;
  t_flux_in_arr    = as_double_1d(t_flux_obj,    "t_flux_in",    0); if (!t_flux_in_arr)    goto cleanup;
  qflux_in_arr     = as_double_1d(qflux_obj,     "qflux_in",     0); if (!qflux_in_arr)     goto cleanup;
  dedq_in_arr      = as_double_1d(dedq_obj,      "dedq_in",      0); if (!dedq_in_arr)      goto cleanup;
  dhdt_in_arr      = as_double_1d(dhdt_obj,      "dhdt_in",      0); if (!dhdt_in_arr)      goto cleanup;
  drdt_in_arr      = as_double_1d(drdt_obj,      "drdt_in",      0); if (!drdt_in_arr)      goto cleanup;
  drag_q_in_arr    = as_double_1d(drag_q_obj,    "drag_q_in",    0); if (!drag_q_in_arr)    goto cleanup;
  p_surf_in_arr    = as_double_1d(p_surf_obj,    "p_surf_in",    0); if (!p_surf_in_arr)    goto cleanup;
  wind_in_arr      = as_double_1d(wind_obj,      "wind_in",      0); if (!wind_in_arr)      goto cleanup;
  land_frac_in_arr = as_double_1d(land_frac_obj, "land_frac_in", 0); if (!land_frac_in_arr) goto cleanup;
  q_ca_old_in_arr  = as_double_1d(q_ca_old_obj,  "q_ca_old_in",  0); if (!q_ca_old_in_arr)  goto cleanup;

  carbon_flux_out_arr  = as_double_1d(carbon_flux_obj,  "carbon_flux_out",      1); if (!carbon_flux_out_arr)  goto cleanup;
  q_ca_out_arr         = as_double_1d(q_ca_obj,         "q_ca_out",             1); if (!q_ca_out_arr)         goto cleanup;
  runoff_out_arr       = as_double_1d(runoff_obj,       "runoff_out",           1); if (!runoff_out_arr)       goto cleanup;
  roughness_out_arr    = as_double_1d(roughness_obj,    "roughness_length_out", 1); if (!roughness_out_arr)    goto cleanup;
  surface_temp_out_arr = as_double_1d(surface_temp_obj, "surface_temp_out",     1); if (!surface_temp_out_arr) goto cleanup;
  albedo_out_arr       = as_double_1d(albedo_obj,       "albedo_out",           1); if (!albedo_out_arr)       goto cleanup;
  evap1_out_arr        = as_double_1d(evap1_obj,        "evap1_out",            1); if (!evap1_out_arr)        goto cleanup;
  dedt_out_arr         = as_double_1d(dedt_obj,         "dedt_out",             1); if (!dedt_out_arr)         goto cleanup;
  gc_out_arr           = as_double_1d(gc_obj,           "gc_out",               1); if (!gc_out_arr)           goto cleanup;
  exlpj_out_arr        = as_double_1d(exlpj_obj,        "exlpj_out",            1); if (!exlpj_out_arr)        goto cleanup;

  /* co2_in: None → NULL (driver reads CO2 from file via getco2()) */
  if (co2_obj != Py_None) {
    co2_in_arr = as_double_1d(co2_obj, "co2_in", 0);
    if (!co2_in_arr) {
      goto cleanup;
    }
  }

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
              co2_in_arr ? (const double *)PyArray_DATA(co2_in_arr) : NULL,
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

  Py_XDECREF(co2_in_arr);
  Py_DECREF(prec_in_arr);      Py_DECREF(temp_mean_in_arr);
  Py_DECREF(swdown_in_arr);    Py_DECREF(lwnet_in_arr);
  Py_DECREF(t_flux_in_arr);    Py_DECREF(qflux_in_arr);
  Py_DECREF(dedq_in_arr);      Py_DECREF(dhdt_in_arr);
  Py_DECREF(drdt_in_arr);      Py_DECREF(drag_q_in_arr);
  Py_DECREF(p_surf_in_arr);    Py_DECREF(wind_in_arr);
  Py_DECREF(land_frac_in_arr); Py_DECREF(q_ca_old_in_arr);
  Py_DECREF(carbon_flux_out_arr);  Py_DECREF(q_ca_out_arr);
  Py_DECREF(runoff_out_arr);       Py_DECREF(roughness_out_arr);
  Py_DECREF(surface_temp_out_arr); Py_DECREF(albedo_out_arr);
  Py_DECREF(evap1_out_arr);        Py_DECREF(dedt_out_arr);
  Py_DECREF(gc_out_arr);           Py_DECREF(exlpj_out_arr);
  Py_RETURN_NONE;

cleanup:
  Py_XDECREF(co2_in_arr);
  Py_XDECREF(prec_in_arr);      Py_XDECREF(temp_mean_in_arr);
  Py_XDECREF(swdown_in_arr);    Py_XDECREF(lwnet_in_arr);
  Py_XDECREF(t_flux_in_arr);    Py_XDECREF(qflux_in_arr);
  Py_XDECREF(dedq_in_arr);      Py_XDECREF(dhdt_in_arr);
  Py_XDECREF(drdt_in_arr);      Py_XDECREF(drag_q_in_arr);
  Py_XDECREF(p_surf_in_arr);    Py_XDECREF(wind_in_arr);
  Py_XDECREF(land_frac_in_arr); Py_XDECREF(q_ca_old_in_arr);
  Py_XDECREF(carbon_flux_out_arr);  Py_XDECREF(q_ca_out_arr);
  Py_XDECREF(runoff_out_arr);       Py_XDECREF(roughness_out_arr);
  Py_XDECREF(surface_temp_out_arr); Py_XDECREF(albedo_out_arr);
  Py_XDECREF(evap1_out_arr);        Py_XDECREF(dedt_out_arr);
  Py_XDECREF(gc_out_arr);           Py_XDECREF(exlpj_out_arr);
  return NULL;
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
   "lpjml_init(config_filename, dt_fast=86400) -> dict\n"
   "\n"
   "Initialise LPJmL.  Returns a dict with keys:\n"
   "  ncells    -- number of local land cells on this MPI rank\n"
   "  firstyear -- first simulation year\n"
   "  lastyear  -- last  simulation year\n"
   "  nspinup   -- number of spinup years\n"},

  {"lpjml_update",
   py_lpjml_update,
   METH_VARARGS,
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
   "are modified in-place.  Pass co2_in=None to read CO2 from file.\n"},

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

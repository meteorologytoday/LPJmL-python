/**************************************************************************************/
/**                                                                                \n**/
/**         l p j _ s t a n d a l o n e . c                                        \n**/
/**                                                                                \n**/
/**     Standalone entry point for LPJmL that uses driver.c.                      \n**/
/**     No FMS, no grid remapping, no domain decomposition knowledge required.     \n**/
/**     Climate inputs are read from files as configured in the LPJmL .cjson.     \n**/
/**                                                                                \n**/
/**     All arrays passed to lpj_update_() are in LPJmL native         \n**/
/**     1-D cell order (size = global_number_lpj_cells / ngridcell per PE).       \n**/
/**                                                                                \n**/
/**     FMS-specific coupling variables (t_flux, qflux, dedq, ...) are set to     \n**/
/**     physically neutral defaults; they are not read from files by LPJmL.       \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file   \n**/
/** authors, and contributors see AUTHORS file                                    \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3              \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/              \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                   \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include <string.h>
#include "lpj.h"

/* ---- forward declarations for driver.c ---- */
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

/* ---- default config filename (overridden by argv[1] at runtime) ---- */
#define LPJML_CONFIG_FILENAME "lpjml_fms.cjson"

/* ---- defaults for FMS atmospheric coupling variables ---- */
/* These have no file-based equivalent in LPJmL and must be prescribed. */
#define DEFAULT_T_FLUX   0.0       /* sensible heat flux         [W m-2]   */
#define DEFAULT_QFLUX    0.0       /* water-vapour flux          [kg m-2 s-1] */
#define DEFAULT_DEDQ     1e-2      /* dE/dq moisture sensitivity            */
#define DEFAULT_DHDT     0.0       /* dH/dt                      [W m-2]   */
#define DEFAULT_DRDT     0.0       /* dR/dt                      [W m-2]   */
#define DEFAULT_DRAG_Q   0.01      /* drag coefficient for moisture [m s-1] */
#define DEFAULT_P_SURF   101325.0  /* surface pressure           [Pa]      */

/* ---- helper: allocate and fill a double array ---- */
static double *allocate_array_double(double fill, int n)
{
  int i;
  double *ptr = malloc(n * sizeof(double));
  if (!ptr) { perror("allocate_array_double"); exit(EXIT_FAILURE); }
  for (i = 0; i < n; i++) {
    ptr[i] = fill;
  }
  return ptr;
}

/* ---- main ---- */
int main(int argc, char **argv)
{
  int ncells;
  int firstyear, lastyear, nspinup;
  int dt_fast = NSECONDSDAY;   /* one call per day: goodmorning and goodnight
                                   are always true on every lpj_update_() call */

  /* time variables passed to lpj_update_() */
  int Time_year, Time_month, Time_day;
  int Time_hour = 0, Time_minute = 0, Time_second = 0;

  /* -----------------------------------------------------------------------
   * MPI initialisation
   * ----------------------------------------------------------------------- */
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  /* -----------------------------------------------------------------------
   * Initialise LPJmL
   * argv[1], if provided, overrides the compiled-in default config filename.
   * ----------------------------------------------------------------------- */
  const char *config_filename = (argc > 1) ? argv[1] : LPJML_CONFIG_FILENAME;
  lpj_init_(&dt_fast, config_filename,
                       &ncells, &firstyear, &lastyear, &nspinup);

  /* -----------------------------------------------------------------------
   * Allocate LPJ-native 1-D arrays (one value per local land cell)
   *
   * Climate inputs (prec, temp_mean, swdown, lwnet, wind) are passed as 0
   * because lpj_driver calls dailyclimate() internally, which
   * overwrites these with values read from the files configured in the
   * .cjson.  The zeros are never used for actual physics.
   *
   * FMS coupling variables have no file equivalent and are set to the
   * physically neutral defaults defined above.
   * ----------------------------------------------------------------------- */

  /* inputs */
  double *prec_in      = allocate_array_double(0.0,            ncells);
  double *temp_mean_in = allocate_array_double(0.0,            ncells);
  double *swdown_in    = allocate_array_double(0.0,            ncells);
  double *lwnet_in     = allocate_array_double(0.0,            ncells);
  double *wind_in      = allocate_array_double(0.0,            ncells);
  double *t_flux_in    = allocate_array_double(DEFAULT_T_FLUX,  ncells);
  double *qflux_in     = allocate_array_double(DEFAULT_QFLUX,   ncells);
  double *dedq_in      = allocate_array_double(DEFAULT_DEDQ,    ncells);
  double *dhdt_in      = allocate_array_double(DEFAULT_DHDT,    ncells);
  double *drdt_in      = allocate_array_double(DEFAULT_DRDT,    ncells);
  double *drag_q_in    = allocate_array_double(DEFAULT_DRAG_Q,  ncells);
  double *p_surf_in    = allocate_array_double(DEFAULT_P_SURF,  ncells);
  double *land_frac_in = allocate_array_double(1.0,             ncells); /* full land fraction */
  double *q_ca_old_in  = allocate_array_double(0.0,             ncells); /* updated each step */

  /* outputs */
  double *carbon_flux_out      = allocate_array_double(0.0, ncells);
  double *q_ca_out             = allocate_array_double(0.0, ncells);
  double *runoff_out           = allocate_array_double(0.0, ncells);
  double *roughness_length_out = allocate_array_double(0.0, ncells);
  double *surface_temp_out     = allocate_array_double(0.0, ncells);
  double *albedo_out           = allocate_array_double(0.0, ncells);
  double *evap1_out            = allocate_array_double(0.0, ncells);
  double *dedt_out             = allocate_array_double(0.0, ncells);
  double *gc_out               = allocate_array_double(0.0, ncells);
  double *exlpj_out            = allocate_array_double(0.0, ncells);

  /* -----------------------------------------------------------------------
   * Time loop: year → month → day
   * lpj_update_() derives goodmorning / goodnight / happynewyear
   * / silvester etc. from the time arguments internally.  With dt_fast =
   * NSECONDSDAY and Time_hour = Time_minute = Time_second = 0, both
   * goodmorning and goodnight are always true, so every call does a
   * complete day's worth of physics.
   * ----------------------------------------------------------------------- */
  for (Time_year  = firstyear - nspinup; Time_year  <= lastyear;            Time_year++) {
    for (Time_month = 1;                   Time_month <= NMONTH;              Time_month++) {
      for (Time_day   = 1;                   Time_day   <= ndaymonth[Time_month-1]; Time_day++) {

        // Hi Claude, Add another loop for hourly. Mimicking lpjfmstest.c
        
        /* carry q_ca forward as the previous timestep's output */
        memcpy(q_ca_old_in, q_ca_out, ncells * sizeof(double));

        lpj_update_(
                  &Time_year, &Time_month, &Time_day,
                  &Time_hour, &Time_minute, &Time_second,
                  prec_in, temp_mean_in, swdown_in, lwnet_in,
                  t_flux_in, qflux_in, dedq_in, dhdt_in, drdt_in,
                  drag_q_in, p_surf_in, wind_in,
                  NULL,          /* co2_in: NULL -> CO2 read from file via getco2() */
                  land_frac_in, q_ca_old_in,
                  carbon_flux_out, q_ca_out, runoff_out,
                  roughness_length_out, surface_temp_out,
                  albedo_out, evap1_out, dedt_out, gc_out, exlpj_out);
      }
    }
  }

  /* -----------------------------------------------------------------------
   * Finalise
   * ----------------------------------------------------------------------- */
  lpj_end_();

  free(prec_in);      free(temp_mean_in); free(swdown_in);
  free(lwnet_in);     free(wind_in);      free(t_flux_in);
  free(qflux_in);     free(dedq_in);      free(dhdt_in);
  free(drdt_in);      free(drag_q_in);    free(p_surf_in);
  free(land_frac_in); free(q_ca_old_in);
  free(carbon_flux_out);      free(q_ca_out);
  free(runoff_out);           free(roughness_length_out);
  free(surface_temp_out);     free(albedo_out);
  free(evap1_out);            free(dedt_out);
  free(gc_out);               free(exlpj_out);

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
} /* of 'main' */

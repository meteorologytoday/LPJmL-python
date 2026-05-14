/**************************************************************************************/
/**                                                                                \n**/
/**     d r i v e r . c                                                            \n**/
/**                                                                                \n**/
/**     Pure LPJmL model logic extracted from lpj_climber4.c.                     \n**/
/**     All inputs and outputs are in LPJmL native 1-D cell order                 \n**/
/**     (size config.ngridcell).  No FMS grid geometry, no CPL remapping.         \n**/
/**                                                                                \n**/
/**     Three entry points:                                                        \n**/
/**       lpj_init_()   -- read config, allocate grid, open I/O        \n**/
/**       lpj_update_() -- advance one fast timestep                   \n**/
/**       lpj_end_()    -- flush output, free memory                   \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file   \n**/
/** authors, and contributors see AUTHORS file                                    \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3              \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/              \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                   \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include <time.h>
#include <fenv.h>
#include <math.h>
#include <string.h>
#include "lpj.h"
#include "grass.h"
#include "tree.h"
#include "crop.h"
#include "natural.h"
#include "grassland.h"
#include "biomass_tree.h"
#include "biomass_grass.h"
#include "agriculture.h"
#include "agriculture_tree.h"
#include "agriculture_grass.h"
#include "woodplantation.h"

#define LPJML_CONFIG_FILENAME "lpjml_fms.cjson"
#define NTYPES       3   /* grass, tree, crop */
#define NSTANDTYPES 13   /* land use types as defined in landuse.h */

/*---------------------------------------------------------------------------*/
/* LPJmL model state (static globals, same as lpj_climber4.c)               */
/*---------------------------------------------------------------------------*/

static const char  *progname;
static Outputfile  *output;
static Input        input;
static Config       config;
static Cell        *grid;
static Standtype    standtype[NSTANDTYPES];
static int          npft, ncft;
static Real         co2, cflux_total;
static Flux         flux;
static int          year;
static int          dt_fast;       /* fast timestep [s] passed by caller */
static Real         dt_fast_real;

#ifdef USE_MPI
static MPI_Comm comm;
#endif

/*---------------------------------------------------------------------------*/
/* Internal helpers                                                          */
/*---------------------------------------------------------------------------*/

static void *getmem(size_t size)
{
  void *ptr = malloc(size);
  if (!ptr) { perror("lpj_driver: getmem"); exit(EXIT_FAILURE); }
  return ptr;
}

/*---------------------------------------------------------------------------*/
/* lpj_init_                                                      */
/*                                                                           */
/* Initialise LPJmL: read config, allocate grid, open input/output files.   */
/* Must be called once before the time loop.                                 */
/*                                                                           */
/* Parameters                                                                */
/*   dt_fast_in              [in]  fast timestep length in seconds           */
/*   global_number_lpj_cells [out] total number of LPJ land cells (nall)    */
/*---------------------------------------------------------------------------*/
#ifdef __cplusplus
extern "C"
#endif
void lpj_init_(const int  *dt_fast_in,
                           const char *config_filename,
                           int        *global_number_lpj_cells,
                           int        *firstyear_out,
                           int        *lastyear_out,
                           int        *nspinup_out)
{
  int    argc = 2;
  char **argv;
  int    rc, cell;

  Pfttype scanfcn[NTYPES] =
  {
    {name_grass, fscanpft_grass},
    {name_tree,  fscanpft_tree},
    {name_crop,  fscanpft_crop}
  };

  dt_fast      = *dt_fast_in;
  dt_fast_real = (Real)dt_fast;

  /* build a minimal argv so readconfig() can find the config file */
  argv = malloc((argc + 1) * sizeof(*argv));
  if (!argv) { perror("cannot malloc argv"); exit(EXIT_FAILURE); }
  argv[argc] = NULL;
  argv[0]    = strdup("lpjml");
  argv[1]    = strdup(config_filename);

  enablefpe();

#ifdef USE_MPI
  comm = MPI_COMM_WORLD;
  initmpiconfig(&config, comm);
#else
  initconfig(&config);
#endif

  progname = strippath(argv[0]);
  if (isroot(config)) {
    time_t tbegin;
    time(&tbegin);
    copyright(progname);
    printf("\nRunning for user %s on %s at %s",
           (getuser() == NULL) ? "N/A" : getuser(),
           gethost(), ctime(&tbegin));
    fflush(stdout);
  }

  rc = readconfig(&config, scanfcn, NTYPES, NSTANDTYPES, NOUT,
                  &argc, &argv, lpj_usage);
  failonerror(&config, rc, READ_CONFIG_ERR, "Cannot read configuration");

  if (isroot(config) && argc) {
    fprintf(stderr,
            "WARNING018: Arguments listed after configuration filename "
            "will be ignored.\n");
  }

  standtype[NATURAL]           = natural_stand;
  standtype[SETASIDE_RF]       = setaside_rf_stand;
  standtype[SETASIDE_IR]       = setaside_ir_stand;
  standtype[AGRICULTURE]       = agriculture_stand;
  standtype[MANAGEDFOREST]     = managedforest_stand;
  standtype[GRASSLAND]         = grassland_stand;
  standtype[OTHERS]            = others_stand;
  standtype[BIOMASS_TREE]      = biomass_tree_stand;
  standtype[BIOMASS_GRASS]     = biomass_grass_stand;
  standtype[AGRICULTURE_TREE]  = agriculture_tree_stand;
  standtype[AGRICULTURE_GRASS] = agriculture_grass_stand;
  standtype[WOODPLANTATION]    = woodplantation_stand;
  standtype[KILL]              = kill_stand;

  if (isroot(config)) {
    printconfig(config.npft[GRASS] + config.npft[TREE],
                config.npft[CROP], &config);
  }

  rc = ((grid = newgrid(&config, standtype, NSTANDTYPES,
                        config.npft[GRASS] + config.npft[TREE],
                        config.npft[CROP])) == NULL);
  failonerror(&config, rc, INIT_GRID_ERR,
              "Initialization of LPJ grid failed");

  rc = initinput(&input, grid,
                 config.npft[GRASS] + config.npft[TREE], &config);
  failonerror(&config, rc, INIT_INPUT_ERR,
              "Initialization of input data failed");

  output = fopenoutput(grid, NOUT, &config);
  rc = (output == NULL);
  failonerror(&config, rc, INIT_OUTPUT_ERR,
              "Initialization of output data failed");

  rc = initoutput(output, grid,
                  config.npft[GRASS] + config.npft[TREE],
                  config.npft[CROP], &config);
  failonerror(&config, rc, INIT_OUTPUT_ERR,
              "Initialization of output data failed");

  if (isopen(output, GRID)) {
    writecoords(output, GRID, grid, &config);
  }
  if (isopen(output, TERR_AREA)) {
    writearea(output, TERR_AREA, grid, &config);
  }
  if (isopen(output, LAKE_AREA)) {
    writearea(output, LAKE_AREA, grid, &config);
  }
  if (isopen(output, COUNTRY) && config.withlanduse) {
    writecountrycode(output, COUNTRY, grid, &config);
  }

  if (config.with_anomaly) {
    if (getclimate(input.climate, grid, config.firstyear, &config)) {
      fprintf(stderr, "ERROR104: Simulation stopped in getclimate().\n");
      fflush(stderr);
      abort();
    }
    for (cell = 0; cell < config.ngridcell; cell++) {
      grid[cell].aprec_poem_last = grid[cell].aprec_ref =
        getaprec(input.climate, cell);
    }
  }

  npft = config.npft[GRASS] + config.npft[TREE];
  ncft = config.npft[CROP];

  if (isroot(config)) {
    printf("Simulation begins...\n");
  }

  *global_number_lpj_cells = config.nall;
  *firstyear_out           = config.firstyear;
  *lastyear_out            = config.lastyear;
  *nspinup_out             = config.nspinup;
} /* of 'lpj_init_' */

/*---------------------------------------------------------------------------*/
/* lpj_update_                                                    */
/*                                                                           */
/* Advance LPJmL by one fast timestep.  All arrays are in LPJmL native      */
/* 1-D order, size config.ngridcell.  No grid remapping is performed here.  */
/*                                                                           */
/* co2_in may be NULL when CO2 is read from file (config.co2_from_file is   */
/* false, or input.climate->co2.data is non-NULL).                           */
/*---------------------------------------------------------------------------*/
#ifdef __cplusplus
extern "C"
#endif
void lpj_update_
(
  const int    const *Time_year,
  const int    const *Time_month,
  const int    const *Time_day,
  const int    const *Time_hour,
  const int    const *Time_minute,
  const int    const *Time_seconds,
  /* inputs -- LPJmL native 1-D, size config.ngridcell */
  const double const *prec_in,
  const double const *temp_mean_in,
  const double const *swdown_in,
  const double const *lwnet_in,
  const double const *t_flux_in,
  const double const *qflux_in,
  const double const *dedq_in,
  const double const *dhdt_in,
  const double const *drdt_in,
  const double const *drag_q_in,
  const double const *p_surf_in,
  const double const *wind_in,
  const double const *co2_in,    /* may be NULL when CO2 comes from file */
  const double const *land_frac_in,
  const double const *q_ca_old_in,
  /* outputs -- LPJmL native 1-D, size config.ngridcell */
  double *carbon_flux_out,
  double *q_ca_out,
  double *runoff_out,
  double *roughness_length_out,
  double *surface_temp_out,
  double *albedo_out,
  double *evap1_out,
  double *dedt_out,
  double *gc_out,
  double *exlpj_out
) {
  
  int i, s, cell;

  /* ---- time flags (mirror logic from lpj_climber4.c) ---- */
  {
    year = *Time_year;   /* must be first: used throughout the block */
    const int month      = *Time_month - 1;   /* dailyclimate.c: month (0..11) */
    const int dayofmonth = *Time_day   - 1;   /* 0..N  */
    static int yesterday = -1;

    const Bool goodmorning   = yesterday != dayofmonth;
    const Bool goodnight     = (*Time_hour * 3600 + *Time_minute * 60
                                + *Time_seconds + dt_fast) >= NSECONDSDAY;
    const Bool newmoon       = goodmorning && 0 == dayofmonth;
    const Bool monthend      = goodnight
                               && dayofmonth == ndaymonth[month] - 1;
    const Bool happynewyear  = newmoon && 0 == month;
    const Bool silvester     = (goodnight && yesterday != -1
                                && 11 == month
                                && 30 == dayofmonth);
    Bool use_spinup = FALSE;

    int climate_year, depos_year, landuse_year, wateruse_year;
    Stand *stand;

    feenableexcept(FE_INVALID);

    /* ---- start-of-year tasks ---- */
    if (!config.co2_from_file || input.climate->co2.data != NULL) {
      int year_co2;
      if (config.fix_co2 && year > config.fix_co2_year) {
        year_co2 = config.fix_co2_year;
      } else {
        year_co2 = year;
      }
      if (getco2(input.climate, &co2, year_co2, &config)) {
        fprintf(stderr, "ERROR015: Invalid year %d in getco2().\n", year_co2);
      }
    }

    if (happynewyear)
    {
      printf("LPJmL_CO2: %g\n", co2);

      if (config.fix_climate && year > config.fix_climate_year) {
        if (config.fix_climate_shuffle) {
          if (isroot(config)) {
            climate_year = config.fix_climate_interval[0]
              + (int)((config.fix_climate_interval[1]
                       - config.fix_climate_interval[0] + 1)
                      * erand48(config.seed));
          }
#ifdef USE_MPI
          MPI_Bcast(&climate_year, 1, MPI_INT, 0, config.comm);
#endif
        } else {
          climate_year = config.fix_climate_interval[0]
            + (year - config.fix_climate_year)
            % (config.fix_climate_interval[1]
               - config.fix_climate_interval[0] + 1);
        }
      } else {
        climate_year = year;
      }

      if (getclimate(input.climate, grid, climate_year, &config))
      {
        fprintf(stderr, "ERROR104: Simulation stopped in getclimate().\n");
        fflush(stderr);
        abort();
      }

      if (config.fix_deposition) {
        if (config.fix_deposition_with_climate) {
          depos_year = climate_year;
        } else if (year > config.fix_deposition_year) {
          if (config.fix_deposition_shuffle) {
            if (isroot(config)) {
              depos_year = config.fix_deposition_interval[0]
                + (int)((config.fix_deposition_interval[1]
                         - config.fix_deposition_interval[0] + 1)
                        * erand48(config.seed));
            }
#ifdef USE_MPI
            MPI_Bcast(&depos_year, 1, MPI_INT, 0, config.comm);
#endif
          } else {
            depos_year = config.fix_deposition_interval[0]
              + (year - config.fix_deposition_year)
              % (config.fix_deposition_interval[1]
                 - config.fix_deposition_interval[0] + 1);
          }
        } else {
          depos_year = year;
        }
      } else {
        depos_year = year;
      }

      if (getdeposition(input.climate, grid, depos_year, &config))
      {
        if (isroot(config))
        {
          fprintf(stderr,
                  "ERROR104: Simulation stopped in getdeposition().\n");
          fflush(stderr);
        }
        abort();
      }

      if (input.landuse != NULL) {
        calc_seasonality(grid, npft, ncft, &config);
        if (config.withlanduse == CONST_LANDUSE
            || config.withlanduse == ALL_CROPS) {
          landuse_year = config.landuse_year_const;
        } else {
          landuse_year = year;
        }
        if (config.withlanduse == CONST_LANDUSE) {
          wateruse_year = config.landuse_year_const;
        } else {
          wateruse_year = year;
        }
        if (getlanduse(input.landuse, grid, landuse_year, year,
                       ncft, &config)) {
          fprintf(stderr,
                  "ERROR104: Simulation stopped in getlanduse().\n");
          fflush(stderr);
          abort();
        }
        if (config.reservoir) {
          allocate_reservoir(grid, year, &config);
        }
      }

      if (input.wateruse != NULL && input.landuse != NULL) {
        if (getwateruse(input.wateruse, grid, wateruse_year, &config)) {
          fprintf(stderr,
                  "ERROR104: Simulation stopped in getwateruse().\n");
          fflush(stderr);
          abort();
        }
      }

      if (config.ispopulation) {
        if (readpopdens(input.popdens, year, grid, &config)) {
          fprintf(stderr,
                  "ERROR104: Simulation stopped in getpopdens().\n");
          fflush(stderr);
          abort();
        }
      }

      if (config.fire == SPITFIRE || config.fire == SPITFIRE_TMAX) {
        if (gethumanignition(input.human_ignition, year, grid, &config)) {
          if (isroot(config)) {
            fprintf(stderr,
                    "ERROR104: Simulation stopped in "
                    "gethumanignition().\n");
            fflush(stderr);
          }
          abort();
        }
      }

      if (config.prescribe_landcover != NO_LANDCOVER) {
        if (readlandcover(input.landcover, grid, year, &config)) {
          fprintf(stderr,
                  "ERROR104: Simulation stopped in readlandcover().\n");
          fflush(stderr);
          abort();
        }
      }

      if (year >= config.outputyear) {
        openoutput_yearly(output, year, &config);
      }

      if (year == input.climate->firstyear) {
        use_spinup = TRUE;
        printf("use temperature values from restart file\n");
      }
    } /* if (happynewyear) */

    /* ---- begin iterateyear_river() ---- */
    {
      static Dailyclimate fast;
      static Bool         intercrop;
      static int          dayofyear;
      Real                norg_soil_agr, nmin_soil_agr, nveg_soil_agr;
      static Real         popdens = 0;

      if (happynewyear)
      {
        popdens = 0;
        intercrop = getintercrop(input.landuse);
        for (cell = 0; cell < config.ngridcell; cell++)
        {
          grid[cell].daily.tmax = -1e38;
          grid[cell].daily.tmin =  1e38;
          grid[cell].daily.prec = grid[cell].daily.swdown =
            grid[cell].daily.lwnet = grid[cell].daily.windspeed =
            grid[cell].daily.temp = grid[cell].VPD =
            grid[cell].VPD2 = 0;
          initoutputdata(&grid[cell].output, ANNUAL, year, &config);
          grid[cell].balance.asoiltemp = 0;
          grid[cell].balance.surface_storage =
            grid[cell].balance.adischarge = 0;
          grid[cell].discharge.afin_ext = 0;
          if (!grid[cell].skip)
          {
            init_annual(grid + cell, ncft, &config);
            if (input.landuse != NULL)
            {
              if (grid[cell].lakefrac < 1)
              {
                if (!config.natural_only
                    && (year > config.firstyear - config.nspinup
                        || config.from_restart)) {
                  landusechange(grid + cell, npft, ncft,
                                intercrop, year, &config);
                } else if (grid[cell].ml.dam) {
                  landusechange_for_reservoir(grid + cell, npft, ncft,
                                              intercrop, year, &config);
                }
              }
            }
            initgdd(grid[cell].gdd, npft);
          }
        }
        dayofyear = 1;
      } /* if (happynewyear) */

      if (newmoon)
      {
        for (cell = 0; cell < config.ngridcell; cell++)
        {
          grid[cell].discharge.mfin = grid[cell].discharge.mfout =
            grid[cell].ml.mdemand = 0.0;
          grid[cell].output.mpet = 0;
          if (grid[cell].ml.dam) {
            grid[cell].ml.resdata->mprec_res = 0;
          }
          initoutputdata(&((grid + cell)->output), MONTHLY, year, &config);
          if (!grid[cell].skip) {
            initclimate_monthly(input.climate, &grid[cell].climbuf,
                                cell, month, grid[cell].seed);
          }
        }
      } /* if (newmoon) */

      if (goodnight) {
        printf("Date:  %04d-%02d-%02d   ",
               *Time_year, *Time_month, *Time_day);
      }

      if (goodmorning)
      {
        for (cell = 0; cell < config.ngridcell; cell++)
        {
          if (!grid[cell].skip)
          {
            grid[cell].discharge.drunoff = 0;
            grid[cell].output.evapotrans_flux_coupler = 0;
            grid[cell].output.dcflux = grid[cell].output.dwflux = 0;
            foreachstand(stand, s, grid[cell].standlist) {
              stand->interception = 0;
              stand->wet = newvec(Real, getnpft(&stand->pftlist));
              for (i = 0; i < getnpft(&stand->pftlist); i++) {
                stand->wet[i] = 0;
              }
            }
            initoutputdata(&(grid[cell].output), DAILY, year, &config);
          }
        }
      } /* if (goodmorning) */

      /* ---- fast-timestep cell loop: read inputs, call update_fast ---- */
      for (cell = 0; cell < config.ngridcell; cell++)
      {
        grid[cell].lad_landfrac = land_frac_in[cell];
        if (!grid[cell].skip)
        {
          fast.temp      = temp_mean_in[cell];
          grid[cell].daily.temp     += fast.temp;
          grid[cell].fasttemp        = temp_mean_in[cell];
          fast.prec      = prec_in[cell];
          grid[cell].fastprec        = prec_in[cell];
          fast.windspeed = wind_in[cell];
          grid[cell].daily.windspeed += fast.windspeed;
          fast.swdown    = swdown_in[cell];
          grid[cell].daily.swdown   += fast.swdown;
          fast.lwnet     = lwnet_in[cell];
          grid[cell].daily.lwnet    += fast.lwnet;
          fast.t_flux    = t_flux_in[cell];
          fast.qflux     = qflux_in[cell];
          fast.q_ca_old  = q_ca_old_in[cell];
          grid[cell].fastqflux       = qflux_in[cell];
          fast.dedq      = dedq_in[cell];
          grid[cell].fastdedq        = dedq_in[cell];
          if (grid[cell].daily.tmax < fast.temp) {
            grid[cell].daily.tmax = fast.temp;
          }
          if (grid[cell].daily.tmin > fast.temp) {
            grid[cell].daily.tmin = fast.temp;
          }
          fast.dhdt      = dhdt_in[cell];
          fast.drdt      = drdt_in[cell];
          fast.drag_q    = drag_q_in[cell];
          fast.p_surf    = p_surf_in[cell];

          if (config.co2_from_file && NULL == input.climate->co2.data) {
            co2 = co2_in[cell];
          }

          grid[cell].output.evapotrans_flux_fast = 0;
          update_fast(grid + cell, fast, dt_fast_real, dayofyear,
                      month, &config, use_spinup);
          grid[cell].daily.prec += fast.prec;

          getoutput(&grid[cell].output, TFLUX, &config) +=
            fast.t_flux * dt_fast_real / NSECONDSDAY;
          getoutput(&grid[cell].output, SWDOWN, &config) +=
            fast.swdown * dt_fast_real / NSECONDSDAY;
          getoutput(&grid[cell].output, LWNET, &config) +=
            fast.lwnet * dt_fast_real / NSECONDSDAY;
        }
      } /* fast-timestep cell loop */
      use_spinup = FALSE;

      /* ---- collect fast-timestep outputs ---- */
      for (cell = 0; cell < config.ngridcell; cell++)
      {
        if (!grid[cell].skip)
        {
          Real temp = 0.0, fracsum = 0.0;

          q_ca_out[cell] = grid[cell].q_ca;
          if (isnan(q_ca_out[cell]))
          {
            printf("ERR_CLIMBER: QCA NAN\n");
            q_ca_out[cell] = 0;
          }

          if (roughnesslength(grid[cell].standlist) < 0.001) {
            printf("ERR_CLIMBER: roughness smaller than 0.001\n");
            roughness_length_out[cell] = 0.001;
          } else if (roughnesslength(grid[cell].standlist) > 3) {
            printf("ERR_CLIMBER: roughness larger than 3\n");
            roughness_length_out[cell] = 3;
          } else {
            roughness_length_out[cell] =
              roughnesslength(grid[cell].standlist);
          }
          grid[cell].roughness_length = roughness_length_out[cell];
          if (isnan(roughness_length_out[cell])) {
            printf("ERR_CLIMBER: roughness NaN\n");
            roughness_length_out[cell] = 0.01;
          }
          if (roughness_length_out[cell] < 0.001) {
            printf("ERR_CLIMBER: roughness smaller than 0.001\n");
            roughness_length_out[cell] = 0.001;
          }

          foreachstand(stand, s, grid[cell].standlist) {
            temp     += stand->soil.surface_temp * stand->frac;
            fracsum  += stand->frac;
          }

          if (grid[cell].laketemp < -150)
          {
            printf("ERR_CLIMBER: laketemp too small %g\n",
                   grid[cell].laketemp);
            grid[cell].laketemp = -150;
          }
          if (grid[cell].laketemp > 150)
          {
            printf("ERR_CLIMBER: laketemp too large %g\n",
                   grid[cell].laketemp);
            grid[cell].laketemp = 150;
          }

          if (grid[cell].lakefrac + grid[cell].ml.reservoirfrac < 1) {
            surface_temp_out[cell] = temp / fracsum;
          } else {
            surface_temp_out[cell] = grid[cell].laketemp;
          }

          getoutput(&grid[cell].output, BETA, &config) +=
            grid[cell].albedo * dt_fast_real / NSECONDSDAY;
          getoutput(&grid[cell].output, SOIL_SURFACE_TEMP, &config) +=
            surface_temp_out[cell] * dt_fast_real / NSECONDSDAY;

          if (surface_temp_out[cell] < -150) {
            printf("ERR_CLIMBER: final soiltemp too small %g\n",
                   surface_temp_out[cell]);
            surface_temp_out[cell] = -150;
          }
          if (surface_temp_out[cell] > 80) {
            printf("ERR_CLIMBER: final soiltemp too large %g\n",
                   surface_temp_out[cell]);
            surface_temp_out[cell] = 80;
          }
          if (isnan(surface_temp_out[cell])) {
            printf("ERR_CLIMBER: TEMP NAN\n");
            surface_temp_out[cell] = 0;
          }
          grid[cell].surface_temp = surface_temp_out[cell];

          if (grid[cell].albedo < 0.05) {
            printf("ERR_CLIMBER: albedo smaller than 0.05\n");
            albedo_out[cell] = 0.05;
          } else if (grid[cell].albedo > 0.9) {
            printf("ERR_CLIMBER: albedo larger than 0.9\n");
            albedo_out[cell] = 0.9;
          } else if (isnan(grid[cell].albedo)) {
            printf("ERR_CLIMBER: albedo NAN\n");
            albedo_out[cell] = 0.2;
          } else {
            albedo_out[cell] = grid[cell].albedo;
          }
          if (isnan(albedo_out[cell])) {
            printf("ERR_CLIMBER: ALBEDO NAN\n");
            albedo_out[cell] = 0.05;
          }
          grid[cell].albedo = albedo_out[cell];

          evap1_out[cell]  = grid[cell].evap1;
          dedt_out[cell]   = grid[cell].dedt;
          gc_out[cell]     = grid[cell].gc;
          exlpj_out[cell]  = grid[cell].exlpj;
        }
      } /* output collection loop */

      /* ---- daily (goodnight) tasks ---- */
      if (goodnight) {
        Real temp_mean_g = 0, prec_mean_g = 0, qca_mean_g = 0,
             albedo_mean_g = 0, global_area = 0;
        Real water_total = 0, water_total_global = 0;

        for (cell = 0; cell < config.ngridcell; cell++) {
          if (!grid[cell].skip) {
            global_area += grid[cell].coord.area;
            getoutput(&grid[cell].output, PREC, &config) +=
              grid[cell].daily.prec;
            grid[cell].aprec_poem += grid[cell].daily.prec;

            grid[cell].daily.temp     *= dt_fast_real / NSECONDSDAY;
            grid[cell].daily.swdown   *= dt_fast_real / NSECONDSDAY;
            grid[cell].daily.lwnet    *= dt_fast_real / NSECONDSDAY;
            grid[cell].daily.windspeed *= dt_fast_real / NSECONDSDAY;

            dailyclimate(&grid[cell].daily, input.climate,
                         &grid[cell].climbuf, grid + cell, cell,
                         dayofyear, month, dayofmonth,
                         config.with_anomaly);

            if (config.ispopulation) {
              popdens = getpopdens(input.popdens, cell);
            }

            update_daily(grid + cell, co2, popdens, grid[cell].daily,
                         dayofyear, npft, ncft, year, month,
                         intercrop, &config);

            getoutput(&grid[cell].output, TEMP, &config) +=
              grid[cell].daily.temp;
            temp_mean_g    += grid[cell].daily.temp  * grid[cell].coord.area;
            prec_mean_g    += grid[cell].daily.prec  * grid[cell].coord.area;
            albedo_mean_g  += grid[cell].albedo      * grid[cell].coord.area;
            qca_mean_g     += grid[cell].q_ca        * grid[cell].coord.area;
          }
        }

        printf(" temp: %-.2f  prec: %-.3f  albedo: %-.3f  q_ca: %-.3f\n",
               temp_mean_g / global_area, prec_mean_g / global_area,
               albedo_mean_g / global_area, qca_mean_g / global_area);

        if (config.river_routing) {
          if (input.landuse != NULL || input.wateruse != NULL) {
            withdrawal_demand(grid, &config);
          }
          drain(grid, month, &config);
          if (input.landuse != NULL || input.wateruse != NULL) {
            wateruse(grid, npft, ncft, month, &config);
          }
        }

        if (config.withdailyoutput && year >= config.outputyear) {
          fwriteoutput(output, grid, year, dayofyear - 1, DAILY,
                       npft, ncft, &config);
        }

        /* collect runoff and water balance */
        for (cell = 0; cell < config.ngridcell; cell++) {
          if (!grid[cell].skip) {
            foreachstand(stand, s, grid[cell].standlist) {
              water_total += soilwater(&(stand->soil)) * stand->frac
                             * grid[cell].coord.area;
            }
            water_total += grid[cell].discharge.dmass_lake
                           + grid[cell].discharge.dmass_river;
          }
          runoff_out[cell] =
            grid[cell].discharge.dfout / 86400.0
            * (grid[cell].discharge.next < 0);
          grid[cell].discharge.delta = 0;
        }

#ifdef USE_MPI
        MPI_Allreduce(&water_total, &water_total_global, 1,
                      MPI_DOUBLE, MPI_SUM, config.comm);
#else
        water_total_global = water_total;
#endif
        if (isroot(config)) {
          printf("Total water mass in LPJmL day %d: %.12e l\n",
                 dayofyear, water_total_global);
        }

        dayofyear++;

        for (cell = 0; cell < config.ngridcell; cell++) {
          grid[cell].daily.tmax = -1e38;
          grid[cell].daily.tmin =  1e38;
          grid[cell].daily.prec = grid[cell].daily.swdown =
            grid[cell].daily.lwnet = grid[cell].daily.windspeed =
            grid[cell].daily.temp = grid[cell].VPD = 0;
        }
      } /* if (goodnight) */

      /* ---- monthly tasks ---- */
      if (monthend) {
        for (cell = 0; cell < config.ngridcell; cell++) {
          if (!grid[cell].skip) {
            update_monthly(grid + cell,
                           getmtemp(input.climate, &grid[cell].climbuf,
                                    cell, month),
                           getmprec(input.climate, &grid[cell].climbuf,
                                    cell, month),
                           month, &config);
          }
        }
        if (year >= config.outputyear) {
          fwriteoutput(output, grid, year, month, MONTHLY,
                       npft, ncft, &config);
        }
      } /* if (monthend) */

      /* ---- annual tasks (silvester) ---- */
      if (silvester) {
        for (cell = 0; cell < config.ngridcell; cell++) {
          if (!grid[cell].skip) {
            grid[cell].landcover =
              (config.prescribe_landcover != NO_LANDCOVER)
              ? getlandcover(input.landcover, cell) : NULL;
            update_annual(grid + cell, npft, ncft, year, TRUE,
                          intercrop, &config);
            check_fluxes(grid + cell, year, cell, &config);

            if (config.equilsoil) {
              if ((year - (config.firstyear - config.nspinup
                           + param.veg_equil_year
                           - param.equisoil_years))
                  % param.equisoil_interval == 0
                  && (year - (config.firstyear - config.nspinup
                              + param.veg_equil_year
                              - param.equisoil_years))
                     / param.equisoil_interval >= 0
                  && (year - (config.firstyear - config.nspinup
                              + param.veg_equil_year
                              - param.equisoil_years))
                     / param.equisoil_interval < param.nequilsoil) {
                equilveg(grid + cell, npft + ncft);
              }

              if (year == (config.firstyear - config.nspinup
                           + param.veg_equil_year)) {
                equilsom(grid + cell, npft + ncft,
                         config.pftpar, TRUE);
              }

              if ((year - (config.firstyear - config.nspinup
                           + param.veg_equil_year))
                  % param.equisoil_interval == 0
                  && (year - (config.firstyear - config.nspinup
                              + param.veg_equil_year))
                     / param.equisoil_interval > 0
                  && (year - (config.firstyear - config.nspinup
                              + param.veg_equil_year))
                     / param.equisoil_interval < param.nequilsoil) {
                equilsom(grid + cell, npft + ncft,
                         config.pftpar, FALSE);
              }

              if (param.equisoil_fadeout > 0) {
                if (year == (config.firstyear - config.nspinup
                             + param.veg_equil_year
                             + param.equisoil_interval
                             * param.nequilsoil)) {
                  equilveg(grid + cell, npft + ncft);
                }
                if (year == (config.firstyear - config.nspinup
                             + param.veg_equil_year
                             + param.equisoil_interval
                             * param.nequilsoil
                             + param.equisoil_fadeout)) {
                  equilsom(grid + cell, npft + ncft,
                           config.pftpar, FALSE);
                }
              }
            }

            if (config.withlanduse) {
              getnsoil_agr(&norg_soil_agr, &nmin_soil_agr,
                           &nveg_soil_agr, grid + cell);
              getoutput(&grid[cell].output,
                        DELTA_NORG_SOIL_AGR, &config) +=
                norg_soil_agr;
              getoutput(&grid[cell].output,
                        DELTA_NMIN_SOIL_AGR, &config) +=
                nmin_soil_agr;
              getoutput(&grid[cell].output,
                        DELTA_NVEG_SOIL_AGR, &config) +=
                nveg_soil_agr;
            }
          }

          grid[cell].balance.surface_storage =
            grid[cell].discharge.dmass_lake
            + grid[cell].discharge.dmass_river;
          grid[cell].aprec_poem_last = grid[cell].aprec_poem;
          if (grid[cell].ml.dam) {
            grid[cell].balance.surface_storage +=
              reservoir_surface_storage(grid[cell].ml.resdata);
          }
        }

        if (year >= config.outputyear) {
          fwriteoutput(output, grid, year, 0, ANNUAL,
                       npft, ncft, &config);
        }
      } /* if (silvester) */

    } /* end iterateyear_river() */

    /* ---- end-of-year flux accounting ---- */
    if (silvester) {
      if (year >= config.outputyear) {
        closeoutput_yearly(output, &config);
      }

      cflux_total = flux_sum(&flux, grid, &config);
      if (isroot(config)) {
        printflux(flux, cflux_total, year, &config);
        if (isopen(output, GLOBALFLUX)) {
          fprintcsvflux(output->files[GLOBALFLUX].fp.file, flux,
                        cflux_total, co2,
                        config.outnames[GLOBALFLUX].scale, year, &config);
        }
        if (output->files[GLOBALFLUX].issocket) {
          send_flux_coupler(&flux, config.outnames[GLOBALFLUX].scale,
                            year, &config);
        }
        fflush(stdout);
      }
      if (iswriterestart(&config) && year == config.restartyear) {
        fwriterestart(grid, npft, ncft, year,
                      config.write_restart_filename, FALSE, &config);
      }
    }

    yesterday = dayofmonth;

  } /* scoping block for time flags */

  fedisableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
} /* of 'lpj_update_' */

/*---------------------------------------------------------------------------*/
/* lpj_end_                                                       */
/*                                                                           */
/* Flush and close all output, free LPJmL memory.                           */
/*---------------------------------------------------------------------------*/
#ifdef __cplusplus
extern "C"
#endif
void lpj_end_(void) {
  fcloseoutput(output, &config);
  if (isroot(config)) {
    puts((year > config.lastyear)
         ? "Simulation ended." : "Simulation stopped.");
  }
  freeinput(input, &config);
  freegrid(grid, config.npft[GRASS] + config.npft[TREE], &config);
  if (isroot(config)) {
    printf((year > config.lastyear) ? "%s successfully" : "%s errorneously",
           progname);
    printf(" terminated, %d grid cells processed.\n", config.total);
  }
  freeconfig(&config);
} /* of 'lpj_end_' */

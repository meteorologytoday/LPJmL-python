import os
import sys

import numpy as np

import LPJmLCoupler

# Days per month — matches ndaymonth[] in LPJmL (no leap years).
NDAYMONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# Default values for FMS atmospheric coupling variables.
# These have no file-based equivalent in LPJmL and must be prescribed.
DEFAULT_T_FLUX  = 0.0       # sensible heat flux          [W m-2]
DEFAULT_QFLUX   = 0.0       # water-vapour flux            [kg m-2 s-1]
DEFAULT_DEDQ    = 1e-2      # dE/dq moisture sensitivity
DEFAULT_DHDT    = 1.0       # dH/dt                        [W m-2]
DEFAULT_DRDT    = 1.0       # dR/dt                        [W m-2]
DEFAULT_DRAG_Q  = 0.001     # drag coefficient for moisture [m s-1]
DEFAULT_P_SURF  = 101325.0  # surface pressure             [Pa]

def main():
    config = sys.argv[1] if len(sys.argv) > 1 else "lpjml_fms.cjson"

    # --- initialise -----------------------------------------------------------
    info = LPJmLCoupler.lpjml_init(config, dt_fast=3600)
    ncells    = info["ncells"]
    firstyear = info["firstyear"]
    lastyear  = info["lastyear"]
    nspinup   = info["nspinup"]
    rank      = info["rank"]
    if rank == 0:
        print(f"ncells={ncells}  years={firstyear - nspinup}..{lastyear}")

    # --- allocate input arrays ------------------------------------------------
    prec_in      = np.zeros(ncells)
    temp_mean_in = np.zeros(ncells) + 15.0
    swdown_in    = np.zeros(ncells) + 100.0
    lwnet_in     = np.zeros(ncells) + 100.0
    wind_in      = np.zeros(ncells) + 3.0
    co2_in       = np.zeros(ncells)

    t_flux_in    = np.full(ncells, DEFAULT_T_FLUX)
    qflux_in     = np.full(ncells, DEFAULT_QFLUX)
    dedq_in      = np.full(ncells, DEFAULT_DEDQ)
    dhdt_in      = np.full(ncells, DEFAULT_DHDT)
    drdt_in      = np.full(ncells, DEFAULT_DRDT)
    drag_q_in    = np.full(ncells, DEFAULT_DRAG_Q)
    p_surf_in    = np.full(ncells, DEFAULT_P_SURF)
    land_frac_in = np.ones(ncells)    # full land fraction
    q_ca_old_in  = np.zeros(ncells)   # updated each timestep from q_ca_out

    # --- allocate output arrays -----------------------------------------------
    carbon_flux_out      = np.zeros(ncells)
    q_ca_out             = np.zeros(ncells)
    runoff_out           = np.zeros(ncells)
    roughness_length_out = np.zeros(ncells)
    surface_temp_out     = np.zeros(ncells)
    albedo_out           = np.zeros(ncells)
    evap1_out            = np.zeros(ncells)
    dedt_out             = np.zeros(ncells)
    gc_out               = np.zeros(ncells)
    exlpj_out            = np.zeros(ncells)

    # --- time loop ------------------------------------------------------------
    for year in range(firstyear - nspinup, lastyear + 1):
        for month_idx, ndays in enumerate(NDAYMONTH):
            month = month_idx + 1        # 1-based
            for day in range(1, ndays + 1):
                for hour in range(24):

                    if rank == 0:
                        print(f"PYTHON LEVEL: {year:04d}-{month:02d}-{day:02d} {hour:02d}:00")
                    # Feed previous timestep's q_ca back as input.
                    np.copyto(q_ca_old_in, q_ca_out)

                    LPJmLCoupler.lpjml_update(
                        year=year, month=month, day=day,
                        hour=hour, minute=0, seconds=0,
                        prec_in=prec_in,
                        temp_mean_in=temp_mean_in,
                        swdown_in=swdown_in,
                        lwnet_in=lwnet_in,
                        t_flux_in=t_flux_in,
                        qflux_in=qflux_in,
                        dedq_in=dedq_in,
                        dhdt_in=dhdt_in,
                        drdt_in=drdt_in,
                        drag_q_in=drag_q_in,
                        p_surf_in=p_surf_in,
                        wind_in=wind_in,
                        co2_in=co2_in,
                        land_frac_in=land_frac_in,
                        q_ca_old_in=q_ca_old_in,
                        carbon_flux_out=carbon_flux_out,
                        q_ca_out=q_ca_out,
                        runoff_out=runoff_out,
                        roughness_length_out=roughness_length_out,
                        surface_temp_out=surface_temp_out,
                        albedo_out=albedo_out,
                        evap1_out=evap1_out,
                        dedt_out=dedt_out,
                        gc_out=gc_out,
                        exlpj_out=exlpj_out,
                    )

    if rank == 0:
        print("Calling lpjml_end")

    # --- finalise -------------------------------------------------------------
    LPJmLCoupler.lpjml_end()

if __name__ == "__main__":
    main()

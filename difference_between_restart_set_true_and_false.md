# Difference between `"restart" : true` and `"restart" : false` in LPJmL

Reference code: `lpjml_poem_59_standalone_driver`

---

## 1. Config parsing (`src/lpj/fscanconfig.c:1178`)

`"restart"` maps directly to `config->from_restart`:

```c
fscanbool2(file, &config->from_restart, "restart");

if (config->from_restart) {
    config->restart_filename = read_from_cjson(...);  // path to .lpj file
} else {
    if (config->sim_id == LPJML_FMS)
        return ERROR230;  // FMS sim_id requires a restart file
    config->restart_filename = NULL;
}
```

The key outcome: `config->restart_filename` is either a file path or `NULL`.

---

## 2. Grid/cell initialization (`src/lpj/newgrid.c:185`, `398`)

### `restart: false` — fresh start

```c
if (config->restart_filename == NULL && !config->ischeckpoint)
{
    file_restart = NULL;
    bfile = NULL;
    config->initsoiltemp = TRUE;   // soil temperatures computed from scratch
}
```

Then for each grid cell:

```c
if (config->restart_filename == NULL && !config->ischeckpoint)
{
    setseed(grid[i].seed, ...);
    grid[i].standlist = newlist(0);       // empty stand list
    grid[i].gdd       = newgdd(npft);     // zeroed growing degree days
    // zero sowing months and growing seasons
    addstand(&natural_stand, grid+i);     // one bare-ground natural stand
    initsoil(stand, soilpar, ...);        // soil pools set to zero/defaults
    new_climbuf(&grid[i].climbuf, ...);   // empty climate buffer
}
```

### `restart: true` — load from file

The block above is skipped entirely. Instead, `openrestart()` (or `bstruct_openrestart()` for the newer binary format) reads the full ecosystem state for each cell from the binary `.lpj` file:
- PFT vegetation state (biomass, leaf area, phenology, …)
- Soil carbon and nitrogen pools
- Soil temperatures and moisture
- Climate buffer (needed for fire, phenology)

---

## 3. Land use change gate (`src/driver/driver.c:490`)

```c
if (!config.natural_only
    && (year > config.firstyear - config.nspinup
        || config.from_restart))
{
    landusechange(...);
}
```

- `restart: false` → `landusechange()` is suppressed for the very first spinup year (`year == firstyear - nspinup`); active for all subsequent years.
- `restart: true` → `landusechange()` is called from year 0 onwards (the `|| config.from_restart` short-circuits).

(When `landuse: "no"` is set, `input.landuse == NULL` and the outer guard skips this entirely anyway.)

---

## 4. `equilsoil`

`equilsoil` is read separately from the config after the restart branch. In the spinup cjson it is controlled by a compile-time `#ifdef FROM_RESTART`:

```json
#ifdef FROM_RESTART
  "equilsoil" : false,   // already equilibrated — skip
#else
  "equilsoil" : true,    // needs soil carbon equilibration during spinup
#endif
```

When `restart: false` and `equilsoil: true`, LPJmL runs an additional soil carbon equilibration pass during the spinup to bring slow soil pools to steady state.

---

## Summary table

| Property | `restart: false` | `restart: true` |
|---|---|---|
| `config->restart_filename` | `NULL` | path to `.lpj` file |
| `config->initsoiltemp` | `TRUE` (computed from scratch) | `FALSE` (from file) |
| Vegetation state | Single empty natural stand (bare ground) | Full PFT state from restart file |
| Soil pools | Zero / defaults | Equilibrated values from restart file |
| Soil temperatures | Recomputed during spinup | Loaded from restart file |
| Climate buffer | Empty | Loaded from restart file |
| Land use at first spinup year | Suppressed | Active |
| Typical `equilsoil` | `true` | `false` |
| Typical `nspinup` | Large (e.g. 3000) | Small or 0 |

### When to use each

- **`restart: false`** — starting a brand-new simulation with no prior state. LPJmL needs a long spinup (`nspinup`) with `equilsoil: true` to build vegetation cover, soil carbon, and soil temperatures from bare ground.
- **`restart: true`** — continuing from a previously saved state. The ecosystem is already near equilibrium; a short or zero spinup is sufficient.

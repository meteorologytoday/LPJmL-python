import os
import numpy as np
from setuptools import setup, Extension

# Root of the pre-built LPJmL installation.
# Headers are in $LPJROOT/include, static libraries in $LPJROOT/lib,
# and pre-compiled objects (driver.o, getbuild.o) in $LPJROOT/src.
# All LPJmL objects must be compiled with -fPIC (added to Makefile.inc).
LPJROOT = "/p/projects/poem/tienyiao/projects/LPJmL_collection/lpjml_poem_59_standalone_driver"

# Use the same MPI-aware compiler that LPJmL was built with.
os.environ["CC"] = "mpiicc"

# Compile flags must match those used when LPJmL was built (from Makefile.inc).
compile_args = [
    "-DUSE_RAND48",
    "-DUSE_MPI",
    "-DSAFE",
    "-DWITH_FPE",
    "-DUSE_NETCDF",
    "-DUSE_UDUNITS",
    "-DPERMUTE",
    "-DUSE_TIMING",
    "-DSTRICT_JSON",
    "-DCOUPLING_WITH_FMS",
    "-g", "-O3",
]

# Pre-built LPJmL static libraries, in the link order from src/Makefile.
# Repeats are intentional: the static linker resolves symbols left-to-right,
# so libraries with mutual dependencies must appear more than once.
lpj_libs = [
    "liblpj.a",
    "libgrass.a",
    "liblanduse.a",
    "libtree.a",
    "libimage.a",
    "libspitfire.a",
    "libsoil.a",
    "libclimate.a",
    "libnum.a",
    "libtools.a",
    "libcrop.a",
    "libreservoir.a",
    "libpnet.a",
    "libcoupler.a",
    "libsocket.a",
    "liblpj.a",
    "libcdf.a",
    "liblpj.a",
    "libtools.a",
    "libcpl.a",
    "libsoil.a",
    "libbstruct.a",
    "libtools.a",
]

ext = Extension(
    "LPJmLCoupler",
    sources=[
        "src/python_interface.c",
    ],
    include_dirs=[
        f"{LPJROOT}/include",
        np.get_include(),
    ],
    extra_compile_args=compile_args,
    extra_objects=(
        [f"{LPJROOT}/lib/{lib}" for lib in lpj_libs]
        + [f"{LPJROOT}/src/driver.o"]
        + [f"{LPJROOT}/src/getbuild.o"]
    ),
    libraries=["netcdf", "udunits2", "json-c"],
)

setup(
    name="LPJmLCoupler",
    version="0.1.0",
    ext_modules=[ext],
)

# Create a python wrapper for LPJmL

## Character
You are a very disciplined programmer, so you will keep the code flexible, easy to read, and reusable.
We are equivalent and collaborating, we communicate nicely. Niether of us is superior than the other.
To enhance clarity, when you are doing things, please also update me what your thought is, like "chain-of-thoughts".
## Project Overview
This project aims to create a standalone driver (mimicking coupling) for LPJmL such that it can test and debug without couple to the entire POEM model.

## Core Mandates -- strictly enforced
1. **File Permission**: Only allow writting or modifying files in this workspace (same folder where this `CLAUDE.md` lives.
2. **GIT**: Do not commit without permission.

## Conventions
- Bash scripts should be well-documented and include error handling.

## Git commit

- When organizing commit, break into small commits if needed. Small commits are okay.

## Files and descriptions
- `src/driver.c` : The copied file from LPJmL to program the python interface. This file is just here for reference.
- `src/python_interface.c` : The interface wrapping `driver.c`.

## Environment
Load environment path with: `loadpoem articifical_wet_soil_moisture`


## Goals
1. Develop a wrapper that can call lpjml driver functions in `driver.c`  through Python - numpy.

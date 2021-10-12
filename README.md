# SPINNER
SPINNER(Structure Prediction of Inorganic crystals using Neural Network potentials with Evolutionary and Random searches)

If you use SPINNER, please cite this paper: S. Kang et al. Accelerated identification of equilibrium structures of multicomponent inorganic crystals using machine learning potentials (arXiv:2107.02594).

Here we describe minimal instructions for running SPINNER.
If you want more information such as tuning parameters, please visit our online manual(https://spinner-csp.readthedocs.io)

## Requirement
- LAMMPS 29Oct2020 or later
- randSpg
- SIMPLE-NN `81761d0` or later (if needed)

## Install
1. LAMMPS

``` bash
  cd src
  make mpi
  make mode=shlib mpi
```

2. randSpg

```bash
  mkdir build
  cd build
  cmake ..
  make -j3
```
 
3. `CMakeLists.txt`.
``` bash
SET ( LMP_PATH /path/to/lammps )
SET ( RANDSPG_PATH /path/to/randspg )

```

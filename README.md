# SPINNER
SPINNER(Structure Prediction of Inorganic crystals using Neural Network potentials with Evolutionary and Random searches)

If you use SPINNER, please cite this paper: S. Kang et al. Accelerated identification of equilibrium structures of multicomponent inorganic crystals using machine learning potentials (arXiv:2107.02594).

Here we describe minimal instructions for running SPINNER.
If you want more information such as tuning parameters, please visit our online manual(https://spinner-csp.readthedocs.io)

## Installation

### Download SPINNER
```
git clone https://github.com/MDIL-SNU/SPINNER.git
```

### Install randSpg

```
  mkdir build
  cd build
  cmake ..
  make -j3
```

```

### Install LAMMPS

```
  cd src
  make mpi
  make mode=shlib mpi
```

## Usage
To use SPINNER, 1 file (XXX.yaml) and 2 directories (input directory and src) are required.

### input file (XXX.yaml; XXX is named by the user)
Parameter list to control SPINNER code is listed in XXX.yaml. 
The simplest form of input.yaml is described below:
```YAML
# XXX.yaml

  input_dir:        input_directory_name
  output_dir:       output_directory_name
  initial_volume:   326.0
  structure:
      generation: 400
  material:
      Li:  4
      Hf:  2
      N:   4
```

### input directory
In input directory (input_dir: in input file), LAMMPS potential file should be located. (Potential file have to be named potential.) Potential should be SIMPLE-NN format (https://github.com/MDIL-SNU/SIMPLE-NN).

### src directory
src directory should be in the running directory. You can copy and paste the src directory to each running directory or you can run multiple calculations in one running folder.

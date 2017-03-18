[![DOI](https://zenodo.org/badge/85336350.svg)](https://zenodo.org/badge/latestdoi/85336350)
# Timeseries Tools
Fast, optimized tools for working with timeseries data. Calculates weighted frequency power spectra of un-evenly sampled data, filter and extracts frequencies. Written in Fortran 90 using OpenMP for parallelisation.

Written by Dr. Rasmus Handberg, Department of Physics and Astronomy, Aarhus University.


## Installation instructions:

### Using existing installation at SAC, Aarhus University

If you are and have access to the computer systems here, the easiest way to use the programs is simply to use the pre-compiled binaries available on the systems.
```
/usr/users/f041080/tidsserie/bin/machine
```
where machine is the system you are running on.

If you wish the program to use all available CPU cores available on your machine (and why wouldn't you?), you need to edit the file `~/.bashrc` or `~/.cshrc` (depending on which shell you use) and add the following line:
```
setenv OMP_NUM_THREADS 8 (csh, tcsh)
export OMP_NUM_THREADS=8 (bash)
```
The entire setup for the computer "ast", with 8 CPU cores, would for example look like this (in bash):
```
export PATH=$PATH:/usr/users/f041080/tidsserie/bin/ast
export OMP_NUM_THREADS=8
```
If you have any questions, feel free to ask me.


#### Matlab
```
setenv MATLABPATH $MATLABPATH:path-to-program/matlab
```

#### Python
```
setenv PYTHON_PATH $PYTHON_PATH:path-to-program/python
```

### Compiling it yourself

1. You should then be able to check out the code from the repo. using the following command:
```
git clone https://github.com/rhandberg/timeseries.git .
```

2. Run `make` followed by `make install`.
This will compile the program (by default using gfortran) and put the compiled binaries into the `bin` directory. If you wish to use a different compiler, you can simply modify the file `makefile.preample` where you can specify which compiler (and settings) to use.

3. Add the program to your path like described above.

## TODO
* Use Numerical Recipies trick in CalcAlphaBeta to speed up calculation.
* Make automatic or semi-automatic mode for windowfunction.
* Make oversampling-factor commandline input for spec and window.
* Finish writing install instructions.

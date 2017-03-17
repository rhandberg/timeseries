[![DOI](https://zenodo.org/badge/85336350.svg)](https://zenodo.org/badge/latestdoi/85336350)
# Timeseries Tools
Fast, optimized tools for working with timeseries data. Calculates weighted frequency power spectra of un-evenly sampled data, filter and extracts frequencies. Written in Fortran 90 using OpenMP for parallelisation.

Written by Dr. Rasmus Handberg, Department of Physics and Astronomy, Aarhus University.


## Installation instructions:
1. Run `make`
2. Run `make install`
3. Add the program to your path.

### Adding the program to your paths
```
set path=($path path-to-program/bin)
```

#### Matlab
```
setenv MATLABPATH $MATLABPATH:path-to-program/matlab
```

#### Python
```
setenv PYTHON_PATH $PYTHON_PATH:path-to-program/python
```


## TODO
* Use Numerical Recipies trick in CalcAlphaBeta to speed up calculation.
* Make automatic or semi-automatic mode for windowfunction.
* Make oversampling-factor commandline input for spec and window.
* Finish writing install instructions.

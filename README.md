# Bayes Factors for Model Selection in Multivariate Regression Designs

> linearReg_R2stat.m - a function to compute the Null-based Bayes Factor (BF)
> for a Multivariate Regression design. 
> Written by [Mert Türkol](mailto:mturkol_at_gmail_dot_com), (c) 2019.

linearReg_R2stat computes the BF (ratio between marginal likelihoods) by 
comparing a target regression model to the Null (intercept only). It can be used 
without having access to the full dataset, by utilizing only the ordinary R^2 
(coeff. of determination) test statistic obtained as a result of the underlying 
regression design.

Methodology: Laplace approximation to BF under Zellner-Siow prior 
(on model params) as a mixture of g-priors is utilized. That is, an inverse-gamma  
IG(1/2, s^2*N/2) prior on 'g', where 's' as in (0,1] is the scale factor 
hyper-parameter. Based on user choice, the numerical integration of the 
likelihood is carried out via Vectorized Adaptive Quadrature or 
High-Precision Numerical Integration using Variable-Precision Arithmetic.

This repository provides the MATLAB implementation of the methodology covered in
[1] & [2] as a package in R language, under the name BayesFactor. A few extra 
functionality that do not exist in vanilla BayesFactor package such as 
variable input tolerance(s) for integrator convergence and high-precision
numerical integration are included here.


## Table of Contents

<!-- vim-markdown-toc GFM -->

* [Usage](#usage)
    * [Example](#example)
    * [Inputs](#inputs)
    * [Options](#options)
* [Documentation](#documentation)
    * [Theoretical Background](#theoretical-background)
    * [File Overview](#file-overview)
* [FAQ](#faq)

<!-- vim-markdown-toc -->

## Install

### Requirements

Make sure you have `fftw3` installed.

### Compilation

Make use of the `Makefile` to compile the `lpsd-exec` by typing:

```
$ make
```

## Usage

`lpsd` can be controlled by command line options or interactively. 

```
$ ./lpsd-exec [OPTION...] 
```

Run `lpsd-exec` with the `help` option for a good overview of lpsd's
functionality and options:

```
$ ./lpsd-exec --help
```

### Example

Create a test dataset with:

```
$ make example-data
```

This will create the file `test.dat` in the directory `example` which can be
processed with `lpsd`:

```
$ ./lpsd-exec --input=example/test.dat
```
### Inputs
  (required, in order): 
   N       : #data-points/observations, (scalar int) N >= 3.
   p       : #predictors excluding the intercept, (scalar int) 1 <= p < N-1.
   R2      : Ordinary coeff. of determination, (real, scalar) 0 <= R2 < 1. 
             Corresponds to the proportion of variance accounted for 
             by the predictors, excluding the intercept.
### Options
  (in any order, use DEFAULT value if not provided as argument):
   's'     : prior scale (real, positive, scalar). 0 < s <= 1 
                         OR 
                         (char-array) as in {'medium','wide','ultrawide'}
                                              DEFAULT - 'medium' (0.3535).
   'useVpa': Logical to utilize (in Bf computation)
             (true)  vpaintegral() for High-Precision Numerical Integration 
                     using Variable-Precision Arithmetic; 
             (false) integral() for Vectorized Adaptive Quadrature - DEFAULT.                    
   'lvlTol': Level of tolerance (char-array) for convergence of integrator,
             as in {'Mdefault', 'Rdefault', 'medium', 'low', 'verylow'}
             'Mdefault' (relTol: 1e-6, absTol: 1e-10) - DEFAULT. 
                          default tolerances for "integral()" func. in MATLAB;
             'Rdefault' (relTol: eps('double')^0.25, eps('double')^0.25) 
                          default tolerances for "integrate()" func. in R;
             'medium'   (relTol: 1e-10, absTol: 1e-12);
             'low'      (relTol: 50*eps, absTol: 1e-14);
             'verylow'  (relTol: 5*eps, absTol: 1e-15).
   'relTol': Overwriting relative tolerance value (if) input by the user
             (real, scalar double) [0, Inf]. DEFAULT - 1e-6.
   'absTol': Overwriting absolute tolerance value (if) input by the user
             (real, scalar double) [0, Inf]. DEFAULT - 1e-10.
   'simple': Logical to return 
             (true)  the raw Bayes factor 'Bf'; 
             (false) log(Bf) in order to prevent possible overflow - DEFAULT.

## Documentation

### Theoretical Background

For a more theoretical understanding read the following publications:

  * [https://doi.org/10.1080/00273171.2012.734737](http://dx.doi.org/10.1080/00273171.2012.734737)
  * 

### File Overview

| File                 | Description                           |
| -------------------- | ------------------------------------- |
| `linearReg_R2stat.m` | Main function to compute Bayes Factor |
| `ask.c`       | Manages user input in interactive mode   |
| `calibrate.c` |                                          |
| `CHANGELOG`   | Changelog                                |
| `config.c`    | Configure lpsd at runtime via a textfile |
| `debug.c`     | Debugging                                |
| `errors.c`    | error messages                           |
| `genwin.c`    | compute window functions                 |
| `goodn.c`     |                                          |
| `IO.c`        | handle all input/output for `lpsd.c`     |
| `libargp.a`   | static library for argument parsing      |
| `lpsd`        | Executable                               |
| `lpsd.c`      |                                          |
| `lpsd.cfg`    | Configuration file                       |
| `lpsd-exec.c` |                                          |
| `Makefile`    | To build the executable                  |
| `misc.c`      |                                          |
| `netlibi0.c`  |                                          |
| `README.md`   | This README.md                           |
| `StrParser.c` |                                          |
| `tics.c`      |                                          |

## FAQ

* Hot can I enable `lpsd` configuration via config file?

To enable configuration of lpsd via `lpsd.cfg`, do one of the following under
linux, use the command

```
$ export LPSDCFN=/your/dir/lpsdconfigfilename
```

**Example:**

```
$ export LPSDCFN=/home/micha/lpsd.cfg
```

under DOS, use the command

```
$ set LPSDNCFN=c:\your\dir\lpsdconfigfilename
```

Example:

```
$ set LPSDCFN=c:\tools\lpsd.cfg
```

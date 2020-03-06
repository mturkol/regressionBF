# Bayes Factors for Model Selection in Multivariate Regression

> [linearReg_R2stat](/linearReg_R2stat.m) - a function to compute the Null-based 
> Bayes Factor (BF) for a Multivariate Regression design. 
> Written by [Mert Türkol](mailto:mturkol_at_gmail_dot_com), (c) 2019.

linearReg_R2stat computes the BF (ratio between marginal likelihoods) by 
comparing a target regression model to the Null (intercept only). It can be used 
without having access to the full dataset, by utilizing only the ordinary R^2 
(coeff. of determination) test statistic obtained as a result of the underlying 
regression design.

Methodology: Laplace approximation to BF under Zellner-Siow prior 
(on model params) as a mixture of g-priors is utilized. That is, an inverse-gamma 
IG(1/2, s^2*N/2) prior on 'g', where 's' as in (0,1] is the scale factor 
hyper-parameter. Based on user choice, numerical integration of the likelihood 
is carried out via Vectorized Adaptive Quadrature or High-Precision Numerical 
Integration using Variable-Precision Arithmetic.

This repository provides the MATLAB implementation of the methodology covered in
[1](http://dx.doi.org/10.1080/00273171.2012.734737) & [2] which resulted in a package in R-language, under the name BayesFactor. 
A few extra functionality that do not exist in vanilla BayesFactor package such 
as variable input tolerance(s) for integrator convergence and high-precision 
numerical integration are included here.


## Table of Contents

<!-- vim-markdown-toc GFM -->

* [Usage](#usage)
    * [Examples](#examples)
    * [Inputs](#inputs) 
    * [Options](#options)    
* [Documentation](#documentation)
    * [Theoretical Background](#theoretical-background)
    * [File Overview](#file-overview)

<!-- vim-markdown-toc -->

## Usage

Call `linearReg_R2stat` without input arguments to see its help documentation: 

```
>> linearReg_R2stat()
```

Call `linearReg_R2stat` with the the function handle '@getOptions' to return 
the default options setting within a Struct:

```
>> [DefOpts] = linearReg_R2stat(@getOptions)
```

Compute and return Bayes factor 'Bf' using 'Options' as a single Struct or 
name-value pairs passed as varargin: 

```
>> linearReg_R2stat(N, p, R2, varargin) 
```

Please see [linearReg_R2stat](/linearReg_R2stat.m) for a detailed list of input 
arguments and run Options. 

### Examples

Return the default Options in a Struct:

```
>> [DefOpts] = linearReg_R2stat(@getOptions) 
DefOpts = 
  struct with fields:

        s: 0.353553390593274
   useVpa: 0
   lvlTol: 'Mdefault'
   relTol: 1e-06
   absTol: 1e-10
   simple: 0 
```

Compute the Bayes Factor by reproducing the first entry of Table 1 \[[1](http://dx.doi.org/10.1080/00273171.2012.734737) \] : 
- Mandatory inputs: #observations (175), #predictors (4), R-squared (0.7109) ; 
- using a prior scale value of '1'; 
- utilizing the 'integral()' function for numerical integration with the  
  same relative and absolute tolerance values as in the built-in  
  'integrate()' function found in R;  
- to return the raw (simple) Bf value along with the Options of computation 

```
>> [rawBf, Opts] = linearReg_R2stat(175, 4, 0.7109, ...
                   's', 1, 'lvlTol', 'Rdefault', 'simple', true)
                   
rawBf =
        3.54443032945501e+41
Opts = 
  struct with fields:

        s: 1
   useVpa: 0
   lvlTol: 'Rdefault'
   relTol: 0.0001220703125
   absTol: 0.0001220703125
   simple: 1 
```

Alternatively, following the mandatory inputs, we could pass a Struct for the 
desired run Options ( i.e. using Model #6 with the covariates 'Local+Parasites' as in 
[1](http://dx.doi.org/10.1080/00273171.2012.734737) ): 
```
>> [Bf60, RouderOpts] = linearReg_R2stat(175, 2, 0.2429, ...
                        struct('s', 1, 'lvlTol', 'Rdefault', 'simple', true) )
                   
Bf60 =
          122607538.194857
RouderOpts = 
  struct with fields:

         s: 1
    useVpa: 0
    lvlTol: 'Rdefault'
    relTol: 0.0001220703125
    absTol: 0.0001220703125
    simple: 1 
```

### Inputs
  (required, in order): 
  | Input    | Description                                                      |
  | -------- | ---------------------------------------------------------------- | 
  | `N`      | #data-points/observations, (scalar int) N >= 3.                  |                
  | `p`      | #predictors excluding the intercept, (scalar int) 1 <= p < N-1.  |
  | `R2`     | Ordinary coeff. of determination, (real, scalar) 0 <= R2 < 1.    |    
  |          | Corresponds to the proportion of variance accounted for by the predictors, excluding the intercept. | 

### Options
  (in any order, use DEFAULT value if not provided as argument):
  | Parameter    | Description, Default Value. Type/Value Tags                                            |
  | ------------ | ---------------------------------------------------------------------------------------|
  | `'s'`        | prior scale factor, DEFAULT - 'medium'.                                                |
  |              | (real, positive, scalar) (0, 1]                                                  |
  |              |          OR                                                                            |
  |              | (char-array) as in {'medium','wide','ultrawide'}                                      |
  | `'useVpa'`   | Logical to utilize Variable-Precision Arithmetic in Bf computation, DEFAULT - 'false'. |
  |              | true  - vpaintegral() for High-Precision Numerical Integration                        |
  |              | false - integral() for Vectorized Adaptive Quadrature                                 |                   
  | `'lvlTol'`   | Level of tolerance for convergence of integrator, DEFAULT - 'Mdefault'.                |
  |              | (char-array) as in {'Mdefault', 'Rdefault', 'medium', 'low', 'verylow'}                |
  |              | 'Mdefault' - (relTol: 1e-6, absTol: 1e-10)                                              |
  |              | 'Rdefault' - (relTol: eps('double')^0.25, absTol: eps('double')^0.25)                   |
  |              | 'medium'   - (relTol: 1e-10, absTol: 1e-12)                                             |
  |              | 'low'      - (relTol: 50*eps, absTol: 1e-14)                                            |
  |              | 'verylow'  - (relTol: 5*eps, absTol: 1e-15)                                             |
  | `'relTol'`   | Relative tolerance value, DEFAULT - 1e-6.                                              |
  |              | (real, scalar double) [0, Inf]                                                        |
  | `'absTol'`   | Absolute tolerance value, DEFAULT - 1e-10.                                             |
  |              | (real, scalar double) [0, Inf]                                                        |
  | `'simple'`   | Logical to return the Bf in raw form, DEFAULT - false.                                 | 
  |              | true  -  the raw Bayes factor 'Bf'                                                     |
  |              | false -  log(Bf) in order to prevent possible overflow;                                 |

## Documentation

### Theoretical Background

For a more theoretical understanding read the following publications:

  * \[1\] [https://doi.org/10.1080/00273171.2012.734737](http://dx.doi.org/10.1080/00273171.2012.734737)
  * \[2\] [10.1198/016214507000001337Liang](https://www.tandfonline.com/doi/abs/10.1198/016214507000001337)

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

function [varargout] = linearReg_R2stat(varargin)
% Name       : linearReg_R2stat.m
% Version    : 1.0
% Author     : Mert Turkol
% Date       : 05/08/2019
% Copyright  : GNU General Public License v2.0 (GPLv2)
% Description: This function computes the Bayes factor (Bf) against the null 
%              hypothesis in a multiple regression design, using the ordinary 
%              R^2 (coefficient of determination) statistic. It can be used 
%              without having access to the full dataset, by utilizing only the 
%              test statistic obtained as a result of the underlying regression 
%              design.
% Requires   : logUtility.m, dinvgamma.m
% 
% For details about the model under consideration, see the references cited at
% the end of this help documentation.
% 
% Syntax : 
%                     linearReg_R2stat() :
%                       Without input arguments, display the help documentation
%
%   [Options]       = linearReg_R2stat(@getOptions) :
%                       Return only the default settings for options (Struct)
%
%   [bf10, Options] = linearReg_R2stat(N, p, R2, varargin) :
%                       Compute and return Bayes factor 'Bf', and the 'Options' 
%                       used (see OPTIONS section below)
% ==============================================================================
%   INPUTS (required, in order): 
%    N       : #data-points/observations, (scalar int) N >= 3.
%    p       : #predictors excluding the intercept, (scalar int) 1 <= p < N-1.
%    R2      : Ordinary coeff. of determination, (real, scalar) 0 <= R2 < 1. 
%              Corresponds to the proportion of variance accounted for 
%              by the predictors, excluding the intercept.
%
%   OPTIONS (in any order, use DEFAULT value if not provided as argument): 
%    's'     : prior scale (real, positive, scalar). 0 < s <= 1 
%                          OR 
%                          (char-array) as in {'medium','wide','ultrawide'}
%                                               DEFAULT - 'medium' (0.3535).
%    'useVpa': Logical to utilize (in Bf computation)
%              (true)  vpaintegral() for High-Precision Numerical Integration 
%                      using Variable-Precision Arithmetic; 
%              (false) integral() for Vectorized Adaptive Quadrature - DEFAULT.                    
%    'lvlTol': Level of tolerance (char-array) for convergence of integrator,
%              as in {'Mdefault', 'Rdefault', 'medium', 'low', 'verylow'}
%              'Mdefault' (relTol: 1e-6, absTol: 1e-10) - DEFAULT. 
%                           default tolerances for "integral()" func. in MATLAB;
%              'Rdefault' (relTol: eps('double')^0.25, eps('double')^0.25) 
%                           default tolerances for "integrate()" func. in R;
%              'medium'   (relTol: 1e-10, absTol: 1e-12);
%              'low'      (relTol: 50*eps, absTol: 1e-14);
%              'verylow'  (relTol: 5*eps, absTol: 1e-15).
%    'relTol': Overwriting relative tolerance value (if) input by the user
%              (real, scalar double) [0, Inf]. DEFAULT - 1e-6.
%    'absTol': Overwriting absolute tolerance value (if) input by the user
%              (real, scalar double) [0, Inf]. DEFAULT - 1e-10.
%    'simple': Logical to return 
%              (true)  the raw Bayes factor 'Bf'; 
%              (false) log(Bf) in order to prevent possible overflow - DEFAULT.
%
%   OUTPUT :
%    bf10    : Bayes factor against the intercept-only null(real, scalar double)
%    Options : Struct with fieldnames corresponding to the OPTIONS above
%   
%   EXAMPLE USAGE:
%   Return the default Options within a Struct :
%
%   >> [DefOpts] = linearReg_R2stat(@getOptions) 
%   DefOpts = 
%     struct with fields:
%
%           s: 0.353553390593274
%      useVpa: 0
%      lvlTol: 'Mdefault'
%      relTol: 1e-06
%      absTol: 1e-10
%      simple: 0 
%
%   Compute the Bayes Factor by imitating the first entry of Table 1 [1] :
%   - #observations: 175, #predictors: 4, R-squared value: 0.7109;
%   - using a prior scale value of '1';
%   - utilizing the 'integral()' function for numerical integration with the 
%     same relative and absolute tolerance values as in the built-in 
%     'integrate()' function found in R;  
%   - to return the raw (simple) Bf value along with the Options of computation
% 
%   >> [rawBf, Opts] = linearReg_R2stat(175, 4, 0.7109, ...
%                      's', 1, 'lvlTol', 'Rdefault', 'simple', true)
%   rawBf =
%           3.54443032945501e+41
%   Opts = 
%     struct with fields:
%
%           s: 1
%      useVpa: 0
%      lvlTol: 'Rdefault'
%      relTol: 0.0001220703125
%      absTol: 0.0001220703125
%      simple: 1 
% ==============================================================================
% Comments: The algorithm behind 'linearReg_R2stat.m' was developed based on 
%           /BayesFactor_0.9.12-4.2/BayesFactor/R/linearReg_R2stat.R
% References: 
%   [1] Jeffrey N. Rouder & Richard D. Morey (2012). Default Bayes Factors 
%       for Model Selection in Regression, Multivariate Behavioral Research, 
%       47:6, 877-903
%       Link: http://dx.doi.org/10.1080/00273171.2012.734737
%   [2] Feng Liang, Rui Paulo, German Molina, Merlise A Clyde & Jim O Berger
%       (2008) Mixtures of g Priors for Bayesian Variable Selection, Journal  
%       of the American Statistical Association, 103:481, 410-423, 
%       DOI: 10.1198/016214507000001337Liang  
%       Link: https://www.tandfonline.com/doi/abs/10.1198/016214507000001337
%   [3] Morey, R. D., Rouder, J. N., Jamil, T., & Morey, M. R. D. (2015). 
%       Package 'bayesfactor'. 
%       URL: http://cran/r-projectorg/web/packages/BayesFactor/BayesFactor.pdf
% ============================================================================== 
if isempty(varargin), help linearReg_R2stat; return; end
narginchk(1, Inf);

if isa(varargin{1}, 'function_handle')
  VALID.OPER_LIST = {'getOptions'};
  
  funH = functions(varargin{1});
  funName = funH.function;
  switch lower(funName)

    case 'getoptions'
      
      if length(varargin) > 1
        funIn = varargin(2:end);
      else
        funIn = {};
      end
    funH = @parseFunInput;
    funArgOut = 1;
    
    otherwise
     
      error('linearReg_R2stat:InvalidFunctionCall',...
      ['''@%s'' is not a valid functionality within %s.\n', ...
       'Available functions to call with a function handle are: \n',...
       repmat('@%s\n', 1, numel(VALID.OPER_LIST) )], ...
       funName, mfilename, VALID.OPER_LIST{:});
  end
  
  nargoutchk(0,funArgOut);
  varargout{funArgOut} = funH(funIn); 
  return
end

nargoutchk(0,2);
if length(varargin) >= 3
  N = varargin{1};
  p = varargin{2};
  R2 = varargin{3};
  if length(varargin) > 3
    varargin = varargin(4:end);
  else
    varargin = {};
  end
else
  error('linearReg_R2stat:NotEnoughInput', ...
        ['Not enough input!\nTo compute a Bayes factor, %s ', ...
         'requires at least 3 inputs: (N, p, R2). User put in only %s.'], ...
         mfilename, num2str( length(varargin) ) );  
end

validNp = @(num) isscalar(num) && isPositiveIntegerValuedNumeric(num); 
validR2 = @(r2) isscalar(r2) && isnumeric(r2) && isreal(r2) && ...
          (r2 >= 0) && (r2 < 1);
assert( (N >=3) && validNp(N), ...
      'linearReg_R2stat:IllegalNoSample', ...
      '#data-points "N" must be a scalar integer greater than or equal to 3.'); 
assert( ( p < (N-1) ) && validNp(p), ...
      'linearReg_R2stat:IllegalNoPredictor', ...
      '#predictors "p" must be a scalar integer: 1 <= p < N-1.(N: #data-points)');
assert( validR2(R2), ...
      'linearReg_R2stat:IllegalR2value', ...
      'Ordy. coeff. of determination "R2" must be a real scalar within [0,1)');

if ~isempty(varargin) && cellfun(@isstruct,varargin(1))
  varargin = varargin{1};
elseif ~isempty(varargin) && ~cellfun(@isstruct, varargin(1) )
  varargin = varargin( ~cellfun(@isstruct, varargin ) );
end
% Call input-parser nested func. to parse the inputs of linearReg_R2stat
[Options, pInput] = parseFunInput(varargin);

% Extract the analysis parameters from the parsed input
s = Options.s; %scale factor
relTol = Options.relTol; %limit of rel. tolerance for integration convergence
absTol = Options.absTol; %limit of abs. tolerance for integration convergence
%name of the function to call for numerical integration
if ~Options.useVpa
  method = 'integral';
else
  method = 'vpaintegral';
end
simple = Options.simple; %false: return log(bf10), true: return exp(log(bf10))
%% Compute approximation to posterior mode of g
% Liang et al Eq. A.3, assuming a=b=0. For details, see [2]. 
g3 = ( R2 - 1 ) * ( p + 3 ); % g^3
% The following expression is algebraically equal to (R2 - 1)*(p+3).
% The advantage in finite precision arithmetic is that it avoids subtracting 
% nearly equal numbers.
if ( g3 == 0 )
	g3 = R2 * (0.5 + (0.5 - 1.0 / R2)) * ( p + 3);
end
% g2 = ( N - p - 4 - 2 .* (1 - R2) ); % g^2
g2 = N + 2 * R2 - (p + 6); % g^2
g1 = N * (2 - R2) - 3 ; % g
g0 = N;

gPolyCoeff = [ g3 g2 g1 g0 ]; %coefficients of the cubic polynomial
sol = roots(gPolyCoeff); %  sol = polyroot(c(g0, g1, g2, g3))
% Find and pick the real solution
[~, minImSqrIdx] = min(imag(sol).^2);
modeg = real(sol(minImSqrIdx)); % modeg = Re(sol[which.min(Im(sol)^2)])

if(modeg <= 0)
  solNew = [];
  if ~verLessThan('symbolic', 'R2012b')
  % Try vpasolve to to find a positive root
    k = []; syms k;  
    gPolyVpa = g3.*k.^3 + g2.*k.^2 + g1.*k + g0;
    %epsZero = eps(0);
    solNew = vpasolve(gPolyVpa == 0, k, [realmin inf]);
  end
  
  if ~isempty(solNew) && isreal(solNew) && double(solNew)>0 
    modeg = double(solNew); 
  else
    modeg = N ./ 20;
  end 
end
%% Integrate the likelihood w.r.t. to priors to obtain Bayes Factor
% Utilizes Laplace approximation to BF under Zellner-Siow prior as a mixture of 
% g-priors. [2]. That is, an inverse-gamma (1/2, s^2*N/2) prior on g.
% The original likelihood (integrand) expression (w/ scale factor 's')
% is in Eqn. 11 [1].

logConst = integrandRegressionU( 0, s.^2, true, 0, log(modeg) );
% log const. definiton within the original R-code as found in
% /BayesFactor_0.9.12-4.2/BayesFactor/R/linearReg_R2stat.R (line 82):
% log.const = integrand.regression.u(0, N, p , R2, rscaleSqr=rscale^2, ...
%                                                  log=TRUE, shift=log(modeg))
funHandle = (@(U) integrandRegressionU( U, s.^2, ...
                                       false, logConst, log(modeg) ) );                                        
if strcmp(method, 'integral')% use MATLAB's integral() to integrate    
  h = integral(funHandle, -Inf, +Inf, 'RelTol', relTol, 'AbsTol', absTol);

else % use MATLAB's vpaintegral() to integrate
  %(if MATLAB's toolbox version allows - checked earlier by parseFunInput)
  %G = []; syms G;
  U = []; syms U;
  
  %logConst = integrandRegression( log(modeg), N, p , R2, s.^2, true, 0 );
  %funHandle = (@(G) integrandRegression( G, N, p, R2, s.^2, false, logConst) );
  
  %H = symfun(funHandle(G), G);
  H = symfun(funHandle(U), U);
  %h = vpaintegral(H, G, 0, +Inf, 'RelTol', relTol, 'AbsTol', absTol, ...
  %                               'MaxFunctionCalls', Inf ); 
  h = vpaintegral(H, U, -Inf, +Inf, 'RelTol', relTol, 'AbsTol', absTol, ...
                                 'MaxFunctionCalls', Inf);
end
  % Integration function call as found in the Bayes Factor R-Package, 
  % /BayesFactor_0.9.12-4.2/BayesFactor/R/linearReg_R2stat.R (line 48):
  %
  % h=integrate(integrand.regression.u,lower=-Inf,upper=Inf,N=N,p=p,R2=R2, ...
  %             rscaleSqr=rscale^2,log.const=log.const,shift=log(modeg))
bf10 = log(h) + logConst; % := log(bf_M0)
if simple
  if bf10 > log(realmax)
    warning('linearReg_R2stat:LargerThanRealMax', ...
            ['The computed value of log(bf) is larger than log(realmax). ', ...
            'Taking its exponential will return Inf!']);
  end
  bf10 = exp(bf10);
end
if strcmp(method, 'vpaintegral')
  bf10 = double(bf10);
end
varargout{1} = bf10;
varargout{2} = Options;
%% Nested Function: integrandRegressionU
% Computes the integrand using a Laplace approximation to BF under Zellner-Siow 
% prior (Appendix A [2]). Vectorized so it can be called as a func. handle or 
% symbolic func. Claimed in the original BF package to be numerically more 
% stable than integrandRegression.
  function [I] = integrandRegressionU( u, rscaleSqr, returnLog, logConst, shift)
    narginchk(1, 5); %narginchk(4, 8);
    if nargin < 2, rscaleSqr = 1; end %was nargin < 5
    if nargin < 3, returnLog = false; end %was nargin < 6
    if nargin < 4, logConst = 0; end %was nargin < 7
    if nargin < 5, shift = 0; end %was nargin < 8
    u = u + shift; % log transformation -> u := log(g)
    
    if ~isa(u, 'sym')
      % utilize logUtility functions to handle possible precision issues
      log1pExpTerm1 = arrayfun(@(logUtilX1) ...
                      logUtility('log1pExp', logUtilX1), u);
      log1pExpTerm2 = arrayfun(@(logUtilX2) ...
                      logUtility('log1pExp', logUtilX2), u + log(1.0 - R2) );
    else 
      log1pExpTerm1 = log( 1.0 + exp(u) );
      log1pExpTerm2 = log( 1.0 + exp( u + log(1-R2) ) );
    end
    a = .5 .* ( (N - p - 1 ) .* ...
        log1pExpTerm1 - (N - 1) .* ...
        log1pExpTerm2 );
    % Equivalent expression within the vectorized func. in BayesFactor R                  
    % a = .5 .* ( (N - p - 1 ) .* ...
    %     log1pExp( u ) - (N - 1) .* ...
    %     log1pExp( u + log(1 - R2) ) );
    
    % Inverse-gamma (shape of 1/2, scale-param of s^2/2) mixture of g's
    shape = .5;
    scale = rscaleSqr .* N ./ 2;
    % log density of the inv. gamma distro, as a function of log(u)
    logDensityIgam = dinvgamma(u, shape, scale, 0, true, true);
    logI = a + logDensityIgam - logConst + u; %log of the integral(likHood)
    if returnLog
      I = logI;
    else
      I = exp(logI);
    end
  end %end-of-nested func.: integrandRegressionU
  %% The original "integrand.regression.u" function definition in R
  % as found in Bayes Factor (0.9.12-4.2):
  % # This is a more numerically stable version of the integrand, as a f(log(g))
  % integrand.regression.u=Vectorize(...
  %                        function(u, N, p, R2, rscaleSqr=1, log=FALSE, ...
  %                        log.const=0, shift = 0){
  % u = u + shift
  % a = .5 * ((N - p - 1 ) * log1pExp(u) -
  %          (N - 1) * log1pExp(u + log(1 - R2)))
  % shape=.5
  % scale=rscaleSqr*N/2
  % log.density.igam <- dinvgamma(u, shape, scale, log=TRUE, logx=TRUE)
  % ans = a + log.density.igam - log.const + u
  % ifelse(log,ans,exp(ans))
  % },"u")
%% Nested Function: integrandRegression
  function [I] = integrandRegression( g, rscaleSqr, returnLog, logConst)
    narginchk(1, 4); %narginchk(4, 7);
    if nargin < 2, rscaleSqr = 1; end %was nargin < 5
    if nargin < 3, returnLog = false; end %was nargin < 6
    if nargin < 4, logConst = 0; end %was nargin < 7
    
    if ~isa(g, 'sym')
      a = .5 .* ( (N - p - 1 ) .* log1p(g) ...
                - (N - 1) .* log1p( g.* (1 - R2) ) );
    else
      a = .5 .* ( (N - p - 1 ) .* log(1 + g) ...
          - (N - 1) .* log(1 + g.* (1 - R2) ) );
    end
    % Inverse-gamma (shape of 1/2, scale-param of s^2/2) mixture of g's
    shape = .5;
    scale = rscaleSqr .* N ./ 2;
    % log density of the inv. gamma distro, as a function of log(u)
    logDensityIgam = dinvgamma(g, shape, scale, 0, true);
    logI = a + logDensityIgam - logConst; %log of the integral
    if returnLog
      I = logI;
    else
      I = exp(logI);
    end 
  end %end-of-nested func.: integrandRegression
  %% The original "integrand.regression" function definition in BF R Package
  % integrand.regression=Vectorize(function(g, N, p, R2, rscaleSqr=1, ...
  %                                         log=FALSE, log.const=0){
  % a = .5 * ((N - p - 1 ) * log(1 + g) - (N - 1) * log(1 + g * (1 - R2)))
  % shape=.5
  % scale=rscaleSqr*N/2
  % log.density.igam <- dinvgamma(g, shape, scale, log=TRUE)
  % ans = a + log.density.igam - log.const
  % ifelse(log,ans,exp(ans))
  % },"g")
%% Nested Function: parseFunInput
  % Parses the input arguments of the caller func. (linearReg_R2stat)
  function [pOptions, pInput] = parseFunInput(funVarargin)
  DEFAULT.S = 'medium';
  DEFAULT.USE_VPA_INTEGRAL = false;
  DEFAULT.LVL_TOL = 'Mdefault'; %MATLAB's default tolerances for integrator
  DEFAULT.REL_TOL = 1e-6;
  DEFAULT.ABS_TOL = 1e-10;
  DEFAULT.RETURN_SIMPLE = false;
  
  % Convert any string-type varargin to char-array (leave non-strings untouched)
  try funVarargin=controllib.internal.util.hString2Char(funVarargin); catch,end
  pInput = inputParser;
  %% Input Validation funcs
  isLogical = @(x) isscalar(x) && islogical(x);
  isCharVec = @(x) ischar(x) && ( size(x, 1) == 1 ) && isvector(x);
  validS = @(scale) isCharVec(scale) || ...
                 ( ~isCharVec(scale) && isscalar(scale) && isreal(scale) && ...
                 isnumeric(scale) && (scale > 0) && (scale <= 1) ); 
  validUseVpa = @(useVpa) isLogical(useVpa); 


  % Use inline-if anonymous function to add parameters based on MATLAB version
  iifParser  = @(varargin) varargin{2*find([varargin{1:2:end}], 1,'first')}();
  addParamFun = @(param, paramName, paramValidator) iifParser...
                ( ~verLessThan('matlab', '8.2.0.29'),  ...
                  @() addParameter( param, paramName, paramValidator{:} ),...
                ( true ), ...
                  @() addParamValue( param, paramName, paramValidator{:} ) );
                   
  addParamFun(pInput, 's', {DEFAULT.S, validS});
  addParamFun(pInput, 'useVpa', {DEFAULT.USE_VPA_INTEGRAL, validUseVpa});
  addParamFun(pInput, 'lvlTol', {DEFAULT.LVL_TOL, isCharVec});
  % default rel. and abs. tol's of the integrator are evaluated automatically.
  addParamFun(pInput, 'relTol', {DEFAULT.REL_TOL});
  addParamFun(pInput, 'absTol', {DEFAULT.ABS_TOL});
  addParamFun(pInput, 'simple', {DEFAULT.RETURN_SIMPLE, isLogical});
  %% Parse the inputs and carry out further input-error checking
  if isstruct(funVarargin)
    parse(pInput, funVarargin); %parse(pInput, N, p, R2, funVarargin);
  else    
    parse(pInput, funVarargin{:}); %parse(pInput, N, p, R2, funVarargin{:}); 
  end
  %% Evaluate parameter-input values passed by the user. 
  % Override necessary ones and assign into pOptions struct
  pOptions = struct('s', pInput.Results.s, ...
                   'useVpa', [], ...
                   'lvlTol', [], ...
                   'relTol', [], ...
                   'absTol', [], ...
                   'simple', pInput.Results.simple);
    if isCharVec(pOptions.s)
      expectedPriorType = {'medium', 'wide', 'ultrawide'};
      try
        pOptions.s = validatestring(pOptions.s, ...
                                   expectedPriorType, ...
                                   mfilename, 's');
      catch
        pOptions.s = DEFAULT.S;
      end
      %  For the 's' (rscaleCont) argument, several named values are recongized: 
      %  "medium" (DEF), "wide", and "ultrawide", which correspond r scales of 
      %  sqrt(2)/4, 1/2, and sqrt(2)/2, respectively. These values were chosen
      %  to yield consistent Bayes factors with anovaBF. 
      %  Ref: BayesFactor Version 0.9.12-4.2 manual,'regressionBF',pp. 47, [3] 
      
        switch pOptions.s
          case 'medium'
            pOptions.s = sqrt(2)/4;
          case 'wide'
            pOptions.s = 1/2;
          case 'ultrawide'
            pOptions.s = sqrt(2)/2;
        end
    end
    
    if pInput.Results.useVpa && ~verLessThan('symbolic', 'R2016b')
      pOptions.useVpa = true;
    else
      pOptions.useVpa = false;
    end 
    
    expectedLvlTol = {'Mdefault', 'Rdefault', 'medium', 'low', 'verylow'};
    try
     pOptions.lvlTol = validatestring(pInput.Results.lvlTol, ...
                                     expectedLvlTol, ...
                                     mfilename, 'lvlTol');
    catch
      pOptions.lvlTol = DEFAULT.LVL_TOL;
    end
    
    switch pOptions.lvlTol
      case 'Mdefault' % Default tolerances in MATLAB's integral(),vpaintegral()
        pOptions.relTol = 1e-6;
        pOptions.absTol = 1e-10;
      case 'Rdefault' % Default tolerances of "integrate()" function in R
        pOptions.relTol = eps('double')^0.25;
        pOptions.absTol = eps('double')^0.25;        
      case 'medium'
        pOptions.relTol = 1e-10;
        pOptions.absTol = 1e-12;
      case 'low'
        pOptions.relTol = 50*eps; % 1.11022302462516e-14
        pOptions.absTol = 1e-14;
      case 'verylow'
        pOptions.relTol = 5*eps; 
        pOptions.absTol = 1e-15;
    end
    % Overwrite relTol and absTol with the values (if) input by the user 
    if ~ismember('relTol', pInput.UsingDefaults)
      pOptions.relTol = pInput.Results.relTol;
    end
    if ~ismember('absTol', pInput.UsingDefaults)
      pOptions.absTol = pInput.Results.absTol;
    end   
  end
end
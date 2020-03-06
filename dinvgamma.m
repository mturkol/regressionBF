function [derLogDensityIgam] = dinvgamma(x, shape, ...
                                     scale, derOrder, isLog, isLogX)
% Name       : dinvgamma.m
% Version    : 1.0
% Author     : Mert Turkol
% Date       : 05/08/2019
% Copyright  : GNU General Public License v2.0 (GPLv2)
% Description: This function computes and returns the derivative of a desired 
%              order of the log density of an inverse gamma (IG) distribution  
%              parametrized by 'shape' (a) and 'scale' (b), i.e. w/ the density
%                f_IG(x; a, b) = ( b^a ) / ( gamma(a) * x^(a+1) ) * exp( -b/x )
%              The desired derivative is computed w.r.t to 'x' or 'log(x)'. 
%
%              The implementation presented here is fully vectorized so that the
%              output can be automatically obtained for an array of the same 
%              sized vector input 'x', using the same (scalar) arguments 
%              within the rest of the function input list.
% Requires   : lngamma.m [2]
% TO DO      : check and rework 'd2dinvgamma1_logx' 
% Syntax     : 
%   dinvgamma() :
%     Without input arguments, display the help documentation
%
%   [derLogDensityIgam] = dinvgamma(x, shape) :
%     Compute as a function of 'x', the (exponential of) log density of 
%     IG(x; shape, scale = 1), .  
%
%   [derLogDensityIgam] = dinvgamma(x, shape, scale) :
%     Compute as a function of 'x', the (exponential of) log density of  
%     IG(x; shape, scale).  
%
%   [derLogDensityIgam] = dinvgamma(x, shape, scale, derOrder) :
%     Compute as a function of 'x', the 'derOrder'th derivative of the 
%     (exponential of) log density of IG(x; shape, scale). 
%     Zeroth order (derOrder = 0) corresponds to computing the density itself.
%
%   [derLogDensityIgam] = dinvgamma(x, shape, scale, derOrder, isLog) :
%     Compute as a function of 'x', the 'derOrder'th derivative of log density 
%     of IG(x; shape, scale) IF isLog:=true; 
%                            otherwise its exponential (as in the above).
%
%   [derLogDensityIgam] = dinvgamma(x, shape, scale, derOrder, isLog, isLogX) :
%     Same as above, except the derivative is taken w.r.t. 
%       'log(x)', IF isLogX:=true,
%       'x'     , otherwise (isLogX:=false).
% ==============================================================================
% Comments:
%   Main portion of the MATLAB-code implemented here (until the definition of 
%   nested func's) is based on /BayesFactor_0.9.12-4.2/BayesFactor/R/common.R 
%   (starting at line 270). For comparison, the original R code from the 
%   BayesFactor package [1] is given at the end.
% 
%   The MATLAB code implemented here utilizes Paul Godfrey's "gammaln" function. 
%   It was renamed to "lngamma" in this package in order to prevent the vanilla 
%   a potential clash with the "gammaln" function that comes as default with 
%   MATLAB distribution. Paul Godfrey's original routine (within the 
%   "Special Functions math library" on MATLAB File Exchange) implements 
%   "an excellent Lanczos series approximation for the complex ln(Gamma) function"
%   https://www.mathworks.com/matlabcentral/fileexchange/978-special-functions-math-library)
%   For a description of the algorithm implemented there, see [2].
% References: 
%   [1] Morey, R. D., Rouder, J. N., Jamil, T., & Morey, M. R. D. (2015). 
%       Package 'bayesfactor'. 
%       http://cran/r-projectorg/web/packages/BayesFactor/BayesFactor.pdf
%   [2] Godfrey, P. (2001). 
%       A note on the computation of the convergent Lanczos complex Gamma approximation.
%       http://my.fit.edu/~gabdo/gamma.txt
% ==============================================================================
if nargin < 1, help dinvgamma; return; end 
narginchk(2, 6); %at least 'x' and 'shape' are required
MAX_DER = 2; %currently up to 2nd order derivative have been implemented
%% Initial input-checking  
isValidDer = @(der) isscalar(der) && ... 
                  isPositiveIntegerValuedNumeric(der, true) && ...
                  ( der <= MAX_DER ) ;

if nargin < 3, scale = 1; end

if (shape <= 0 || scale <= 0) %if (shape <= 0 | scale <= 0)
  error('dinvgamma:negativeParameter', ...
  'Shape or scale parameter negative.')
end

if nargin < 4, derOrder = 0; end
if ~isValidDer(derOrder)
  error('dinvgamma:invalidDerivative', ...
  ['Order of derivative of the log density (inverse gamma distribution) ', ...
  'must be an integer between 0 and %d.'], MAX_DER);
end

if nargin < 5, isLog = false; end
if ~islogical(isLog)
  error('dinvgamma:incorrectType',...
 'input isLog must be a logical.'); 
end
% if ~isscalar(isLog)
  % error('dinvgamma:incorrectInputDimension',...
 % 'Input isLog must be a scalar.'); 
% end 
if nargin < 6, isLogX = false; end
if ~islogical(isLogX)
  error('dinvgamma:incorrectType',...
 'input isLogX must be a logical.'); 
end
% if ~isscalar(isLogX)
  % error('dinvgamma:incorrectInputDimension',...
  % 'Input isLogX must be a scalar.'); 
% end 
%% Define and vectorize the parameters
shape = zeros(1, length(x)) + shape;
scale = zeros(1, length(x)) + scale;
%% Compute the derivative of the log density of the inv. gamma distro 
if derOrder == 0
  if isLogX
  %dinvgamma1_logx: log density of the inverse gamma distribution, as a function of log(x)
    derLogDensity = arrayfun(@(X, A, B) ...
                    dinvgamma1_logx(X, A, B), x, shape, scale);
  else
  %dinvgamma1: log density of the inverse gamma distribution
    derLogDensity = arrayfun(@(X, A, B) ...
                    dinvgamma1(X, A, B), x, shape, scale);
  end

elseif derOrder == 1
  if isLogX
  % ddinvgamma1_logx: first derivative of the log density of the inverse gamma distribution, 
  % as function of log(x)
    derLogDensity = arrayfun(@(X, A, B) ...
                    ddinvgamma1_logx(X, A, B), x, shape, scale);    
  else
  % ddinvgamma1: first derivative of the log density of the inverse gamma distribution
    derLogDensity = arrayfun(@(X, A, B) ...
                    ddinvgamma1(X, A, B), x, shape, scale);
  end
  
elseif derOrder == 2  
  if isLogX
  % d2dinvgamma1_logx: second derivative of the log density of the inverse gamma distribution, 
  % as a function of log(x)
    derLogDensity = arrayfun(@(X, A, B) ...
                    d2dinvgamma1_logx(X, A, B), x, shape, scale);     
  
  else  
  % d2dinvgamma1: second derivative of the log density of the inverse gamma distribution
    derLogDensity = arrayfun(@(X, A, B) ...
                    d2dinvgamma1(X, A, B), x, shape, scale);    
  end
  
end 

%% Return the output in desired form  
if isLog
  derLogDensityIgam = derLogDensity;
else
  derLogDensityIgam = exp(derLogDensity);
end
%
% For comparison, snippets of the original function definitions in C++ code are 
% given in comments on top of each nested function defined in this MATLAB code.
%% Nested functions equivalent of those found in: 
% BayesFactor_0.9.12-4.2/BayesFactor/src/dinvgamma.cpp [1].
  % double dinvgamma1_Rcpp(const double x, const double a, const double b){
  % return a * log( b ) - lgamma( a ) - ( a + 1 ) * log( x ) - b / x ;
  % }
  function [derLogDensity] = dinvgamma1(x, a, b)
    derLogDensity = a .* log( b ) - real( lngamma( a ) ) - ( a + 1 ) .* log( x ) - b./x; 
  end %end-of-nested function: dinvgamma1

  % double dinvgamma1_logx_Rcpp(const double x, const double a, const double b){
  % return a * log( b ) - lgamma( a ) - ( a + 1 ) * x - b * exp( -x ) ;
  % }
  function [derLogDensity] = dinvgamma1_logx(x, a, b)
    derLogDensity = a .* log( b ) - real( lngamma( a ) ) - ( a + 1 ) .* x - b .* exp( -x ) ;
  end %end-of-nested function: dinvgamma1_logx

  % double ddinvgamma1_Rcpp(const double x, const double a, const double b){
  % return -( a + 1 ) / x + b / ( x * x ) ;
  % }
  function [derLogDensity] = ddinvgamma1(x, a, b)
    derLogDensity = -( a + 1 ) ./ x + b ./ ( x .* x ) ;
  end %end-of-nested function: ddinvgamma1
  
  function [derLogDensity] = ddinvgamma1_logx(x, a, b)
    % derLogDensity = -( a + 1 ) .* exp( -x ) + b .* exp( -2.*x );
    derLogDensity = -( a + 1 ) + b .* exp( -x );
  end %end-of-nested function: ddinvgamma1_logx
  
  % double d2dinvgamma1_Rcpp(const double x, const double a, const double b){
  % return ( a + 1 ) / ( x * x ) - 2 * b / (x * x * x) ;
  % }
  function [derLogDensity] = d2dinvgamma1(x, a, b)
    derLogDensity = ( a + 1 ) ./ ( x .* x ) - 2 .* b ./ (x .* x .* x) ;
  end %end-of-nested function: d2dinvgamma1
  
  function [derLogDensity] = d2dinvgamma1_logx(x, a, b)
  %%% WRONG!
    derLogDensity = ( a + 1 ) .* exp( -2.*x ) - 2 .* b .* exp( -3.*x ) ;
  end %end-of-nested function: d2dinvgamma1_logx
  
end %end-of-dinvgamma() 
%
%% The original source-code of the "dinvgamma" function written in R for BayesFactor package, 
% found under /BayesFactor_0.9.12-4.2/BayesFactor/R/common.R (starting at line 270).
%
% dinvgamma = function (x, shape, scale = 1, log = FALSE, logx = FALSE)
% {
    % if (shape <= 0 | scale <= 0) {
        % stop("Shape or scale parameter negative in dinvgamma().\n")
    % }
    % shape = rep(0, length(x)) + shape
    % scale = rep(0, length(x)) + scale
    % if(logx){
      % log.density = mapply(dinvgamma1_logx_Rcpp, x = x, a = shape, b = scale)
    % }else{
      % log.density = mapply(dinvgamma1_Rcpp, x = x, a = shape, b = scale)
    % }
    % if(log){
      % return(log.density)
    % }else{
      % return(exp(log.density))
    % }
% }
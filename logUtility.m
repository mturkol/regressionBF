function [logUtilOut] = logUtility(funName, funIn)
% Name       : logUtility.m
% Version    : 1.0
% Author     : Mert Turkol
% Date       : 05/08/2019
% Copyright  : GNU General Public License v2.0 (GPLv2)
% Description: This function computes accurately various logarithm/exponential 
%              related expressions which may suffer from loss of precision when 
%              the most direct method of evalution could cause issues like 
%              under/over-flow. 
% Expressions: log(1 + exp(x))
%              log(exp(x) + exp(y))
%              exp(x) - 1                         
% TO DO      : Add  logExpXminusExpY, Logit(p), LogitInverse(x), 
%                   LogLogitInverse(x), LogitInverseDifference(x, y), 
%                   ComplementaryLogLog(p), ComplementaryLogLogInverse(x)
% Syntax : 
%                  logUtility() :
%                     Without input arguments, display the help documentation
%
%   [logUtilOut] = logUtility(funName, funIn) :
%                     Compute and return the expression for 'funName', 
%                     using its required input(s) 'funIn'.  
%
%   [logUtilOut] = logUtility(funName, 'query') :
%                     Print the number of required input arguments to evaluate 
%                     the expression for 'funName'. 
% ==============================================================================
%  INPUTS (required, in order): 
%    funName : name of the function to evaluate the expression of, (char-array)
%              'log1pexp' -> log( 1 + exp(x) )
%              'logexpxplusexpy' -> log( exp(x) + exp(y) ) 
%              'expminusone' -> exp(x) - 1
%    funIn   : array of input(s) required to compute 'funName', (vector double)
%              OR
%              'query' -> print #input required to compute 'funName'
%
%  OUTPUT :
%    logUtilOut: Result of the evaluated expression 'funName', (scalar double)
%                OR
%                Character output for #inputs required to evaluate 'funName'
% ==============================================================================
% Comments: The code developed in 'logUtility.m' borrows ideas from: 
%   [1] John Cook (2014). Avoiding Overflow, Underflow, and Loss of Precision.
%       http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
%   [2] Morey, R. D., Rouder, J. N., Jamil, T., & Morey, M. R. D. (2015). 
%       Package 'bayesfactor'. 
%       http://cran/r-projectorg/web/packages/BayesFactor/BayesFactor.pdf
% ============================================================================== 
if nargin < 1, help logUtility; return; end
%% Initial input-checking
narginchk(1,2)

isCharVec = @(arg) ischar(arg) && isvector(arg) && size(arg, 1) == 1;
% isFiniteRealDblVec =  @(arg) isvector(arg) && isa(arg, 'double') && ...
                             % all(isfinite(arg)) && isreal(arg);
isDblVec = @(arg) isa(arg, 'double') && isvector(arg) && all( ~isnan(arg) );
%isSymVec = @(arg) isvector(arg) && isa(arg, 'sym') && all( ~isnan(arg) );
if ~isCharVec(funName) 
  error('logUtility:IncorrectType',...
       ['Input ''funName'' must be a character array, not a %s.'], ...
       class(funName) ); 
end
if nargin < 2, funIn = 'query'; end 

%if isSymVec(funIn), condSym = isreal(funIn); end  
condQry = strcmpi(funIn,'query'); 
condDbl = isDblVec(funIn);
if ( ~condQry && ~condDbl) 
  error('logUtility:IncorrectType',...
       ['Input ''funIn'' must be ', ...
        'either the char. array ''query'' or a vector of double(s).'] ); 
end
%% Define the function-handle and the required #input for the desired funName  
switch lower(funName)
  % return log(1 + exp(x)), preventing cancellation and overflow [1]
  case 'log1pexp'
    funName = 'log1pExp';
    numInReq = 1;
    logUtilH = @log1pExp;
    %logUtilOut = logUtilH(x);
    
  % Compute log(exp(x) + exp(y))    
  case 'logexpxplusexpy'
    funName = 'logExpXplusExpY';
    numInReq = 2;
    logUtilH = @logExpXplusExpY;
    %logUtilOut = logUtilH(x, y);
    
  % Calculate exp(x) - 1.
  % The most direct method is inaccurate for very small arguments.    
  case 'expminusone'
    funName = 'ExpMinusOne';
    numInReq = 1;
    logUtilH = @ExpMinusOne;
    %logUtilOut = logUtilH(x);
    
  otherwise
    error('logUtility:InexistentFunction',...
    '''%s'' does not exist as a utility function to call within ''logUtility()''.', ...
    funName);
end
%% Error-checking on the existence of function-name and proper #inputs
if condQry %when the user asks the #input arguments required for funName
  fprintf('''%s'' utility function within ''logUtility()'' requires %s input(s).\n', ...
          funName, num2str(numInReq) );
  return
end

numIn = numel(funIn);
if (numIn ~= numInReq) 
  error('logUtility:IncorrectArgumentCount',...
       ['''%s'' utility function within ''logUtility()'' requires %s input(s).\n', ...
        'User provided #input-argument(s): %s '], ...
        funName, num2str(numInReq), num2str(numIn) );
end
%   isAllScalar = all( arrayfun( @(idx) isscalar( funIn(idx) ), 1:numIn ) );
%   if ~isAllScalar
%     error('');
%   end
%% Call the desired nested-function via fun. handle and return its output
logUtilOut = logUtilH(funIn);
  
  %% Nested Function: log1pExp  
  function [funOut] = log1pExp(funIn)
    x = funIn;
    LOG_DBL_EPSILON = log( eps('double') );
    LOG_ONE_QUARTER = log(0.25);

    if (x > -LOG_DBL_EPSILON)
      % log(exp(x) + 1) == x to machine precision
      funOut = x;
    elseif (x > LOG_ONE_QUARTER)
      funOut = log( 1.0 + exp(x) );
    else
      % Prevent loss of precision that would result from adding small argument to 1.
      funOut = log1p( exp(x) );
    end
  end %end-of-nested fun: log1pExp
  %% Nested Function: logExpXplusExpY 
  function [funOut] = logExpXplusExpY( funIn )
    x = funIn(1); y = funIn(2);
    funOut = x + log1pExp( y - x );
  end %end-of-nested fun: logExpXplusExpY
  %% Nested Function: ExpMinusOne
  function [funOut] = ExpMinusOne(funIn)
    x = funIn;
    w = exp(x);
    
    % The following expression is algebraically equal to exp(x) - 1.
    % The advantage in finite precision arithmetic is that
    % it avoids subtracting nearly equal numbers.       
    rexp = w * ( 0.5 + ( 0.5 - 1.0/w ));
    
    % Use rational approximation for small arguments.
    if( abs(x) < 0.15 )
      p1 =  0.914041914819518e-09;
      p2 =  0.238082361044469e-01;
      q1 = -0.499999999085958e+00;
      q2 =  0.107141568980644e+00;
      q3 = -0.119041179760821e-01;
      q4 =  0.595130811860248e-03;      
      
      rexp = x*(((p2*x + p1)*x + 1.0)/((((q4*x + q3)*x + q2)*x + q1)*x + 1.0));

    % For large negative arguments, direct calculation is OK. 
    elseif ( x <= -0.15 )
      
      rexp = w - 1.0;
      
    end  
   
    funOut = rexp;

  end %end-of-nested fun: ExpMinusOne
  
end %end-of-logUtility
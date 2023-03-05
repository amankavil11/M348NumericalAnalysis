%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NEWTON'S METHOD
%   
%  Solves the problem
%    f(x) = 0
%  using Newton's method. For a known true solution calculates errors.
%   
%  To run, type 'newton_err' (without quotes) on the command line.
%  The main function is newton:
% 
%  [state, x, errors, iter] = newton(@f, @df, x0, tolerance, maxIteration, debug);
%   
%  Inputs:
%    @f            Handle to function f
%    @df           Handle to the derivative of function f
%    x0            The initial guess at the solution.
%    tolerance     The convergence tolerance (must be > 0).
%    maxIteration  The maximum number of iterations that can be taken.
%    debug         Boolean for printing out information on every iteration.
%  Outputs:
%    x             The solution
%    errors        Array with errors at each iteration
%    iter          number of iterations to convergence
%  Return:

%    state         An error status code.
%      SUCCESS     Sucessful termination.
%      WONT_STOP   Error: Exceeded maximum number of iterations.
%      BAD_ITERATE Error: The function had a vanishing derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

global SUCCESS WONT_STOP BAD_ITERATE
SUCCESS = 0; 
WONT_STOP = 1;
BAD_ITERATE = 2;

global x_true
% func  = @(x) x*x-2;
% dfunc = @(x) 2*x;
% x_true = sqrt(2);
dfunc  = @(x) ((324*x^5)+(225*x^4)-(408*x^3)-(207*x^2)+(70*x)+(16));
func = @(x) ((54*x^6)+(45*x^5)-(102*x^4)-(69*x^3)+(35*x^2)+(16*x)-4);
x_true = -0.666667;

% Input
diary AP3IterHistory.txt

disp('Solves the problem f(x) = 0 using Newton''s method')
x0 = input('Enter guess at root: ');
tol = input('Enter tolerance: ');
maxIter = input('Enter maxIteration: ');
debug = input('Monitor iterations? (1/0): ');

% Solve for Newton method

[s, x, errors, iter] = newton(func, dfunc, x0, tol, maxIter, debug);

% Report results
format longe
switch s
    case SUCCESS
        fprintf('The root is %.6f\n', x)
        fprintf('The number of iterations is %d\n', iter)
    case WONT_STOP
        fprintf('ERROR: Failed to converge in %d iterations!\n', maxIter);
    case BAD_ITERATE
        fprintf('ERROR: Obtained a vanishing derivative!\n');
        return
    otherwise
        disp('ERROR: Coding error!')
        return
end

diary off

errors = errors(1:iter+1);
disp('Press a button to continue ...')
pause
errors
format short


function [state, x, errors, iter] = newton(f, df, x0, tolerance, maxIteration, debug)

global SUCCESS WONT_STOP BAD_ITERATE x_true

% format string
   prec = 12; % number of significant digits to output
   fmt = sprintf('Iter %%d: x = %%.%dg, error = %%.%dg, error_ratio_quad = %%.%dg, error_ratio_linear = %%.%dg\n',prec,prec,prec,prec);
   eps = 1e-20;   
   errors = zeros(maxIteration+1,1);
   x = x0;
   err = abs(x - x_true);
   errors(1) = err;
   
   if (debug)
       fprintf('Guess: x = %.8g, error = %.8g\n', x, err);
   end
      
  % Newton loop
  for itn = 1:maxIteration
      
      dfx = df(x);
      if (abs(dfx) < eps)
          state = BAD_ITERATE;
          iter = itn; 
          return  
      end
      dx = -f(x)/dfx;
      dx = 2*dx; % Modified Newton, multiplicity 2
      x = x+dx;
      err = abs(x - x_true);
      errors(itn+1) = err; 
      
      if (debug)
%           fprintf('Iter %d: x = %f, dx = %.8g, error = %.8g\n', itn, x, dx, err);
          err_quad = err/((errors(itn))^2);
          err_lin = err/(errors(itn));
          fprintf(fmt, itn, x,err, err_quad, err_lin)
      end
      
      % Check error tolerance
      if (abs(dx) <= tolerance)
          iter = itn; 
          state = SUCCESS;
          return
      end
      
  end
  state = WONT_STOP;
  iter = itn; 
  return 
end




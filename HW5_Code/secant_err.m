%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SECANT METHOD
%   
%  Solves the problem
%    f(x) = 0
%  using Secant method. For a known true solution calculates errors.
%   
%  To run, type 'secant_err' (without quotes) on the command line.
%  The main function is secant:
% 
%  [state, x, errors] = secant(@f, x0, x1, tolerance, maxIteration, debug);
%   
%  Inputs:
%    @f            Handle to function f
%    x0            The initial guess at the solution.
%    x1            The second guess at the solution.
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
func  = @(x) (54*x^6+45*x^5-102*x^4-69*x^3+35*x^2+16*x-4); 
x_true = -1.381298;

% Input

disp('Solves the problem f(x) = 0 using Secant method')
x0 = input('Enter guess 0 at root: ');
x1 = input('Enter guess 1 at root: ');
tol = input('Enter tolerance: ');
maxIter = input('Enter maxIteration: ');
debug = input('Monitor iterations? (1/0): ');

% Solve for Secant method

[s, x, errors, iter] = secant(func, x0, x1, tol, maxIter, debug);

% Report results
format longe
switch s
    case SUCCESS
        fprintf('The root is %.6g\n', x)
        fprintf('The number of iterations is %d\n', iter)
    case WONT_STOP
        fprintf('ERROR: Failed to converge in %d iterations!\n', maxIter);
    case BAD_ITERATE
        fprintf('ERROR: Obtained a vanishing derivative\!\n');
        return
    otherwise
        disp('ERROR: Coding error!')
        return
end


format short


function [state, x, errors, iter] = secant(f, x0, x1, tolerance, maxIteration, debug)

global SUCCESS WONT_STOP BAD_ITERATE x_true

% format string
   prec = 6; % number of significant digits to output
   fmt = sprintf('Iter %%d: x = %%.%dg, dx = %%.%dg, error = %%.%dg, superlinear error ratio = %%.%dg, linear error ratio = %%.%dg\\n',prec,prec,prec,prec,prec);
   eps = 1e-20; 

   errors = zeros(maxIteration+2,1);
   x = x1; fx0 = f(x0);
   errors(1) = abs(x0 - x_true);
   errors(2) = abs(x1 - x_true);
   
   if (debug)
       fprintf('Guess 0: x = %f, error = %.6g\n', x0, errors(1));
       fprintf('Guess 1: x = %f, error = %.6g\n', x1, errors(2));       
   end
      
  % Secant loop
  for itn = 1:maxIteration
      
      fx = f(x);
      if (abs(fx-fx0)<eps)
          state = BAD_ITERATE;
          iter = itn; 
          return  
      end
      dx = -f(x)*(x-x0)/(fx-fx0);
      x0 = x; fx0 = fx;
      x = x+dx;
      err = abs(x - x_true);
      errors(itn+2) = err; 
      
      if (debug)
          if itn > 2
            fprintf(fmt, itn, x, dx, err, (errors(itn)/(errors(itn-1)*errors(itn-2))), errors(itn+2)/errors(itn+1) );
          else
              fprintf(fmt, itn, x, dx, err, 0, errors(itn+2)/errors(itn+1) );
          end
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BISECTION METHOD
% 
%  Solves the problem
%    f(x) = 0
%  using the bisection algorithm.
%  
%  To run, type 'bisection' (without quotes) on the command line.
%  The main function is bisect:
% 
%  [state, x] = bisect(@f, a, b, tolerance, maxIteration, debug);
% 
%  Inputs:
%    @f            Handle to function f
%    a,b           The initial bounding interval, with a root between.
%    tolerance     The convergence tolerance (must be > 0).
%    maxIteration  The maximum number of iterations that can be taken.
%    debug         Boolean for printing out information on every iteration.
%  Outputs:
%    x             The solution.
%  Return:
%    state         An error status code.
%      SUCCESS     Sucessful termination.
%      WONT_STOP   Error: Exceeded maximum number of iterations.
%      BAD_DATA    Error: The interval may not bracket a root.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

global SUCCESS WONT_STOP BAD_DATA
SUCCESS = 0; 
WONT_STOP = 1;
BAD_DATA = 2;

diary iterHistory.txt

func = @(x) (cos(x) - sin(x));    % Example 1.2 in Sauer

% Input

disp('Solves the problem f(x) = 0 on interval [a,b] using the bisection algorithm')
a = input('Enter a: ');
b = input('Enter b: ');
tol = input('Enter tolerance: ');
maxIter = input('Enter maxIteration: ');
debug = input('Monitor iterations? (1/0): ');

% Solve for a root

[s, x] = bisect(func, a, b, tol, maxIter, debug);

% Report results

switch s
    case SUCCESS
        fprintf('The root is %.6g\n', x)
        fprintf('f(%.6g) = %.6g \n', x, func(x))
    case WONT_STOP
        fprintf('ERROR: Failed to converge in %d iterations!\n', maxIter);
    case BAD_DATA
        disp('ERROR: Unsuitable interval!');
    otherwise
        disp('ERROR: Coding error!')
end

diary off

function [err, x] = bisect(f, a, b, tolerance, maxIteration, debug)

global SUCCESS WONT_STOP BAD_DATA

% format string
prec = 6; % number of significant digits to output
fmt = sprintf('Iter %%d: x = %%.%dg, dx = %%.%dg, a = %%.%dg, b = %%.%dg, f(x) = %%.%dg\\n', ...
    prec,prec,prec,prec,prec+4);
  % Swap a and b if necessary so a < b
  if (a > b) 
      c = a; a = b; b = c; 
  end
  fa = f(a); fb = f(b);

  % Make sure there is a root between a and b
  if (sign(fa)*sign(fb) > 0.0) 
      err = BAD_DATA; 
      return
  end

  % Bisection iteration loop
  dx = b-a;
  for iter = 0:maxIteration
      
      dx = dx/2;
      x = a + dx;
      fx = f(x);
      
      if (debug)
%           fprintf('Iter %d: x = %f, dx = %f, a = %f, b = %f, f(x) = %f\n',  ...
%               iter, x, dx, a, b, fx);
          fprintf(fmt, iter, x, dx, a, b, fx);
      end
      
      % Check error tolerance
      if (dx <= tolerance)
          err = SUCCESS;
          return
      end
      
      if (sign(fa)*sign(fx) > 0.0)
          a = x; fa = fx;
      else
          b = x; fb = fx;
      end
      
  end
  err = WONT_STOP;
  return 
end



  
  


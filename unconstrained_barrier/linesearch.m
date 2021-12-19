function [xkp1,fxkp1,niter,tk] = linesearch(fun,fxk,gradfxk,xk,pk,tkmax,beta,c,nitermax)
% LINESEARCH Computes the value of decision variables xkp1 at iteration k+1
% by carrying out a back-tracking line search approach, trying to achieve a
% sufficient decrease of the cost function f(xk+tk*pk) where pk is the
% search direction. pk must be locally a descent direction (i.e. gradfxk'*pk<0).
%
%   INPUTS:
%           fun         =   cost function
%           fxk         =   cost function evaluated at xk
%           gradfxk     =   gradient of the cost function with respect to
%                           x, evaluated at xk
%           xk          =   value of the decision variables at the current
%                           iterate
%           pk          =   search direction
%           tkmax       =   maximum step size
%           beta        =   ratio tk_ip1/tk_i during the back-tracking line
%                           search (beta in (0,1))
%           c           =   ratio between acheived decrease and decrease
%                           predicted by the first-order Taylor expansion,
%                           sufficient to exit the back-tracking line
%                           search algorithm
%           nitermax    =   maximum number of iterations
%
%   OUTPUTS:
%           xkp1        =   obtained value of xk+tk*pk
%           fxkp1       =   cost function evaluated at xkp1
%           tk          =   step length
%           niter       =   number of iterations employed to satisfy Armijo
%                           condition

niter       =       0;
tk          =       tkmax;
xkp1        =       xk+tk*pk;
fxkp1       =       fun(xkp1);

while   isnan(fxkp1) || (fxkp1>fxk+tk*c*gradfxk'*pk && niter<nitermax)
    tk      =       beta*tk;
    xkp1    =       xk+tk*pk;
    fxkp1   =       fun(xkp1);
    niter   =       niter+1;
end
    
    
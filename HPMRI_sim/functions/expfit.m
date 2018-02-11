function [estimates, fitted, rss, rmse] = expfit(xdata, ydata, start_point)

% EXPFIT Fits an exponential decay curve:  y = A*exp(-lambda*x)
%
% Syntax:
%   [ESTIMATES, FITTED, RSS, RMSE] = EXPFIT(XDATA, YDATA)
%
% Inputs:
%   XDATA = vector of independent data
%   YDATA = vector of dependent data
%
% Outputs:
%   ESTIMATES = estimated parameters [A, lambda]
%   FITTED = vector of fitted dependent data
%   RSS = sum of the squared residuals.
%   RMSE = root mean squared error.
%
% See also:  expfit2
%
% Based on a MATLAB example.
% 10/26/10, Dave J. Niles, dave.j.niles@gmail.com

% Check for a starting point.
if ~exist('start_point','var'),
    start_point = rand(1, 2);
end

% Call fminsearch.
model = @expfun;
options = optimset('MaxIter',1e15,'MaxFunEvals',1e15);
estimates = fminsearch(model, start_point,options);
fitted = estimates(1).*exp(-estimates(2).*xdata);
rss = sum((fitted-ydata).^2);
rmse = sqrt(rss)./length(fitted);

% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [Error, FittedCurve] = expfun(params)
        A = params(1);
        lambda = params(2);
        FittedCurve = A .* exp(-lambda * xdata);
        ErrorVector = FittedCurve - ydata;
        Error = sum(ErrorVector .^ 2);
    end
end

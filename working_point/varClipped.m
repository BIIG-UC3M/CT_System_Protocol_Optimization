function [v E] = varClipped(gX,maxB,h,n)

% 
% function v = varClipped(gX,maxB,h,n)
%
% Function to estimate the variance assuming a gaussian clipped model,
% lets test its validity and use it to simulate the dynamic range by
% not use simulations, just analytical expressions...
% Model is accurately described @ "Econometric Analysis" Greene, 2002 and
% it is referred as the Gaussian Censored model (chapter 22)
% 
%

% Temporary terms

sigma_sq = h*gX + n;              % Variance of the original (non-clipped) signal
sigma    = sqrt(sigma_sq);        % Square root of variance
alpha    = (maxB-gX)/sigma;       % Alpha parameter
pdf      = normpdf(alpha);        % Prob for the given evaluation point
cdf      = normcdf(alpha);        % Accumulated prob for this point
lambda   = -(pdf/cdf);            % Auxiliary parameter
delta    = lambda*(lambda-alpha); % Auxiliary parameter

% Expression for the variance to return
v = sigma_sq*cdf*((1-delta)+((alpha-lambda)^2)*(1-cdf));

% Now the mean
E = (1-cdf)*maxB + cdf*(gX+sigma*lambda);

end
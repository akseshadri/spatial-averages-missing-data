function [Fdist,xdist] = getcumdist(x,N)
% GETCUMDIST  Empirical cumulative distribution function.
%
%   [Fdist, xdist] = getcumdist(x, N) computes the empirical CDF of the
%   data vector x, evaluated at N evenly spaced points spanning slightly
%   beyond the range of x.
%
%   Inputs:
%     x  — Data vector (e.g. optimal weights beta_i).
%     N  — Number of evaluation points for the CDF.
%
%   Outputs:
%     xdist — 1 x N vector of evaluation points.
%     Fdist — 1 x N vector of CDF values: Fdist(i) = fraction of elements
%             in x that are <= xdist(i).
%
%   Used in PFig4.m to plot cumulative distributions of optimal weights.

range = max(x) - min(x);
delta = 1e-5;
xdist = linspace(min(x)-delta*range,max(x)+delta*range,N);

Fdist = zeros(1,N);

for i = 1:N
    xi = xdist(i);
    Fdist(i) = sum(x<=xi)/numel(x);
end

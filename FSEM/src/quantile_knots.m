%	Create knots at sample quantiles.  If boundstab == 1, then nknots + 2
%	knots are created and the first and last are deleted.  This
%	mitigates the extra variability of regression spline estimates near
%	the boundaries.
%	
%		INPUT (required)
%	x = independent variable.  (The knots are at sample quantiles of x.)
%	nknots = number of knots
%
%		INPUT (optional)
%	boundstab = parameter for boundary stability (DEFAULT is 0)
%
%	USAGE: [knots] = quantile_knots(x, nknots, boundstab) ;
%
% 09/07/2014
% Shikai Luo

function [knots] = quantile_knots(x, nknots, boundstab)

if nargin < 3 
	boundstab = 0;
end

n = length(x);
xsort = sort(x);   

location = n*(1 : (nknots + 2*boundstab))' / (nknots + 1 + 2*boundstab);
knots=xsort(round(location));
knots=knots((1 + boundstab) : (nknots + boundstab));

end
%  bsplineMat  Computes values or derivative values of B-spline basis functions
%  Arguments:
%  x         ... Argument values for which function values are computed
%  breaks    ... Increasing knot sequence spanning argument range
%  norder    ... Order of B-spline (one greater than degree) max = 19
%                Default 4.
%  nderiv    ... Order of derivative required, default 0.
%  sparsewrd ... if 1, return in sparse form
%  Return:
%  bsplineM  ... length(x) times number of basis functions matrix
%                 of bspline values

% 09/22/2014
% Shikai Luo
% knots = [min(x); quantileKnots(x, nknots); max(x)];


function [bsplineM] = bsplineMat(x, breaks, norder, nderiv, sparsewrd)

sizex = size(x);
if sizex(1) > 1 && sizex(2) > 1
    error('Argument X is not a vector.');
end
x = x(:);

n = length(x);

if nargin < 5, sparsewrd = 1;  end
if nargin < 4, nderiv    = 0;  end
if nargin < 3, norder    = 4;  end

if nargin < 2
    error('BREAKS not supplied as the second argument.');
end

nbreaks = length(breaks);
if nbreaks < 2
    error('Number of knots less than 2.');
end

if norder < 1
    error('Order of basis less than one.');
end

if nderiv < 0
    error('NDERIV is negative');
end

if nderiv >= norder
    error('NDERIV cannot be as large as order of B-spline.');
end

if min(diff(x)) < 0
    [x,isrt] = sort(x);
    sortwrd = 1;
else
    sortwrd = 0;
end

if x(1) - breaks(1) < -1e-10 || x(n) - breaks(nbreaks) > 1e-10
    disp([x(1), x(n)])
    error ('Argument values out of range.')
end

k   = norder;        
km1 = k-1;
nb  = length(breaks);
nx  = length(x);      
nd  = nderiv+1;       
ns  = nb - 2 + k;  
if ns < 1
    fprintf('There are no B-splines for the given input.\n')
    bsplineM = [];
    return
end
onenx = ones(nx,1);
onenb = ones(k, 1);
onens = ones(ns,1);

if size(breaks,1) > 1
    breaks = breaks';
end
knots  = [breaks(1)*ones(1,km1), breaks, breaks(nb)*ones(1,km1)]';
nbasis = length(knots) - k;

knotslower      = knots(1:nbasis);
[~,index] = sort([knotslower', x']);
pointer         = find(index > nbasis) - (1:length(x));
left            = max([pointer; k*onenx']);

temp = [1, zeros(1,km1)];
b    = temp(ones(nd*nx,1),:);
nxs  = nd*(1:nx);

for j=1:k-nd
    saved = zeros(nx,1);
    for r=1:j
        leftpr   = left + r;
        tr       = knots(leftpr) - x;
        tl       = x - knots(leftpr-j);
        term     = b(nxs,r)./(tr+tl);
        b(nxs,r) = saved + tr.*term;
        saved    = tl.*term;
    end
    b(nxs,j+1)  = saved;
end

for jj=1:nd-1
    j = k - nd + jj;
    saved = zeros(nx,1);
    nxn   = nxs - 1;
    for r=1:j
        leftpr   = left + r;
        tr       = knots(leftpr) - x;
        tl       = x - knots(leftpr-j);
        term     = b(nxs,r)./(tr+tl);
        b(nxn,r) = saved + tr.*term;
        saved    = tl.*term;
    end
    b(nxn,j+1)  = saved;
    nxs = nxn;
end

for jj=nd-1:-1:1
    j = k - jj;
    temp = (jj:nd-1).'*onenx' + ones(nd-jj,1)*nxn;
    nxs = reshape(temp,(nd-1-jj+1)*nx,1);
    for r=j:-1:1
        leftpr     = left + r;
        temp       = ones(nd-jj,1)*(knots(leftpr) - knots(leftpr-j)).'/j;
        b(nxs,r)   = -b(nxs,r)./temp(:);
        b(nxs,r+1) =  b(nxs,r+1) - b(nxs,r);
    end
end

index = find(x < breaks(1) | x > breaks(nb));
if ~isempty(index)
    temp = (1-nd:0).'*ones(1,length(index))+nd*ones(nd,1)*index(:).';
    b(temp(:),:) = zeros(nd*length(index),k);
end

width = max([ns,nbasis]) + km1 + km1;
cc    = zeros(nx*width,1);
index = (1-nx:0).'*onenb' + nx*(left.'*onenb' + onenx*(-km1:0));
cc(index) = b(nd*(1:nx),:);

bsplineM = reshape(cc((1-nx:0).'*onens' + nx*onenx*((1:ns))), nx, ns);

if sortwrd
    temp = bsplineM;
    bsplineM(isrt,:) = temp;
end

if sparsewrd
    bsplineM = sparse(bsplineM);
end

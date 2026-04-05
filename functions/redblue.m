function cmap = redblue(n)
if nargin < 1
    n = 256;
end
half = floor(n/2);

r = [linspace(0,1,half) ones(1,n-half)];
g = [linspace(0,1,half) linspace(1,0,n-half)];
b = [ones(1,half) linspace(1,0,n-half)];

cmap = [r(:) g(:) b(:)];
end
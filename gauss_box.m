% pcex 1.3.1
% Gaussian random number generator using Box-Muller

% Inputs   n     desired number of of random numbers (MUST BE EVEN)
%          s     seed for the uniform random variable (MUST BE NEGATIVE)
%          mu    (optional) mean of the distribution DEFAULT = 0
%          sig   (optional) standard deviation DEFAULT = 1

% Output   x     n-array storing n random values from the desired distribution
% =================================================================================
function x = gauss_box(s,n,varargin)
    reshape = 0;
    if mod(n,2) ~= 0
        n = n + 1;
        reshape = 1;
    end
    
    v = ran2(s,n);
    n2 = n/2;
    g = zeros(n,1);
    vlog = sqrt(-2 * log(v(1:n2)));
    vpi = 2 * pi * v(n2+1:n);
    g(1:n2) = vlog .* cos(vpi);
    g(n2+1:n) = vlog .* sin(vpi);
    
    if nargin == 4
        g = g * varargin{2} + varargin{1};
    end
    
    if reshape
        x = g(1:end-1);
    else
        x = g;
    end
end
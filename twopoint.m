% pcex 1.3.1
% Two-point random number generator

% Inputs   n     desired number ofrandom numbers
%          s     seed for the uniform random variable
%          x1, x2, p1    parameters for two-point RNG
%          p1    probability of having value x1

% Output   x     nx1 array storing n random values from the desired distribution
% =================================================================================
function t = twopoint(s,n,x1,x2,p1)    
    v = ran2(s,n);    
    t = zeros(n,1);
    t(v < p1) = x1;
    t(v >= p1) = x2;    
end
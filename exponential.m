% pcex 1.3.1
% Exponential random number generator using inverse transform X = -ln(1-U)/lambda, 
% where X is the exponential random value corresponding to the uniform random value U, 
% and lambda is the parameter of the desired exponential function.

% Inputs   n     desired number ofrandom numbers
%          s     seed for the uniform random variable
%          lambda    parameter for exponential distribution

% Output   x     nx1 array storing n random values from the desired distribution
% ==================================================================================
function x = exponential(s,n,lambda)
    v = ran2(s,n);    
    x = -log(v) / lambda; %instead of 1-v, here we use only v to save computation
end


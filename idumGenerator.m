% Generate pairs of uniform negative random numbers to be used as seeds 
% (idum) in other random number applications

% Inputs    seed1, seed2    seeds for the uniform random number generator
%                           on (0,1), must be negative
%           n               number of pairs of seeds desired, >= 1
%           lb              lower bound in absolute value of the interval (lb,0) from which
%                           the random numbers should be drawn from; must
%                           be positive

% Output    s               nx2 matrix storing pairs of seeds
% 27.02.2015
%==========================================================================
function s = idumGenerator(s,n,lb)    
    if (n < 1) || (s >= 0) || (lb <= 0)
        error('n positive; s negative; lb positive')
    end
    s = floor((ran2(s,n) / 2 - 1) * lb);
end
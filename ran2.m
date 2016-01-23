% Uniform random number generator. 
% Numerical Recipes in C, v.1999, 272.

% Long period (> 2 x 10^18) random number generator of L'Ecuyer with
% Bays-Durham shuffle and added safeguards.
% Returns uniform random deviate(s) between 0.0 and 1.0, exclusive.

% input   idum   initial seed, must be NEGATIVE
%         n      (optional) number of random values desired DEFAULT = 1

% output  num    array storing random value(s) from U(0,1)
% =========================================================================
function num = ran2(varargin)   
    idum = varargin{1};
    %----------------------------------------------------------------------
    % initialize static/persistent variables/constants
%     persistent idum2 iy iv IM1 IM2 AM IMM1 IA1 IA2 IQ1 IQ2 IR1 IR2 NTAB NDIV RNMX
    % This isempty check initializes the persistent variables if ran2() is 
    % called the first time; otherwise use persistent values from previous 
    % calls. Since idum is negative only for the very first call, here we
    % also the shuffle.
%     if isempty(idum2) 
    % initialize constants
    IM1 = 2147483563; 
    IM2 = 2147483399; 
    AM = 1.0/IM1;
    IMM1 = IM1-1; 
    IA1 = 40014;
    IA2 = 40692;
    IQ1 = 53668;
    IQ2 = 52774;
    IR1 = 12211;
    IR2 = 3791;
    NTAB = 32;
    NDIV = 1 + floor(IMM1/NTAB); % use floor() to preserve integer division
    RNMX = 1.0 - eps; % this uses matlab's eps
        iv = zeros(NTAB,1); 
        idum = max(-idum,1);
        idum2 = idum;
        for j = NTAB+8:-1:1
            k = floor(idum/IQ1);
            idum = IA1 * (idum - k*IQ1) - k*IR1; 
            idum = idum + IM1 * (idum < 0);
            if j <= NTAB
                iv(j) = idum;
            end
        end
        iy = iv(1); 
%     end
    %======================================================================

    % main program to compute a pseudo random number from the uniform
    % distribution
    if nargin == 1
        n = 1;
    else
        n = varargin{2};
    end
    
    num = zeros(n,1);
    
    for i = 1:n
    k = floor(idum/IQ1);  % floor() preserves integer division
    idum = IA1 * (idum - k*IQ1) - k*IR1;
    idum = idum + IM1 * (idum < 0);
    
    k = floor(idum2/IQ2); 
    idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
    idum2 = idum2 + IM2 * (idum2 < 0);

    j = floor(iy/NDIV) + 1; % + 1 to make index correspond to Matlab's 1:NTAB rule
    iy = iv(j) - idum2;
    iv(j) = idum;
    iy = iy + IMM1 * (iy < 1);    

    num(i) = min(AM * iy, RNMX);
    end
end
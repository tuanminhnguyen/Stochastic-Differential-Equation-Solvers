% Generate Wiener increments dW_t = W_t - W_{t-1}
% Recall Wiener process:
%             \Delta W_t ~ N(0,sqrt(\Delta))

% Inputs     t0,tf    [t0,tf] on which the Wiener process is generated
%            m         number of intervals
%            s         seed for the normal random number generator

% Output     dW         vector storing the increments dW_t for each
%                       end points according to interv and Delta (array length m)
%==========================================================================
function dW = WienerIncrement(t0,tf,m,s)
Delta = (tf - t0) / m;
dW = gauss_box(s,m,0,sqrt(Delta)); %mean 0, standard deviation sqrt(Delta)
end
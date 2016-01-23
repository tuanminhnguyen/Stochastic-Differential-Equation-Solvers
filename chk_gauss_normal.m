% check Gaussian normality of data, using histogram, mean and variance

% input    v    vector or matrix storing the data
% ouput    mu   mean
%          std  standard deviation
% =========================================================================
function [mu,var] = chk_gauss_normal(v)
    v = v(:);
    histogram(v,-6:10E-2:6);
    mu = sum(v) / (length(v) - 1);
    var = sum((v - mu).^2) / (length(v) - 1);
end


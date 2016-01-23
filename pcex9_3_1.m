% pcex9.3.1
% Simulate N trajectories of the Ito process X satisfying the conditions in
% chapter 9 (2.1) pg 308 and return error vector.

% 25.02.2015
%==========================================================================
global IFUNC
IFUNC = 1;
nintervals  = [16    32    64   128];
% delta = (2 * ones(4,1)) .^ ((-4:-1:-7)');
N=100;
K = length(nintervals);
err = zeros(K,1);  % error vector, each element correspond to one step size
seedsK = idumGenerator(-2,K,100);
X0 = 1.0;
a = 1.5; b = 0.1;
T = 1;
for i = 1:K
    t = linspace(0,1,nintervals(i)+1);        
    seedsN = idumGenerator(seedsK(i),N,100); 
    X = zeros(nintervals(i)+1,N);            % matrix to store the exact solution
    Y = zeros(nintervals(i)+1,N);            % approximate solution
    for j = 1:N
        dW = WienerIncrement(0,1,nintervals(i),seedsN(j));
        X(:,j) = exactItoSoln(X0,t,dW);
        Y(:,j) = eulerMaruyama(X0,nintervals(i),dW);        
    end   
    Sum = sum(abs(X(end,:) - Y(end,:)));
    err(i) = Sum/N;
end

delta = 1./nintervals;
fprintf('delta  |')
fprintf('%8.4f', delta)
fprintf('\n')
fprintf('-------------------------------------------\n')
fprintf('err    |')
fprintf('%8.4f', err)
fprintf('\n')
fprintf('N = %d\n', N)
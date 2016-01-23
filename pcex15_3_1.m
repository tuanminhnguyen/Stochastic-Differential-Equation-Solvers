global IFUNC
IFUNC = 3;
nnodes = [8 16 32 64];  %[128 256 512 1024];
M = 20; N = 100;
X0 = 0.1;
CI = zeros(2,4);
err = zeros(M,1);
for j = 1:4    
    X = zeros(nnodes(j) + 1,N);
    Y1 = zeros(nnodes(j) + 1,N);
    L = nnodes(j)/2;
    Y2 = zeros(L + 1,N);
    t = linspace(0,1,nnodes(j)+1);
    for i = 1:M
        for k = 1:N
            Delta = 1/nnodes(j);
            sqrtDelta = sqrt(1/nnodes(j));
            dW = sqrtDelta.*randn(nnodes(j),1) + 0;        
            X(:,k) = exactItoSoln(X0,t,dW);
            Ym = zeros(nnodes(j) + 1,1);
            Ym(1,1) = X0;
            for n = 1:nnodes(j)
                Ym(n+1,1) = (1 + 1.5 * Delta + 0.01*dW(n))*Ym(n);
            end    
            
            Y1(:,k) = Ym;            
            Delta = Delta * 2;            
            Ym = zeros(L+1,1);
            Ym(1) = X0;
            for n = 1:L
                Ym(n+1,1) = (1 + 1.5 * Delta + 0.01*dW(n))*Ym(n);
            end    
            Y2(:,k) = Ym;            
        end        
        Vi = 2*mean(Y1(end,:)) - mean(Y2(end,:));
        err(i) = Vi - mean(X(end,:));
    end
     [~,~,CI(:,j),~] = normfit(err,0.1);   
end

Delta = 1./nnodes;

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean error \epsilon against time step size |\Delta| in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex15_3_1.jpeg')

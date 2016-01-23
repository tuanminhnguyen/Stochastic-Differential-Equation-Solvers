% x = linspace(0,1,100);
% f = @(x,m) sin(2*m*pi*x);
% figure
% for m = 1:100
%     y = f(x,m);
%     plot(x,y);
%     pause(0.1)
% %     hold on
% end

f = @(x,m) (1+(-1)^(2*m))/((2*m)^2 + 4) * sin(
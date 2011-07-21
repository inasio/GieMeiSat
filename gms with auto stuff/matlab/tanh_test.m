% function tanh_test

% z1 = 1.23;
% z2 = sqrt(i);
% z3 = 3.4123 + i*0.123;
% 
% [tanh(z1), hyp_tanh(z1)]
% [tanh(z2), hyp_tanh(z2)]
% [tanh(z3), hyp_tanh(z3)]


h=figure(1)
clf()
for a=1:3
x = -5:0.0001:5;
y = 1./(sqrt(i*x).*(tanh(sqrt(i*x*a))));
plot(x,real(y)), hold on
end
text(-0.3,0.35,'a = 1','FontSize',15)
text(-0.3,0.49,'a = 2','FontSize',15)
text(-0.3,0.59,'a = 3','FontSize',15)
xlabel('x','fontsize',15)
ylabel('y(x)','fontsize',15)
title('r y(x) = \sqrt[ix]','fontsize',20);
print

% x = 0.001;
% L = 0:0.001:10;
% y = real(1./(sqrt(i)*x.*(tanh(sqrt(i)*x*L))));
% figure(2)
% plot(L,y)

% function y = hyp_tanh(z)
%   a=real(z);
%   b=imag(z);
%   
%   y = (sinh(a)*cos(b) + i*sin(b)*cosh(a))/(cosh(a)*cos(b) + i*sinh(a)*sin(b));
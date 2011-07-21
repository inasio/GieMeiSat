function gms_2d_radial

path_figs = '/home/rozada/Temporary/Figures/GMSpaper/';
name = 'GMS_splitting';
L_init = 0.8; L_end = 5.3;
N = 100;
T_range = [0,100];
Step_r = ceil(N/200);
L_all = linspace(L_init,L_end,200);
H = 1/(N+1);
R = (H:H:1-H)';
D = 1;
Eps = 0.02;
% Eps = 0.002;
Kap = 2.5;
Tau = 1;

R_l = 1/4;
U0 = 0.346*tanh(25*(R-R_l));
U0 = U0 + 0.346*tanh(25*(1-R_l-R));
V0 = 0.287 - 0.03*cos(2*pi*R);

[T0,Y0] = ode15s(@rhs, [0,100], [U0;V0],[],L_init);

Y = Y0(end,:)';
figure(1)
plot(R,Y0(1:N));

% tic
% for I = 2:length(L_all)
%   I/length(L_all)
%   [T,W] = ode15s(@rhs, T_range, Y(:,I-1),[],L_all(I));
% %   W(end,1:N)'
%   Y(:,I) = [W(end,1:N)';W(end,N+1:2*N)'];
% end
% toc
% 
% U = Y(1:N,:);
% V = Y(2*N:end,:);
% 
% figure(2)
% surf(R(3:Step_r:N),L_all,U(3:Step_r:N,:)')
% rotate3d on
% view(10,50)
% shading interp
% print(gcf,'-depsc2','-opengl',[path_figs,name,num2str(log10(1/Rho))],'-r300')
% print(gcf,'-dpdf','-opengl',[path_figs,name,num2str(log10(1/Rho))],'-r300')
% print(gcf,'-djpeg100','-opengl',[path_figs,name,num2str(log10(1/Rho))])

% figure(2)
% surf(X(3:Step_x:N),L,V(3:Step_x:N,:)')
% rotate3d on, view(0,70)
% shading interp

function dfdt = rhs(t,y,l)

u = y(1:N);
v = y(N+1:2*N);

u0 = u(1); u1 = u(N);
v0 = v(1); v1 = v(N);
cu = (Eps/H/l)^2;
cv = D/(Tau*(l*H)^2);

a = - u + u.^2./(v.*(1 + Kap*u.^2)); 
b = u.^2;

dudt = zeros(N,1);
dvdt = zeros(N,1);

dudt(1) = (u0-2*u(1)+u(2))*cu + H*(u(2)-u0)*0.5*cu/(R(1)+H) + a(1);
dudt(N) = (u(N-1)-2*u(N)+u1)*cu + H*(u1-u(N-1))*0.5*cu/(R(N)+H) + a(N);

dvdt(1) = (v0-2*v(1)+v(2))*cv + H*(v(2)-v0)*0.5*cv/(R(1)+H) + (-v(1) + b(1))/Tau;
dvdt(N) = (v(N-1)-2*v(N)+v1)*cv + H*(v1-u(N-1))*0.5*cv/(R(N)+H) + (-v(N) + b(N))/Tau;

for i=2:N-1
  dudt(i) = (u(i-1)-2*u(i)+u(i+1))*cu + H*(u(i+1)-u(i-1))*0.5*cu/(R(i)+H) + a(i);
  dvdt(i) = (v(i-1)-2*v(i)+v(i+1))*cv + H*(u(i-1)-u(i-1))*0.5*cv/(R(i)+H) + (-v(i) + b(i))/Tau;
end

dfdt = [dudt; dvdt];
end
end
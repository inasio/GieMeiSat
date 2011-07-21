function gms_stationary

path_figs = '/home/rozada/Temporary/Figures/GMSpaper/';
name = 'GMS_splitting';
L_init = 2; L_end = 2.3;
N = 1500;
T_range = [0,100];
Step_x = ceil(N/200);
L_all = linspace(L_init,L_end,200);
H = 1/(N+1);
X = (H:H:1-H)';
D = 1;
% Eps = 0.02;
Eps = 0.002;
Kap = 2.5;
Tau = 1;

X_l = 1/4;
U0 = 0.346*tanh(25*(X-X_l));
U0 = U0 + 0.346*tanh(25*(1-X_l-X));
V0 = 0.287 - 0.03*cos(2*pi*X);

[T0,Y0] = ode15s(@rhs, [0,100], [U0;V0],[],L_init);

Y = Y0(end,:)';
tic
for I = 2:length(L_all)
  I/length(L_all)
  [T,W] = ode15s(@rhs, T_range, Y(:,I-1),[],L_all(I));
%   W(end,1:N)'
  Y(:,I) = [W(end,1:N)';W(end,N+1:2*N)'];
end
toc

U = Y(1:N,:);
V = Y(2*N:end,:);

figure(1)
surf(X(3:Step_x:N),L_all,U(3:Step_x:N,:)')
rotate3d on
view(10,50)
shading interp
% print(gcf,'-depsc2','-opengl',[path_figs,name,num2str(log10(1/Rho))],'-r300')
% print(gcf,'-dpdf','-opengl',[path_figs,name,num2str(log10(1/Rho))],'-r300')
% print(gcf,'-djpeg100','-opengl',[path_figs,name,num2str(log10(1/Rho))])

for I=1:size(W,1)
  Param(I) = norm(W(I,1:2*N));
end

% figure(2)
% surf(X(3:Step_x:N),L,V(3:Step_x:N,:)')
% rotate3d on, view(0,70)
% shading interp

cd ~/UBC/Chamba_Michael/patterns_in_growing_domains/code/auto/gms
for I=1:size(U,2)
  NormU(I) = (dot(U(:,I),U(:,I))/N)^.5;
  NormV(I) = (dot(V(:,I),V(:,I))/N)^.5;
end

for I=1:4
    nom = ['branch',num2str(I)];
    nom = load(nom);
    figure(5)
    plot(nom(:,1),nom(:,4),'r')
    hold on
    plot(L_all,max(V),'--')
    figure(6)
    plot(nom(:,1),nom(:,3),'r')
    hold on
    plot(L_all,max(U),'--')
    figure(7)
    plot(nom(:,1),nom(:,5),'r')
    hold on
    plot(L_all,NormU,'--')
    figure(8)
    plot(nom(:,1),nom(:,6),'r')
    hold on
    plot(L_all,V(1,:),'--')
end

% fid = fopen(['GMS_rho_',num2str(Rho),'_N_',num2str(N),'.txt'],'w');
% maxU = max(U); maxV = max(V);
% for I = 1:length(L)
%   fprintf(fid,['%6.6f  %6.6f  %6.6f  %6.6f  %6.6f\n'],...
%     L(I),NormU(I),NormV(I),maxU(I),maxV(I));
% end
% fclose(fid)

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

dudt(1) = (u0-2*u(1)+u(2))*cu + a(1);
dudt(N) = (u(N-1)-2*u(N)+u1)*cu + a(N);

dvdt(1) = (v0-2*v(1)+v(2))*cv + (-v(1) + b(1))/Tau;
dvdt(N) = (v(N-1)-2*v(N)+v1)*cv + (-v(N) + b(N))/Tau;

for i=2:N-1
  dudt(i) = (u(i-1)-2*u(i)+u(i+1))*cu + a(i);
  dvdt(i) = (v(i-1)-2*v(i)+v(i+1))*cv + (-v(i) + b(i))/Tau;
end

dfdt = [dudt; dvdt];
end
end
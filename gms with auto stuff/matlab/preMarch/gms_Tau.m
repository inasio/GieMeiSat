function gms_Tau

% on [12,65] should be unstable, on [66,96] stable, [97,163] unstable
% on 82 (1), and Tau=330, oscillations split the solution

D = 1;
Eps = 0.02;
Kap = 2.5;
% Tau_all = [56:60, 61:0.2:75,76:100,101:3:301];
Tau_all = 380;
N = 300;
H = 1/(N+1);
X = (H:H:1-H)';
T = 0:4000;
noise = 15000000;
skipx = ceil(max(length(X))/250); skipt = ceil(max(T)/250);
% skipx = 1; skipt = 1;

% L1_all = load('../../data/L1param');
L2_all = load('../../data/L2param');
% Int = 12:164;
% Int = 204:266
Int = 82;
Temo = 1; clf()
for I = 1:length(Tau_all)
  Tau = Tau_all(I);
%   L = L1_all(Int-11); 
  L = 1.6;   % at this value there are permanent oscillations
  Nom_u = ['../../data/sol1_u',num2str(Int)];
  Nom_v = ['../../data/sol1_v',num2str(Int)];
  Nom_x = ['../../data/sol1_x',num2str(Int)];
  U0 = load(Nom_u);
  V0 = load(Nom_v);
  X0 = load(Nom_x);
  Uend = spline(X0,U0,X) + max(U0)*randn(N,1)/noise;
  Vend = spline(X0,V0,X) + max(V0)*randn(N,1)/noise;

  options = odeset('RelTol',1e-8,'AbsTol',1e-9);
  [T,Y] = ode15s(@rhs, T, [Uend;Vend],[],L);
  U = Y(:,1:N)';
  V = Y(:,N+1:2*N)';
  Uend = U(:,end);
  Vend = V(:,end);
  Lap = laplacian;
  A = dfg_duv(Uend,Vend);
  [EVE,EVA] = eig(full(A + Lap));
  Eva_sorted_all = sort(diag(EVA),'ascend');
%   Eva_sorted_re = sort(real(diag(EVA)),'descend');
%   Eva_sorted_im = sort(imag(diag(EVA)),'descend');
%   Max_6_re(1:6,Temo) = Eva_sorted_re(1:6);
%   Max_6_im(1:6,Temo) = Eva_sorted_im(1:6);
%   MaxEva_re(Temo) = max(real(diag(EVA)));
%   MaxEva_im(Temo) = max(imag(diag(EVA)));
%   [MaxEva_re, MaxEva_im]
  wild = Eva_sorted_all(1:4);
  plot(wild,'o'), hold on, 
%   axis([-0.001,0.001,-0.05,0.05])
  grid on, pause(0.1)
  [r1,c1] = find(EVA == wild(1));
  [r2,c2] = find(EVA == wild(2));
  [r3,c3] = find(EVA == wild(3));
  [r4,c4] = find(EVA == wild(4));
  [r1, r2, r3, r4]
  Tau
end

figure(2)
subplot(2,1,1)
surf(X(1:skipx:end),T(1:skipt:end),U(1:skipx:end,1:skipt:end)')
rotate3d on, shading interp, view(-90, 90)
subplot(2,1,2)
surf(X(1:skipx:end),T(1:skipt:end),V(1:skipx:end,1:skipt:end)')
shading interp, view(-90, 90)
figure(3)
subplot(2,1,1)
plot(real(EVE(1:N,r1)),'r')
xlabel(num2str(EVA(r1,r1)),'fontsize',12)

subplot(2,1,2)
plot(real(EVE(1:N,r3)))
xlabel(num2str(EVA(r3,r3)),'fontsize',12)


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

function y = laplacian_v
  
  cv = D/(Tau*(L*H)^2);

  e = ones(N,1);
  a = spdiags([e -2*e e], -1:1,  N, N);
  y = cv*a;
  y(1,1) = -cv; y(N,N) = -cv;
end

function y = laplacian_u
  
  cu = (Eps/H/L)^2;

  e = ones(N,1);
  a = spdiags([e -2*e e], -1:1,  N, N);
  y = cu*a;
  y(1,1) = -cu; y(N,N) = -cu;
end

function y = a_nonlin(u,v)
  
  f = - u + u.^2./(v.*(1 + Kap*u.^2));
  g = ( - v + u.^2)/Tau;
  y = [f; g];
end

function y = dfg_duv(u,v)

  dfdu = (2*u)./(v.*(1+Kap*u.^2)) - (2*Kap*u.^3)./(v.*(1+Kap*u.^2).^2) - 1;
  dfdv = -u.^2./(v.^2.*(1 + Kap*u.^2));
  dgdu = 2*u/Tau;
  dgdv = -ones(N,1)/Tau;

  a11 = diag(dfdu);
  a12 = diag(dfdv);
  a21 = diag(dgdu);
  a22 = diag(dgdv);

  x = [a11, a12; a21, a22];
  y = sparse(x);
end

function y = laplacian

  a = laplacian_u;
  b = laplacian_v;
  y = [a, zeros(N); zeros(N), b];
end
end
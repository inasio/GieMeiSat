function gms_eigvals

% on [12,42] should be unstable, on [43,96] stable, [97,163] unstable
path_figs = '/home/rozada/Temporary/Figures/GMSpaper/';
name = 'GMS_initial_cond';
opts.maxit = 10000;
opts.tol = 1e-6;
D = 1;
Eps = 0.02;
Kap = 2.5;
Tau = 1;
N = 150;
H = 1/(N+1);
X = (H:H:1-H)';

L1_all = load('../data/L1param');
% L2_all = load('../data/L2param');
Int = 12:163;
% Int = 204:266
% Int = 85;
Temo = 1;
for I = Int
%   (I-min(Int))/(max(Int) - min(Int))
  L = L1_all(I+1-12);
  I
  Nom_u = ['../data/sol1_u',num2str(I)];
  Nom_v = ['../data/sol1_v',num2str(I)];
  Nom_x = ['../data/sol1_x',num2str(I)];
  U0 = load(Nom_u);
  V0 = load(Nom_v);
  X0 = load(Nom_x);
  U = spline(X0,U0,X);
  V = spline(X0,V0,X);
%   U = ones(N,1)*0.5603;
%   V = ones(N,1)*0.3139;
  
%   YY = [I,norm(laplacian_v*V + (- V + U.^2))/Tau,...
%   norm(laplacian_u*U - U + U.^2./(V.*(1 + Kap*U.^2)))];
  Lap = laplacian;
  A = dfg_duv(U,V);
  [EVE,EVA] = eig(full(A + Lap));
%   [EVE,EVA] = eigs(A+Lap,1,'LR',opts);  
  Eva_sorted = sort(real(diag(EVA)),'descend');
  Max_6(1:6,Temo) = Eva_sorted(1:6);
  MaxEva(Temo) = max(real(diag(EVA)));
  [r,c] = find(EVA == MaxEva(Temo));
%   plot(EVE(1:N,c))
%   pause(0.1)
  Temo = Temo+1;
  
  [T,Y] = ode15s(@rhs, [0,100], [U;V],[],L);
end
figure(1)
plot(Int,MaxEva), grid on
figure(2)
plot(Int,Max_6), grid on
figure(3)
plot(L1_all,Max_6), grid on
figure(4)
plot(L1_all,MaxEva), grid on
xlabel('L','fontsize',15)
ylabel('max(real(\lambda))','fontsize',15)


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
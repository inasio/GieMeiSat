function gms_critical_tau

% on [12,65] should be unstable, on [66,96] stable, [97,163] unstable
% on 82 (1), and Tau=330, oscillations split the solution

D = 1;
Eps = 0.02;
Kap = 2.5;
N = 300;
H = 1/(N+1);
X = (H:H:1-H)';
T = 0:4000;
Noise = 15000000;
L1_all = load('../data/L1param');
Int = 82;
L = L1_all(Int-11)
Nom_u = ['../data/sol1_u',num2str(Int)];
Nom_v = ['../data/sol1_v',num2str(Int)];
Nom_x = ['../data/sol1_x',num2str(Int)];
U0 = load(Nom_u);
V0 = load(Nom_v);
X0 = load(Nom_x);

Tau_Critical = fzero(@eval_tau,280)

function y = eval_tau(Tau)

Uend = spline(X0,U0,X) + max(U0)*randn(N,1)/Noise;
Vend = spline(X0,V0,X) + max(V0)*randn(N,1)/Noise;

options = odeset('RelTol',1e-8,'AbsTol',1e-9);
[T,Y] = ode15s(@rhs, T, [Uend;Vend],[],L,Tau);
U = Y(:,1:N)';
V = Y(:,N+1:2*N)';
Uend = U(:,end);
Vend = V(:,end);
Lap = laplacian(Tau);
A = dfg_duv(Uend,Vend,Tau);
[EVE,EVA] = eig(full(A + Lap));
Eva_sorted_all = sort(diag(EVA),'ascend');
wild = Eva_sorted_all(1:4);
y = real(wild(3));
end

function dfdt = rhs(t,y,l,Tau)

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

function y = laplacian_v(Tau)
  
  cv = D/(Tau*(L*H)^2);

  e = ones(N,1);
  a = spdiags([e -2*e e], -1:1,  N, N);
  y = cv*a;
  y(1,1) = -cv; y(N,N) = -cv;
end

function y = laplacian_u(Tau)
  
  cu = (Eps/H/L)^2;

  e = ones(N,1);
  a = spdiags([e -2*e e], -1:1,  N, N);
  y = cu*a;
  y(1,1) = -cu; y(N,N) = -cu;
end

function y = a_nonlin(u,v,Tau)
  
  f = - u + u.^2./(v.*(1 + Kap*u.^2));
  g = ( - v + u.^2)/Tau;
  y = [f; g];
end

function y = dfg_duv(u,v,Tau)

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

function y = laplacian(Tau)

  a = laplacian_u(Tau);
  b = laplacian_v(Tau);
  y = [a, zeros(N); zeros(N), b];
end
end
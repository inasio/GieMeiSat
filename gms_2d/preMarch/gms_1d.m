function gms_explicit

path_figs = '/home/rozada/Temporary/Figures/GMSpaper/';
name = 'GMS_splitting';
L = 2;
N = 200;
% T = 1:100;
H = 1/(N+1);
X = (H:H:1-H)';
D = 1;
Eps = 0.02;
Kap = 2.5;
Tau = 1;

X_l = 1/4;
[U0,V0] = Initial_condition;
% U0 = 0.346*tanh(25*(X-X_l));
% U0 = U0 + 0.346*tanh(25*(1-X_l-X));
% V0 = 0.287 - 0.03*cos(2*pi*X);

[T0,Y0] = ode45(@rhs2, [0:0.1:10], [U0;V0]);

figure(1)
surf(Y0(:,1:N))
shading interp
figure(2)
plot(X,U0), hold on
plot(X,V0)
plot(X,Y0(1,1:N),'r')
plot(X,Y0(1,N+1:2*N),'c')

function dfdt = rhs2(t,y)

  u = y(1:N);
  v = y(N+1:2*N);

  u0 = u(1); u1 = u(N);
  v0 = v(1); v1 = v(N);
  cu = (Eps/H/L)^2;
  cv = D/(Tau*(L*H)^2);

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
function [u0,v0] = Initial_condition
  h = 0.3;
  xl = 0.3;
  xr = 1-xl;
  w_plus = 3.295209;
  denom = .9/Eps;
  uu = zeros(N,N); vv = zeros(N,N);
  for i = 1:N
    u0(i) = h*w_plus*0.425*(tanh(denom*(X(i)-xl))) + ...
      h*w_plus*0.425*(tanh(denom*(xr-X(i))));
    v0(i) = 0.312 - 0.035*cos(2*pi*X(i));
  end
end

end
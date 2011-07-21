function draw_single_mesa

%% system parameters (constants)
L_init = 1;
L = 2;  % domain length
N = 200;
Step_x = ceil(N/200); 
T_range = [0:1000]; 
H = 1/(N+1);
X = (H:H:1-H)';
D = 1;
Eps = 0.01;
Kap = 2.5;
Tau = 1;

%% initial conditions
X_l = 1/4;
U0 = 0.346*tanh(25*(X-X_l));
U0 = U0 + 0.346*tanh(25*(1-X_l-X));
V0 = 0.287 - 0.03*cos(2*pi*X);

[T0,Y0] = ode15s(@rhs, [0,20], [U0;V0], [], L_init);

% Y = [Y0(end,:)';L_init];
Y = Y0(end,:)';
[T,W] = ode15s(@rhs, T_range, Y, [], L);
U = W(:,1:N)';
V = W(:,N+1:2*N)';

% figure(1)
% surf(X(3:Step_x:N),T,U(3:Step_x:N,:)')
% rotate3d on
% shading interp
left = 0.245; lw = 2;
figure(2)
plot(X, U(:,end)')
hold on
plot(X, 10*(V(:,end)'-0.25))
% axis([-0.2 1.2 -0.3 1.2])
axis off
Ly1 = line([0 0], [-0.2 -.1],'linewidth',lw)
set(Ly1,'color','k')
Ly2 = line([left left], [-0.2 -.1],'linewidth',lw)
set(Ly2,'color','k')
Ly3 = line([1-left 1-left], [-0.2 -.1],'linewidth',lw)
set(Ly3,'color','k')
Ly4 = line([1 1], [-0.2 -.1],'linewidth',lw)
set(Ly4,'color','k')
Lx = line([-0.1 1], [-0.15 -0.15],'linewidth',lw)
set(Lx,'color','k')

text(-.015,-.25, '0','fontsize',15)
text(1-.015,-.25, 'L','fontsize',15)
text(left-.015,-.25, '\chi_l','fontsize',15)
text(1-left-.015,-.25, '\chi_r','fontsize',15)
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
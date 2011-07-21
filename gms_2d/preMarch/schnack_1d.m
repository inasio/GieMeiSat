function schnack_1d

Name = 'full_soln';
N = 50;
L = 8;
X=linspace(0,L(1),N);
% u0 = 0.9 + 0.01*rand(N,1);
% v0 = 1 + 0.01*rand(N,1);
% u0 = 1-.3*sech((x'-L(1))/.5).^2;
% v0 = .2+1.4*sech((x'-L(1)/2)*2).^2;
U0 = zeros(N,1); U0(ceil(end/2)) = 0.9;
V0 = zeros(N,1); V0(ceil(end/2)) = 1;

Y0(:,1) = [U0;V0];

% options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[T1,Y1] = ode15s(@fun_eval, [0:100], Y0);

figure(1)
plot(Y1(end,N+1:2*N))
Y1(end,N+1:2*N)
grid on

%% save variables to disk
% fid = fopen(['/home/rozada/Desktop/Matlab_tests/data_schnak/',name,'_L','.mat'],'w');
% fwrite(fid,L,'double');
% 
% fid = fopen(['/home/rozada/Desktop/Matlab_tests/data_schnak/',name,'_V','.mat'],'w');
% fwrite(fid,Y0(N+1,:),'double');

  function dfdt = fun_eval(t,y)

    u = y(1:N);
    v = y(N+1:2*N);

    a  = 0.9;
    b  = 0.01;  % something like O(\epsilon^2)
    du = 1;
    dv = 0.06;
    % u0 = u(1); u1 = u(N); % from Neumann boundary conditions
    % v0 = v(1); v1 = v(N);
    l2 = L^2;
    h  = 1/(N+1);
    cu = du/(l2*h^2);
    cv = dv/(l2*h^2);

    dudt = zeros(N,1);
    dvdt = zeros(N,1);

    dudt(1) = cu*(u(2) - u(1)) + a - u(1)*v(1)^2;
    dudt(N) = cu*(u(N-1) - u(N)) + a - u(N)*v(N)^2;

    dvdt(1) = cv*(v(2) - v(1)) + b + u(1)*v(1)^2 - v(1);
    dvdt(N) = cv*(v(N-1) - v(N)) + b + u(N)*v(N)^2 - v(N);

    for i=2:N-1
      dudt(i) = cu*(u(i-1)-2*u(i)+u(i+1)) + a - u(i)*v(i)^2;
      dvdt(i) = cv*(v(i-1)-2*v(i)+v(i+1)) + b + u(i)*v(i)^2 - v(i);
    end

    dfdt = [dudt; dvdt];
  end
end
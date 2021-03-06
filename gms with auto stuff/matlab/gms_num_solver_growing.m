function gms_num_solver_growing  % solve the system Au = b implicitly
%% constants

  N = 400;
  L = 2.5; % length of the domain
  L_end = 0.9 % for when domain is growing;
  H = 2*L/(N-1); % 100 grid points 
  X = -L:H:L; X = X';
  Ep = 0.01;
  D = 50; % diffusion coefficient ratio, in the v equation
  DD = Ep*D
  Tau = 25000;
  Kap = 1.0;
  Tend = 1e5; % integration time
  Dt = .1; % time step size, arbitrary for nonlinear problems
  
%% function definitions   
  F_gms = inline('u.^2./(v.*(1+kap*u.^2))','u','v','kap');
  G_gms = inline('u.^2','u','v');
  F_gms_full = inline('-u+u.^2./(v.*(1+kap*u.^2))','u','v','kap');
  G_gms_full = inline('-v+u.^2','u','v');
  Dfu = inline('- 1 + 2*u./(v.*(1+Kap*u.^2).^2)','u','v','Kap');
  Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
  Dgu = inline('2*u','u');
  Dgv = -1;
  
  UV_0(:,1) = 0.37*(tanh(.3*(.4*L-X)/Ep) + tanh(.3*(X+.4*L)/Ep));
  UV_0(:,2) = 0.3251*ones(N,1); % initial conditions


%% compute stationary solution
  [UU0,VV0,TT] = UV_stationary(UV_0,1,200, 1.0, Kap+0.06);
%   [UU,VV,TT] = UV_non_stationary([UU0(:,end),VV0(:,end)],Tau,Tend, Dt);
%   [UU,VV,TT] = UV_stationary([UU0(:,end),VV0(:,end)],Tau,Tend, Dt, Kap);
%   Y0(1:N,1) = UU0(:,end); 
  Y0(N+1:2*N,1) = VV0(:,end);
  Vspline = spline(X,UU0(:,end),X-0.01);
  Y0(1:N,1) = Vspline;
  opts = ['Reltolerance',1e-6];
  [TT,YY] = ode15s(@UV_stationary2, [0:round(Tend/500):Tend], Y0, opts, Kap);
  UU = YY(:,1:N)';

%% plotting
  figure(5)
  surfc(X,TT(1:round(length(TT)/500):end),0.01+UU(:,1:round(length(TT)/500):end)') 
% %   contour(X,TT(1:round(length(TT)/200):end),0.01+UU(:,1:round(length(TT)/200):end)') 
  shading interp, view(0,90)
% %   colormap('gray')
% %   surf(U_stat), shading interp

  function [uu,vv,t] = UV_stationary(uv_0,tau,t_end,dt,kap)  % y = [U,V,T]
    uu = zeros(length(uv_0),round(t_end/dt));
    vv = zeros(length(uv_0),round(t_end/dt));
    uu(:,1) = uv_0(:,1);
    vv(:,1) = uv_0(:,2);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    au = speye(N)*(1/dt+1) - (Ep/H)^2*lap;
    av = speye(N)*(tau/dt+1) - D/H^2*lap;
    
    for i = 2:round(t_end/dt)   
      bu = uu(:,i-1)/dt + F_gms(uu(:,i-1),vv(:,i-1),kap);
      uu(:,i) = au\bu;
      bv = tau*vv(:,i-1)/dt + G_gms(uu(:,i),vv(:,i-1));
      vv(:,i) = av\bv;
    end
    t = dt:dt:i*dt;
  end

  function [uu,vv,t] = UV_non_stationary(uv_0,tau,t_end,dt,kap)  % y = [U,V,T]
    uu = zeros(length(uv_0),round(t_end/dt));
    vv = zeros(length(uv_0),round(t_end/dt));
    uu(:,1) = uv_0(:,1);
    vv(:,1) = uv_0(:,2);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs

    for i = 2:round(t_end/dt)
      t = i*dt; 
      ll = L + (L_end-L)*t/t_end;
      [i/round(t_end/dt),ll]
      h = 2*ll/(N-1);
      au = speye(N)*(1/dt+1) - (Ep/h)^2*lap;
      av = speye(N)*(tau/dt+1) - D/h^2*lap;
      bu = uu(:,i-1)/dt + F_gms(uu(:,i-1),vv(:,i-1),kap);
      uu(:,i) = au\bu;
      bv = tau*vv(:,i-1)/dt + G_gms(uu(:,i),vv(:,i-1));
      vv(:,i) = av\bv;
    end
    t = dt:dt:i*dt;
  end

  function y = UV_non_stationary2(t,uv_0,kap) 
    uu = uv_0(1:N);
    vv = uv_0(N+1:2*N);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
 
    ll = L + (L_end-L)*t/Tend;
    h = 2*ll/(N-1);
    y(1:N,1) = (Ep/h)^2*lap*uu + F_gms_full(uu,vv,kap);
    y(N+1:2*N,1) = D/(Tau*h^2)*lap*uu + G_gms_full(uu,vv)/Tau;
  end

  function y = UV_stationary2(t,uv_0,kap) 
    uu = uv_0(1:N);
    vv = uv_0(N+1:2*N);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    y(1:N,1) = (Ep/H)^2*lap*uu + F_gms_full(uu,vv,kap);
    y(N+1:2*N,1) = (D/(Tau*H^2))*lap*vv + G_gms_full(uu,vv)/Tau;
  end

end
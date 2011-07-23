function gms_simple_1dplot  % plot simple mesa figures
%% constants

Title = ['data/basic_mesa_']
Kappa = [0.5,1,2.5,5];
  for I = 1:length(Kappa)
  N = 400;
  L = 1; % length of the domain
  H = 2*L/(N-1); % 100 grid points 
  X = -L:H:L; X = X';
  Ep = 0.01;
  D = 10; % diffusion coefficient ratio, in the v equation
  Tau = 1;
  Kap = Kappa(I);
  Tend = 100; % integration time
  Dt = .1; % time step size, arbitrary for nonlinear problems
  
%% function definitions   
  F_gms = inline('u.^2./(v.*(1+kap*u.^2))','u','v','kap');
  G_gms = inline('u.^2','u','v');
  Dfu = inline('- 1 + 2*u./(v.*(1+Kap*u.^2).^2)','u','v','Kap');
  Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
  Dgu = inline('2*u','u');
  Dgv = -1;
  
  UV_0(:,1) = 0.37*(tanh(.3*(.4*L-X)/Ep) + tanh(.3*(X+.4*L)/Ep));
  UV_0(:,2) = 0.3251*ones(N,1); % initial conditions

%% compute stationary solution
  [UU0,VV0,TT] = UV_stationary(UV_0,1,500, 1.0);
  [UU,VV,TT] = UV_stationary([UU0(:,end),VV0(:,end)],Tau,Tend, Dt);
%   Lam_nonstat = EVA_not_stationary([UU(:,end), VV(:,end)], L_range);
  
%% plotting
  figure(1)
  surfc(X,TT(1:round(length(TT)/500):end),UU(:,1:round(length(TT)/500):end)') 
%   contour(X,TT(1:round(length(TT)/200):end),0.01+UU(:,1:round(length(TT)/200):end)') 
  shading interp, view(0,90)
%   colormap('gray')
%   surf(U_stat), shading interp

  figure(2)
  plot(X, UU(:,end)'), hold on
  Usave = UU(:,end)'
  
  save([Title,'kappa',num2str(Kap),'.txt'],'Usave','-ASCII')
  end
save([Title,'X','.txt'],'X','-ASCII') 

  function [uu,vv,t] = UV_stationary(uv_0,tau,t_end,dt)  % y = [U,V,T]
    uu = zeros(length(uv_0),round(t_end/dt));
    vv = zeros(length(uv_0),round(t_end/dt));
    uu(:,1) = uv_0(:,1);
    vv(:,1) = uv_0(:,2);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    au = speye(N)*(1/dt+1) - (Ep/H)^2*lap;
    av = speye(N)*(tau/dt+1) - D/H^2*lap;

%     i = 2; norm_err = 1;
%     while norm_err > 1e-6
    for i = 2:round(t_end/dt)
        bu = uu(:,i-1)/dt + F_gms(uu(:,i-1),vv(:,i-1),Kap);
        uu(:,i) = au\bu;
        bv = tau*vv(:,i-1)/dt + G_gms(uu(:,i),vv(:,i-1));
        vv(:,i) = av\bv;
%         norm_err = norm(uu(:,i) - uu(:,i-1),'inf');
    end
    t = dt:dt:i*dt;
  end

end
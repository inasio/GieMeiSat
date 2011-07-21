function gms_1d_implicit  % solve the system Au = b
% goes with mesa_eigenvals2
%% constants

  N = 100;
  L = .5; % length of the domain
  H = 2*L/(N-1); % 100 grid points 
  X = -L:H:L; X = X';
  Ep = 0.02;
  D = 20; % diffusion coefficient ratio, in the v equation
  DD = Ep*D;
  Tau = 1;
  Kap = 5;
  Tend = 100; % integration time
  Dt = 0.1; % time step size, arbitrary for nonlinear problems
  L_range = [.5:0.05:2];
%   Tau_range = 1:200;
  
%% function definitions   
  F_gms = inline('u.^2./(v.*(1+kap*u.^2))','u','v','kap');
  G_gms = inline('u.^2','u','v');
  Dfu = inline('- 1 + 2*u./(v.*(1+Kap*u.^2).^2)','u','v','Kap');
  Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
  Dgu = inline('2*u','u');
  Dgv = -1;
  
  UV_0(:,1) = 0.37*(tanh(.3*(.4*L-X)/Ep) + tanh(.3*(X+.4*L)/Ep));
  UV_0(:,2) = 0.3251*ones(N,1); % initial conditions

%% compute stationary solution and associated eigenvalues
  [UU,VV,TT] = UV_stationary(UV_0,H);
  [Lam_all, U_stat, V_stat] = EVA_all(UV_0, L_range);
%   Lam_nonstat = EVA_not_stationary([UU(:,end), VV(:,end)], L_range);
  
%% plotting
  figure(1)
%   surf(X,TT,VV'), shading interp
  surf(U_stat), shading interp
  figure(2)
  hold on
  plot(L_range,Lam_all(1,:))
%   hold on
  plot(L_range,Lam_all(2,:),'r')


  function [uu,vv,t] = UV_stationary(uv_0,h)  % y = [U,V,T]
    uu(:,1) = uv_0(:,1);
    vv(:,1) = uv_0(:,2);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    au = speye(N)*(1/Dt+1) - (Ep/h)^2*lap;
    av = Tau*speye(N)*(1/Dt+1) - D/h^2*lap;

    i = 2; norm_err = 1;
    while norm_err > 1e-6
        bu = uu(:,i-1)/Dt + F_gms(uu(:,i-1),vv(:,i-1),Kap);
        uu(:,i) = au\bu;
        bv = Tau*vv(:,i-1)/Dt + G_gms(uu(:,i),vv(:,i-1));
        vv(:,i) = av\bv;
        norm_err = norm(uu(:,i) - uu(:,i-1),'inf'); i = i+1
    end
    t = Dt:Dt:(i-1)*Dt;
  end

  function y = EVA_stationary(uv_end,h)
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    ue = uv_end(:,1); ve = uv_end(:,2);
    a_eig = [diag(Dfu(ue,ve,Kap)),diag(Dfv(ue,ve,Kap));...
    diag(Dgu(ue)),-speye(N)];
    b_eig = [(Ep/h)^2*lap, zeros(N); zeros(N), D/h^2*lap];

    [eve,eva] = eig(full(a_eig + b_eig));
    [eva_sorted, index] = sort(real(diag(eva)),'descend');
    y = [eva(index(1),index(1)); eva(index(2),index(2))];
  end

  function y = EVA_not_stationary(uv_end, l_range)
% from a single stat solution compute it's 2 largest eigvals
% as one of the parameters changes
    for i=1:length(l_range)
      l = l_range(i); % length of the domain
      h = 2*l/(N-1); % 100 grid points 
      y(:,i) = EVA_stationary(uv_end, h);
    end
  end
      
  function [lam, u_stat, v_stat] = EVA_all(uv_0, l_range)
% track a stationary solution and compute it's 2 largest eigvals
    for i=1:length(l_range)
      l = l_range(i); % length of the domain
      h = 2*l/(N-1); % 100 grid points 
      [uu,vv,t] = UV_stationary(uv_0,h);
      uv_0 = [uu(:,end), vv(:,end)];
      u_stat(:,i) = uu(:,end); v_stat(:,i) = vv(:,end);
      lam(:,i) = EVA_stationary(uv_0, h);
    end
  end
    
end


%% save variables to disk
% fid = fopen([title,'_Us','.mat'],'w');
% fwrite(fid,Us,'double');
% fid = fopen([title,'_Vs','.mat'],'w');
% fwrite(fid,Vs,'double');
% fid = fopen([title,'_t','.mat'],'w');
% fwrite(fid,t,'double');
% fid = fopen([title,'_x','.mat'],'w');
% fwrite(fid,x,'double');
% fid = fopen([title,'_L','.mat'],'w');
% fwrite(fid,L,'double');

%% call variables from disk
% clear all % just to test that it's working
% title = 'crank_schnak';
% fid = fopen([title,'_Us','.mat'],'r');
% tempu = fread(fid,'double');
% fid = fopen([title,'_Vs','.mat'],'r');
% tempv = fread(fid,'double');
% fid = fopen([title,'_t','.mat'],'r');
% t = fread(fid,'double');
% fid = fopen([title,'_x','.mat'],'r');
% x = fread(fid,'double');
% fid = fopen([title,'_L','.mat'],'r');
% L = fread(fid,'double');

%% now rescale into matrices and plot them
% N = length(x);
% for j=1:length(t)
%     Us(:,j) = tempu(1+(j-1)*N:j*N);
%     Vs(:,j) = tempv(1+(j-1)*N:j*N);
% end

% skipt = ceil(length(t)/100); skipx = ceil(N/100);

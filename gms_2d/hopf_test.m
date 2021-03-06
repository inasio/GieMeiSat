function hopf_test
  LL = 0.5;
%   X = -LL:0.01:LL;
  N = 50;
  H = 2*LL/(N-1);
  X = linspace(-LL,LL,N)';
  Time = 0:.001:.1;

  K = 2*pi;
  Y0 = cos(K*X);
  [TT,YY] = ode15s(@rhs, Time, Y0);
  V = YY(:,1:N)';
  figure(1)
  surf(X,TT,V'), shading interp
  for I=1:length(TT)
    V_exact(:,I) = exp(-K^2*TT(I))*cos(K*X);
  end
  figure(2)
  surf(X,TT,V_exact'-V'), shading interp
  figure(3)
  surf(X,TT,V_exact'), shading interp
  
  function dvldt = rhs(t,y)
    w = y(1:N);
%     l1 = y(N+1);
%     l2 = y(N+2);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    dvldt(1:N,1) = lap*w/H^2;
%     dvldt(N+1,1) = Mu*spline(X,w,y(N+1));
%     dvldt(N+2,1) = - Mu*spline(X,w,y(N+2));
  end

  function y = heavi(x)
    y = (abs(x)>Delta).*(x>0) + ...
        (abs(x)<=Delta).*(0.5 + x/(2*Delta) + sin(pi*x/Delta)/(2*pi));
  end

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

end
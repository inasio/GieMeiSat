function gms_2d

  N = 200;
  H = 1/(N-1);
  X = 0:H:1;
  D = 1;
  Eps0 = 0.01;
  Kappa = 1.92;
  Tau = 1;
  L = Laplacian_matrix;
  Domain_width = 4;
  
  Dt = 5;
  Tend = 300;
  T = 0:Dt:Tend;
  loadfile = 'data/basic1mesa'; savefile = 'data/testy';
  [U0,V0] = Load_mesas(2);
%   [U0, V0] = Initial_condition;
%   [U0,V0] = Load_data(loadfile);
  Yu_1 = U0; Yv_1 = V0;
  Yu_1 = Yu_1+rand(N^2,1)*0.05;
  Yv_1 = Yv_1+rand(N^2,1)*0.05;
  Lu = Umatrix;
  Lv = Vmatrix;
  Tnow = 0; I = 1; Tol = 1;
  while (Tnow < Tend) & (Tol>1e-3)
    RHSu = Yu_1 + Dt*Yu_1.^2./(Yv_1.*(1+...
      Kappa*Yu_1.^2));
    Yu_2 = Lu\RHSu;
    RHSv = Yv_1 + Dt*Yu_2.^2;
    Yv_2 = Lv\RHSv;
    pause(0.01)
    figure(1)
    mesh(X,X,reshape(Yu_2,N,N))
    shading interp, view(90,90), 
%     axis([0 1 0 1 0 1])
    Tol = norm(Yu_2 - Yu_1,2);
%     Dt = (Dt+0.1)*(Tol(I)<0.001) + Dt*0.8*(Tol(I)>10) + ...
%       Dt*(Tol(I)>0.01)*(Tol(I)<10);
    Tnow = Tnow + Dt;
    [Tnow/Tend, Dt, Tol]
    I = I + 1;
    Yu_1 = Yu_2; Yv_1 = Yv_2; 
  end
  figure(1)
  surf(X,X,reshape(Yu_2,N,N))
  shading interp, view(-76,14), axis([0 1 0 1 0 1])

  figure(3)
  plot(X,U0(1:N),'r'), hold on, plot(X,V0(1:N))
  plot(X,Yu_1(1:N),'r--'), plot(X,Yv_1(1:N),'--')
  Save_data(Yu_1,Yv_1,savefile)
  
  function l = Laplacian_matrix
    m1 = speye(N);
    data = ones(N,3); data(:,2) = -2;
    m2 = spdiags(data,-1:1,N,N);
    l = kron(m1,m2) + kron(m2,m1);
    for i = 1:N
      l(i,N+i) = 2;
      l(N*(N-1)+i,N*(N-2)+i) = 2;
      l(N*i,N*i-1) = 2;
      l(N*(i-1)+1,N*(i-1)+2) = 2;
    end
    l = l/H^2;
  end

  function [u0,v0] = Initial_condition
    h = 0.3;
    xl = 0.2;
    xr = 0.8;
    w_plus = 3.295209;
    denom = 1/Eps0;
    uu = zeros(N,N); vv = zeros(N,N);
    for i = 1:N
      for j = 1:N
        uu(i,j) = h*w_plus*0.425*(tanh(denom*(X(i)-xl))) + ...
          h*w_plus*0.425*(tanh(denom*(xr-X(i))));
        vv(i,j) = 0.312 - 0.035*cos(2*pi*X(i));
      end
    end
    u0 = reshape(uu,N^2,1); v0 = reshape(vv,N^2,1);
  end

  function a = Umatrix
    a = (Dt + 1)*speye(N^2) - Dt*Eps0^2*L/Domain_width^2;
  end

  function a = Vmatrix
    a = (Dt + Tau)*speye(N^2) - Dt*D*L/Domain_width^2;
  end

  function Save_data(u,v,filename)
    fid = fopen([filename 'u'],'w');
    fwrite(fid,u,'double');
    fid = fopen([filename 'u'],'w');
    fwrite(fid,v,'double');
  end
  
  function [u0,v0] = Load_data(filename)
    fid = fopen([filename 'u'],'r');
    u0 = fread(fid,'double');
    fid = fopen([filename 'v'],'r');
    v0 = fread(fid,'double');
  end

  function [u0,v0] = Load_mesas(n)
    [u,v] = Load_data('data/basic1mesaD10');
    ns = round(sqrt(length(u)));
    u1d = u(1:ns); v1d = v(1:ns);
    x1d = linspace(0,1,length(u1d));
    x = mod(X*n,1);
    for i = 1:N
      uu(:,i) = spline(x1d,u1d,x);
      vv(:,i) = spline(x1d,v1d,x);
    end
    u0 = reshape(uu,N^2,1); v0 = reshape(vv,N^2,1);
  end

end
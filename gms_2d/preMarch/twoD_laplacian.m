function twoD_laplacian

  for I = 1:25
    I
    Ene(I) = 20*I;
    N = Ene(I);
    H = 1/(N-1);
    X = 0:H:1;
    L = laplacian_matrix;
    Y = initial_condition;
    Z = L*Y;
    Feo(I) = norm(Z + 2*pi^2*Y,inf);
    A = speye(N^2) - 0.1*L;
    tic;
    A\Y;
    Time(I) = toc;
  end
  figure(1)
  loglog(Ene, Feo)
  hold on
  X = 1:500; Y = X.^(-2);
  loglog(X,Y)
  figure(2)
  loglog(Ene,Time)
  hold on
  X = 1:500; Y = X.^2/1000;
  loglog(X,Y)
%   figure(2)
%   surf(reshape(Z,N,N))
%   shading interp
%   figure(3)
%   surf(reshape(-2*pi^2*Y,N,N))
%   shading interp

  function l = laplacian_matrix

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

  function y = initial_condition
    
    for i = 1:N
      for j = 1:N
        y1(i,j) = cos(pi*X(i))*cos(pi*X(j));
      end
    end
    y = reshape(y1,N^2,1);
  end
end
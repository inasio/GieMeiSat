function react_diff_test

  N = 100;
  H = 1/(N-1);
  X = 0:H:1;
  L = laplacian_matrix;
  Y = initial_condition;
  NumIts = 100;
  for I = 1:NumIts
    Dt = 1/NumIts;
    M = speye(N^2) - Dt*L;
    Y(:,I+1) = M\Y(:,I); I
    surf(reshape(Y(:,I+1),N,N))
    shading interp
%     axis([0 N 0 N -1 1])
%     pause(0.1)
  end
  sosol = norm(Y(:,end) - Y(:,1)*exp(-2*pi^2))

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
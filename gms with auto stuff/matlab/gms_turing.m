function gms_turing

K = 2.5;
Tau = 1;
Ep = 0.02;
D = 1;
L = 1;
Du = Ep^2/L^2;
Dv = D/Tau/L^2;
N = 8;
Kap = 2*N*pi;

Bla = fzero(@stat_sol,1)

Turing = fzero(@det_DA,1)

function f = F_uv(u,v)

  f = - u + u.^2/(v.*(1 + K*u.^2));
end

function g = G_uv(u,v)
  
  g = (-v + u.^2)/Tau;
end

function y = stat_sol(u)
  
  y = u + K*u.^3 - 1;
end

function y = A_lin(u,v)
  
  y = [-1 + 2.*u./((1+K*u.^2).*v) - 2*K*u.^3.*v/((1+...
    K*u.^2).*v).^2, -u.^2./(v.^2.*(1+K*u.^2)); 2*u/Tau,...
    -1/Tau];
end

function y = det_DA(l)
  
  u0 = fzero(@stat_sol,1);
  v0 = u0.^2;
  a = A_lin(u0,v0);
  
  y = l.^4.*det(a) - l.^2*Kap^2*(Ep^2*a(2,2) + ...
    D*a(1,1)/Tau) + Kap^4*Ep^2*D/Tau;
end
end
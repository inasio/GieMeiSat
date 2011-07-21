function threshold_L
% estimate the domain length at which the mesa breaks into two

K = 2.5;
D = 1;
L = 1;

U_0 = 1/sqrt(K);
V_0 = 1/(2*sqrt(K));
U_c = 1.515/sqrt(K);
V_c = 0.4597/sqrt(K);

g_gms(U_c)
Term1 = atanh(sqrt(2*f_gms(U_c))/V_c)
Term2 = - sqrt(2*f_gms(U_c))/g_gms(U_c)
Term3 = - quad(@integrand_3_gms,U_0,U_c)
L = sqrt(D)*(Term1 + Term2 + Term3)
Lil_L = sqrt(D)*(Term2 + Term3);
[.5-Lil_L/2,Lil_L/2+.5]
2*L
% U_gms = U_0 + 0.01:0.01:2;
% for I = 1:length(U_gms)
%   F_gms(I) = quad(@f_integrand_gms,U_0,U_gms(I));
% end
% plot(U_gms,F_gms)

% f_gms(1.502)

function y = g_gms(s)

  y = s.*(1-s.*(1+K*s.^2))./(1+K*s.^2);
end

function y = h_prime_gms(s)

  y = (1 - K*s.^2)./(1 + K*s.^2).^2;
end

function y = g_prime_gms(s)

% esto esta mega chafa
%   syms u;
%   z = u*(1-u*(1+K*u^2))/(1+K*u^2);
%   z_prime = diff(z);
%   y = subs(z_prime,s);
  y = (2*K*s.^2.*(s.*(K*s.^2 + 1) - 1))./(K*s.^2 + 1).^2 -...
    (s.*(3*K*s.^2 + 1))./(K*s.^2 + 1) - (s.*(K*s.^2 + 1) - 1)/(K*s.^2 + 1);
end

function y = f_integrand_gms(s)
  
  y = g_gms(s).*h_prime_gms(s);
end

function y = f_gms(s)
  
  y = zeros(1,length(s));
  for i = 1:length(s)
    y(i) = quad(@f_integrand_gms,U_0,s(i));
  end
end

function y = integrand_3_gms(s)
  
  y = g_prime_gms(s).*sqrt(2*f_gms(s))./(g_gms(s).^2);
end
end
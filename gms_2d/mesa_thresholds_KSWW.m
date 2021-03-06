function mesa_thresholds_KSWW
% plot figure 24 on KSWW, stability of the strip GMS
figure(1), clf

for D = [10, 8, 4]
Ep = 0.0025;
% D = 4;
DD = D*Ep;
Wplus = 3.295209
Kap = 2;
Beta = 1.49882;
% L = 1;
B0 = 0.211376;
Uplus = sqrt(B0/Kap)*Wplus;
V0 = sqrt(B0/Kap);
[X,Fval,Flag] = fsolve(@f_minimizer_fzero,[V0+2.3,Uplus+2.3])  
Uplus = X(2);
VV = X(1);
Wplus = Uplus/sqrt(B0/Kap);
L = 1/(VV*Wplus^2);

Eigenmode = 0.0:0.1:10;
for I = 1:length(Eigenmode)
  
  M = Eigenmode(I);  %eigenmode
 
  Theta_minus = (M^2 + Ep/DD)^(1/2);
  Theta_plus = (M^2 + (Ep/DD)*(1+2*Wplus/(L*(Wplus-2))))^(1/2);
  Theta_minus = M; Theta_plus = M;
  Alpha = 2*Beta*L*DD/Wplus^2;
  Sigma_plus = (Theta_plus*tanh(Theta_plus*L/2) + ...
    Theta_minus*tanh(Theta_minus*(1-L)/2))^-1;
  Sigma_minus = (Theta_plus*coth(Theta_plus*L/2) + ...
    Theta_minus*tanh(Theta_minus*(1-L)/2))^-1;
  
  La_plus(I) = (Ep^2/Alpha)*(-Alpha*M^2 + L*(1-L)/2 - Sigma_plus);
  La_minus(I) = (Ep^2/Alpha)*(-Alpha*M^2 + L*(1-L)/2 - Sigma_minus);
%   La_plus(I) = - Sigma_plus;
%   La_minus(I) = - Sigma_minus;
end

figure(1)
plot(Eigenmode,La_minus,'b'), 
hold on
plot(Eigenmode,La_plus,'--'), 
axis([0 10 -0.0013 0.00005])
% plot(Eigenmode,La2)
end


  function z = f_minimizer_fzero(x)
    v0 = x(1);  
    uplus = x(2);
    f_of_uv = inline('- u + u.^2./(v.*(1 + k*u.^2))','u','v','k');
    z(1) = quad(f_of_uv,0,uplus,[],[],v0,Kap);
    z(2) = f_of_uv(uplus,v0,Kap);
  end
end
function mesa_eigenvals2
%% consider the case M =/= 0, Tau =/= 0, for only one mesa

  param = 0.5:0.05:20
  Lam_p0 = -0.1;
  Lam_m0 = -0.1;
  for I=1:length(param)
%     M = param(I);
    M = 0;
    Ep = 0.02;
    D = 20;
    Tau = 1;
%     Tau = param(I);
    DD = D*Ep;
    Kap = 5;
    Beta = 1.49882;
    B0 = 0.211376;
    Wplus = 3.295209;
    Uplus = sqrt(B0/Kap)*Wplus;
    Uminus = 0;
    V0 = sqrt(B0/Kap)
    LL = param(I);
%     LL = 0.5;  
1
    F_gms = inline('-u + u.^2./(v.*(1+Kap*u.^2))','u','v','Kap');
    G_gms = inline('- v + u.^2','u','v');
    Dfu = inline('- 1 + 2*u./(v.*(1+Kap*u.^2).^2)','u','v','Kap');
    Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
    Dgu = inline('2*u','u');
    Dgv = inline('-1','v');

%     [X,Fval,Flag] = fsolve(@f_minimizer_fsolve,[0.1,1],...
%       optimset('TolFun',1e-9))
%     Uplus = X(2);
%     V0 = X(1);
%     Wplus = Uplus/V0;

    Ll = G_gms(Uminus, V0)*LL/(G_gms(Uminus, V0) - G_gms(Uplus, V0));
    K0 = Beta*V0^2/(quad(Dfv,Uminus,Uplus,[],[],V0,Kap));
    
    Kap_p = - (Dgv(V0)-Dfv(Uplus,V0,Kap)/Dfu(Uplus,V0,Kap)*Dgu(Uplus));
    Kap_m = - (Dgv(V0)-Dfv(Uminus,V0,Kap)/Dfu(Uminus,V0,Kap)*Dgu(Uminus));
    Va = -DD*K0*Ll/(G_gms(Uminus, V0)*LL);
    [X,Fval,Flag] = fsolve(@lambda_pm,[Lam_p0,Lam_m0]);
    Lam_p(I) = X(1);
    Lam_m(I) = X(2);
    isreal(X(1))
    Lam_p0 = X(1); Lam_m0 = X(2);
  end
  figure(2)
  plot(param,Lam_p,'r--')
  hold on
%   axis([0 10 -0.00012 0.00004])
  plot(param,Lam_m,'--')
  hold on
  
  function lam = lambda_pm(x)
    lam_p = x(1);
    lam_m = x(2);
    
    sig_p = sqrt(M^2+Ep/DD*(Kap_p + Tau*lam_p));
    sig_m = sqrt(M^2+Ep/DD*(Kap_m + Tau*lam_m));
%     sig_p = M;
%     sig_m = M;
    om_p = (sig_m*tanh(sig_m*(LL-Ll)) + sig_p*tanh(Ll*sig_p))^(-1);
    om_m = (sig_m*tanh(sig_m*(LL-Ll)) + sig_p*coth(Ll*sig_p))^(-1);
    lam(1) = lam_p - Ep^2/Va*(om_p + (1-LL/Ll)*Ll^2/LL - M^2*Va);
    lam(2) = lam_m - Ep^2/Va*(om_m + (1-LL/Ll)*Ll^2/LL - M^2*Va);
  end

  function z = f_minimizer_fsolve(x)
    v0 = x(1); 
    uplus = x(2);
    f_of_uv = inline('- u + u.^2./(v.*(1 + k*u.^2))','u','v','k');
    z(1) = quad(f_of_uv,0,uplus,[],[],v0,Kap);
    z(2) = f_of_uv(uplus,v0,Kap);
  end
end
function mesa_eigenvals3
%% calculates the M vs L and l vs L graphs for the eigenvals

Param = 0.05:0.05:2;
Yp = 1; Ym = 1;
  for I = 2:length(Param)
  
    Ep = 0.02;
    D = 1;
    Tau = 1;
    DD = D*Ep;
    Kap = 2;
    B0 = 0.211376;
    Wplus = 3.295209;
    Uplus = sqrt(B0/Kap)*Wplus;
    Uminus = 0;
    V0 = sqrt(B0/Kap);
%     LL = 0.5;
    LL = Param(I)

    F_gms = inline('-u + u.^2./(v.*(1+Kap*u.^2))','u','v','Kap');
    G_gms = inline('- v + u.^2','u','v');
    Dfu = inline('- 1 + 2*u./(v.*(1+Kap*u.^2).^2)','u','v','Kap');
    Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
    Dgu = inline('2*u','u');
    Dgv = inline('-1','v');

%     Gplus = Param(I);
    Gplus = G_gms(Uplus, V0);
    Ll = G_gms(Uminus, V0)*LL/(G_gms(Uminus, V0) - Gplus);
    LLl(I) = Ll;
    Beta = 1.49882;
    K0 = Beta*V0^2/(quad(Dfv,Uminus,Uplus,[],[],V0,Kap));
    Va = -DD*K0*Ll/(G_gms(Uminus, V0)*LL);

    Yp(I) = fzero(@critical_modes_p,Yp(I-1));
%     Ym(I) = fzero(@critical_modes_m,Ym(I-1));
  end
  
  plot(LLl,Yp,'r'), hold on, plot(LLl,Ym)

  function y = omega_approx_p(m)
    y = (m*tanh(m*(LL-Ll)) + m*tanh(Ll*m))^(-1);
  end

  function y = omega_approx_m(m)
    y = (m*tanh(m*(LL-Ll)) + m*coth(Ll*m))^(-1);
  end

  function dy = domega_dm_p(m)
    small = 1e-8;
    dy = (omega_approx_p(m+small) - omega_approx_p(m-small))/(2*small);
  end

  function dy = domega_dm_m(m)
    small = 1e-4;
    dy = (omega_approx_m(m+small) - omega_approx_m(m-small))/(2*small);
%     a = LL-ll; b = Ll;
%     dy = -
  end

  function y = critical_modes_p(m)
    y = m.*domega_dm_p(m)/2 - omega_approx_p(m) - (1-LL/Ll)*Ll^2/Ll
  end

  function y = critical_modes_m(m)
    y = m.*domega_dm_m(m)/2 - omega_approx_m(m) - (1-LL/Ll)*Ll^2/Ll;
  end

end
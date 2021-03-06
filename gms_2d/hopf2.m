function hopf2
% estimation of the eigenvalues, use with hopf1.m and 
% gms_num_solver.m

  Kap = 1;
  D = 50; % diffusion coefficient ratio, in the v equation
  Ep = 0.01;
  DD = Ep*D
  G_gms = inline('- v + u.^2','u','v');
  Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
  Beta = 1.49882;
  B0 = 0.211376;
  Wplus = 3.295209;
  Uplus = sqrt(B0/Kap)*Wplus;
  Uminus = 0;
  V0 = sqrt(B0/Kap);
  K0 = Beta*V0^2/(quad(Dfv,Uminus,Uplus,[],[],V0,Kap));
  Delta1 = 0.01; Delta2 = 0.01;
  Tau1 = 1000; Tau2 = 1000; 
  G_minus = G_gms(Uminus,V0); G_plus = G_gms(Uplus,V0);
  Ll_ratio = (G_minus - G_plus)/G_minus
  
  plot_2evals(Kap,Ll_ratio)
%   plot_evals_kappa_ratio(0.4:0.1:1)
  
  function plot_evals_kappa_ratio(kappa_ratio)
     
    for i = 1:length(kappa_ratio)
      kap = kappa_ratio(i);
      uplus = sqrt(B0/kap)*Wplus;
      v0 = sqrt(B0/kap);
      k0 = Beta*v0^2/(quad(Dfv,Uminus,uplus,[],[],v0,kap));
      g_minus = G_gms(Uminus,v0); g_plus = G_gms(uplus,v0);
      l_ratio = (g_minus - g_plus)/g_minus
      lL_Ll = l_ratio;
      Va = -DD*k0/(G_gms(Uminus, v0)*lL_Ll);
      L = 0.1:0.01:8;
      for I = 1:length(L)
        LL = L(I);    
        Ll = LL/lL_Ll;
        Delta1(I+1) = fzero(@Plus_1, Delta1(I),[],LL,Ll);
        Tau1(I+1) = fzero(@Minus_1, Tau1(I), [], Delta1(I+1),LL,Ll);
        Delta2(I+1) = fzero(@Plus_2, Delta2(I),[],LL,Ll);
        Tau2(I+1) = fzero(@Minus_2, Tau2(I), [], Delta2(I+1),LL,Ll);  
      end

      figure(4)
      semilogy(L,-Tau1(2:end)*Va*DD/Ep^3), hold on
      semilogy(L,-Tau2(2:end)*Va*DD/Ep^3,'r--')
    end
%     title(['D = ',num2str(D)'],'fontsize',20), grid on
    xlabel('\kappa','fontsize',20), ylabel('\tau','fontsize',20)
    legend('\lambda_+','\lambda_-')  
%     title(['L/l = ',num2str(l_ratio)])
%     h = legend('\kappa','\lambda');
%     set(h,'fontsize',30)
  end

  function plot_2evals(Kap,LL_Ll)
% solve equations 18a and 18b in the HB notes (find the hopf 
% point in the GMS system

    Va = -DD*K0/(G_gms(Uminus, V0)*LL_Ll);
    L = 0.1:0.01:18;
    for I = 1:length(L)
      LL = L(I);    
      Ll = LL/LL_Ll;
      Delta1(I+1) = fzero(@Plus_1, Delta1(I),[],LL,Ll);
      Tau1(I+1) = Minus_1(Delta1(I+1),LL,Ll);
      Delta2(I+1) = fzero(@Plus_2, Delta2(I),[],LL,Ll);
      Tau2(I+1) = Minus_2(Delta2(I+1),LL,Ll);  
    end

    %% save variables to disk
    Tau_p = -Tau1(2:end)'*Va*DD/Ep^3;
    Tau_m = -Tau2(2:end)'*Va*DD/Ep^3;
    Lam_p = Ep*Delta1(2:end)'*DD./Tau_p;
    Lam_m = Ep*Delta2(2:end)'*DD./Tau_m;
    El = L';
    Title = ['../gms with auto stuff/data/kappa_',num2str(Kap)]

    
%     save([Title,'_lambda_plus','.txt'],'Lam_p','-ASCII')
%     save([Title,'_lambda_minus','.txt'],'Lam_m','-ASCII')
%     save([Title,'_tau_plus','.txt'],'Tau_p','-ASCII')
%     save([Title,'_tau_minus','.txt'],'Tau_m','-ASCII')
%     save([Title,'_L','.txt'],'El','-ASCII')
    
    figure(3)
    semilogy(L,Tau_p), hold on
    semilogy(L,Tau_m,'r--')
    grid on
    xlabel('L','fontsize',18), ylabel('\tau','fontsize',18)
    h = legend('\tau_+','\tau_-')
    set(h,'fontsize',20)
    title(['L/l = ',num2str(LL_Ll)],'fontsize',20)
    
    figure(4)
    plot(L,Lam_p), hold on
    plot(L,Lam_m,'r--')
    grid on
    xlabel('L','fontsize',18), ylabel('\lambda','fontsize',18)
    h = legend('\lambda_+','\lambda_-')
    set(h,'fontsize',20)
    title(['L/l = ',num2str(LL_Ll)],'fontsize',20)
    
  end
  
  
  function y = Plus_1(delta,LL,Ll)
    a = sqrt(i*delta);
    y = (1-LL/Ll)*Ll^2/LL + ...
      real((a*(tanh(a*(LL-Ll)) + tanh(a*Ll)))^(-1));
  end

  function tau = Minus_1(delta,LL,Ll)
    a = sqrt(i*delta);
    tau = - delta/imag((a*(tanh(a*LL-Ll)+tanh(a*Ll)))^(-1));
  end
   
  function y = Plus_2(delta,LL,Ll)
    a = sqrt(i*delta);
    y = (1-LL/Ll)*Ll^2/LL + ...
      real((a*(tanh(a*(LL-Ll)) + coth(a*Ll)))^(-1));
  end

  function tau = Minus_2(delta,LL,Ll)
    a = sqrt(i*delta);
    tau = - delta/imag((a*(tanh(a*LL-Ll)+coth(a*Ll)))^(-1));
  end
   
end
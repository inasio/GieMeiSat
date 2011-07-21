function hopf1
% compute solutions to dv1/dt and dv2/dt based on hopf notes
% compare to full numerical solution in gms_num_solver.m and
% to hopf2.m
  path_figs = '/home/rozada/UBC/Chamba_Michael/patterns_in_growing_domains/code/gms with auto stuff/images/';
  name = 'l1_l2_';
  
  LL = 2.5;
%   X = -LL:0.01:LL;
  N = 300;
  X = linspace(-LL,LL,N)';
  H = 2*LL/(N-1);
  Ep = 0.01;
  TTime = linspace(0,1e4,1500);
  Time = TTime*Ep^2;
  Delta = (Time(2)-Time(1))*5;
  Tau = 25000;
%   Tau_0 = .0001;   % this is tilde tau_0 = Ep^3 tau
  Tau_0 = Ep^3*Tau;
  DD = 50*Ep;
  DD_0 = DD/Tau_0;
  Kap = 1;
  
  Beta = 1.49882; B0 = 0.211376; Wplus = 3.295209;
  Uplus = sqrt(B0/Kap)*Wplus; Uminus = 0; V0 = sqrt(B0/Kap);
%   F_gms = inline('-u + u.^2./(v.*(1+Kap*u.^2))','u','v','Kap');
  G_gms = inline('- v + u.^2','u','v');
  Dfv = inline('- u.^2./(v.^2.*(1+Kap*u.^2))','u','v','Kap');
  G_minus = G_gms(Uminus,V0); G_plus = G_gms(Uplus,V0);
  Ll_eq = LL*G_minus/(G_minus - G_plus);
  K0 = Beta*V0^2/(quad(Dfv,Uminus,Uplus,[],[],V0,Kap));
  Mu = G_minus/(K0*Tau_0);
  
  nudge = 0.1;
  Y0 = V0*Tau_0*(1+0.02*randn(N+2,1))/G_minus; 
  Y0(N+1) = Ll_eq + nudge; Y0(N+2) = -Ll_eq + nudge;
  %% the good stuff (solver) below
  [TT,YY] = ode15s(@rhs, Time, Y0);
  
  V = YY(:,1:N)';
  L1 = YY(:,N+1)
  L2 = YY(:,N+2);
  
  for I = 1:length(TT)
    Int_V_x(I) = sum(YY(I,:))*H;
  end
  Check = (Int_V_x(2:end) - Int_V_x(1:end-1))/(TT(2)-TT(1)) - ...
    LL*(2 + (L2(2:end)' - L1(2:end)')/Ll_eq)
  norm((Int_V_x(2:end) - Int_V_x(1:end-1))/(TT(2)-TT(1)),inf)
  norm(LL*(2 + (L2(2:end)' - L1(2:end)')/Ll_eq),inf)
  norm(Check,inf)
  
  
  figure(2)
  clf()
  plot(L1,TT/Ep^2,'r'), hold on, 
  plot(L2,TT/Ep^2)
  axis([-LL,LL,0,TT(end)/Ep^2])
%   view(90,-90)
  ylabel('time','fontsize',18), xlabel('X','fontsize',18)
%    h = legend('\lambda_+','\lambda_-')
%     set(h,'fontsize',20) 
    title(['L = ',num2str(LL),' \tau = ',num2str(Tau)],'fontsize',20)
  Output = [path_figs,name,'DD',num2str(DD),'Tau',num2str(Tau)];
% % %   print(gcf,'-depsc','-opengl','-r300',[path_figs,name,'DD',num2str(DD),'Tau',num2str(Tau)],'-r300')
% % %   print(gcf,'-dpdf','-opengl','-r300',[path_figs,name,'DD',num2str(DD),'Tau',num2str(Tau)],'-r300')
%   saveas(gcf,Output,'eps')
%   saveas(gcf,Output,'pdf')
  
  
  function dvldt = rhs(t,y)
    w = y(1:N);
    l1 = y(N+1);
    l2 = y(N+2);
    data = ones(N,3); data(:,2) = -2; lap = spdiags(data,-1:1,N,N); 
    lap(1,2) = 2; lap(N,N-1) = 2; % Neumann BCs
    dvldt(1:N,1) = DD_0*lap*w/H^2 + 1 + (heavi(X-l1) - heavi(X-l2))*LL/Ll_eq;
    dvldt(N+1,1) = Mu*spline(X,w,y(N+1));
    dvldt(N+2,1) = - Mu*spline(X,w,y(N+2));
  end

  function y = heavi(x)
    y = (abs(x)>Delta).*(x>0) + ...
        (abs(x)<=Delta).*(0.5 + x/(2*Delta) + sin(pi*x/Delta)/(2*pi));
  end

end
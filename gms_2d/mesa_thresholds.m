function mesa_thresholds

figure(2), clf()
for D = [10, 8, 4]
Eigenmode = 0.0:0.1:10;

for I = 1:length(Eigenmode)

  Ep = 0.0025;
%   D = 8;
  DD = D*Ep;
  Z = exp(2*pi*i*1/2);
  Wplus = 3.295209;
  B0 = 0.211376;
  Kap = 2;
%   Wminus = (1 - sqrt(1-4*B0))/2;
%   Beta = simpson(35);   %
  Beta = 1.49882;
  Uplus = sqrt(B0/Kap)*Wplus;
%   Uminus = sqrt(B0/Kap)*(1-sqrt(1-4*B0))/2;  %wrong
  Uminus = 0;
  V0 = sqrt(B0/Kap);
  L = .5;
  M = Eigenmode(I);  %eigenmode

  A_gms = inline('-u + u.^2./(v.*(1+Kap*u.^2))','u','v','Kap');
  B_gms = inline('- v + u.^2','u','v');

  
%   [X,Fval,Flag] = fsolve(@f_minimizer_fzero,[V0,0.1,Uplus],optimset('TolFun',1e-9))
%   [X1,Fval1,Flag1] = fminsearch(@f_minimizer_fminsearch,[V0-0.1,0.1,Uplus-0.09],optimset('TolFun',1e-9))
%   Uplus = X(3);
%   Uminus = X(2)
%   V0 = X(1);
  Bplus = B_gms(Uplus,V0);
  Bminus = B_gms(Uminus,V0);
  Lsmall = Bminus*L/(Bminus - Bplus);
  K0 = Beta*V0^2/(Ep*quad(@Da_Dv_gms,Uminus,Uplus));
  
%   Y = eta;
  Y(1) = eta(Lsmall,Lsmall); 
  Y(2) = eta(-Lsmall,Lsmall);
  Aa = ((Bminus - Bplus)*Y(1) - Bplus*Lsmall/V0)/D;
  Bb = - (Bminus - Bplus)*Y(2)/D;
%   Lambda1(I) = (-Bplus*Lsmall + (Bplus-Bminus)*...
%     (eta(Lsmall,-Lsmall,M) - eta(Lsmall,Lsmall,M)))*Ep/(D*K0) - Ep^2*M^2
%   Lambda2(I) = (-Bplus*Lsmall + (Bplus-Bminus)*...
%     (-eta(Lsmall,-Lsmall,M) - eta(Lsmall,Lsmall,M)))*Ep/(D*K0) - Ep^2*M^2
  LL1(I) = (Aa + norm(Bb))/(K0) - Ep^2*M^2;
  LL2(I) = (Aa - norm(Bb))/(K0) - Ep^2*M^2;
  LL1(I) = ((Bminus - Bplus)*(Y(1)+Y(2)) - Bplus*Lsmall/V0)/(D*K0) - Ep^2*M^2;
  LL2(I) = ((Bminus - Bplus)*(Y(1)-Y(2)) - Bplus*Lsmall/V0)/(D*K0) - Ep^2*M^2;
%   LL1(I) = (-V0*Wplus^2*(Y(1)+Y(2)) - Bplus*Lsmall/V0)/(D*K0) - Ep^2*M^2;
%   LL2(I) = (-V0*Wplus^2*(Y(1)-Y(2)) - Bplus*Lsmall/V0)/(D*K0) - Ep^2*M^2;
%   LL1(I) = (Y(1)+Y(2));
%   LL2(I) = (Y(1)-Y(2));
end
% figure(1)
% plot(Eigenmode, Lambda1), hold on 
% plot(Eigenmode, Lambda2,'--')

figure(2)

plot(Eigenmode, LL1), hold on 
plot(Eigenmode, LL2,'--')
% axis([0 10 -0.0013 0.0000])
% axis([0 10 -0.00012 0.00004])
end

T1 = quad(@Da_Dv_gms,Uminus,Uplus)
T2 = -Wplus^2*V0/2
T3 = V0
T4 = - 2*Beta/(Ep*V0*Wplus^2)
T5 = (Bminus - Bplus)/(D*K0)
T6 = Ep*L*V0^2*Wplus^2/(2*Lsmall*D*Beta)

XX = 0:0.01:L;
for I=1:length(XX)
  YY(I) = eta(XX(I),Lsmall);
end
figure(3)
plot(XX,YY)
[Lsmall, L]
  



  function y = big_F(w,b)

    f_of_w = inline('- w + w.^2./(1+b.*w.^2)','w','b');
    y = -quad(f_of_w,0,w,1e-10,[],b);
  end

  function y = simpson(n)
    w = linspace(0,Wplus,n);
    h = Wplus/(n-1);
    for i=1:length(w)
      f(i) = sqrt(2*big_F(w(i),B0));
    end
    y = (f(1) + f(n) + 4*sum(f(2:2:n-1)) + 2*sum(f(3:2:n-1)))*h/3;
  end

  function y = eta(x,x0)
    b = exp(2*M)/(2*M*(Z-exp(2*M)));
    c = conj(b);
%     y(1) = c+1/(2*M) + b;   %eta(l;l)
%     y(2) = (c+1/(2*M))*exp(2*M*Lsmall) + b*exp(-2*M*Lsmall); % eta(-l,l)
    if x<x0
      y = (c + 1/(2*M))*exp(-M*(x-x0)) + b*exp(M*(x-x0));
    else
      y = (b + 1/(2*M))*exp(M*(x-x0)) + c*exp(-M*(x-x0));
    end
  end

  function y = Da_Dv_gms(u,v)
    v = V0;
    y = - u.^2./(v.^2.*(1+Kap*u.^2));
  end

  function z = f_minimizer_fzero(x)
    v0 = x(1); 
    uminus = x(2); 
    uplus = x(3);
    f_of_uv = inline('- u + u.^2./(v.*(1 + k*u.^2))','u','v','k');
    z(1) = quad(f_of_uv,uminus,uplus,[],[],v0,Kap);
    z(2) = f_of_uv(uminus,v0,Kap);
    z(3) = f_of_uv(uplus,v0,Kap);
%     y = norm(z);

%   f_of_a = inline('x','x');
%   a0 = x(1);
%   a1 = x(2);
%   z(1) = quad(f_of_a,a0,a1);
%   z(2) = a0.^2 + a1.^2 - 2;
%   y = norm(z);
  end

  function y = f_minimizer_fminsearch(x)
    v0 = x(1); 
    uminus = x(2); 
    uplus = x(3);
    f_of_uv = inline('- u + u.^2./(v.*(1 + k*u.^2))','u','v','k');
    z(1) = quad(f_of_uv,uminus,uplus,[],[],v0,Kap);
    z(2) = f_of_uv(uminus,v0,Kap);
    z(3) = f_of_uv(uplus,v0,Kap);
    y = norm(z);
  end

end
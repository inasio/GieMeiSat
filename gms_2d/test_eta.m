function test_eta

% solve eta'' - m^2*eta = delta(x-x0)

Z = exp(-i*34);
M = 2.3;
H = 0.01;
L = 1;
XX = H:H:L-H;
X0 = 0.5;
for I = 1:length(XX)
  X = XX(I);
  YY(I) = (eta(X+H,X0) - 2*eta(X,X0) + eta(X-H,X0))/H^2 - M^2*eta(X,X0);
end

eta(L,X0) - Z*eta(-L,X0)
(eta(L+H,X0) - eta(L,X0))/H - Z*(eta(-L+H,X0) - eta(-L,X0))/H
plot(XX,YY)


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

end
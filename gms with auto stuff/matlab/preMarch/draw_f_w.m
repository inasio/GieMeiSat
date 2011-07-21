b = 0.25/1.2
w = -.5:0.01:4
f = -w + w.^2./(1+b*w.^2);
plot(w,f)
grid on
xlabel('w','fontsize',20)
ylabel('f(w)','fontsize',20)
line([-.5 4], [0 0],'linewidth',2)
line([0 0], [-.5 1],'linewidth',2)
text(1.3, -.15, 'w_{-}','fontsize',18)
text(3.25, -.15, 'w_{+}','fontsize',18)
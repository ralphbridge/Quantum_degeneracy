v = VideoWriter('DISPERSION.avi');
open(v);
tp=10;
GVD=-1;

vg=2
kl=0.25;
wl=100;
for t=-10:0.5:600
    GVD=1;
 z=-500:0.01:1800;
A=tp*(sqrt(tp^2-1i*GVD*z)).^(-1).*exp(-(t-z/vg).^2./(2*(tp^2-1i*GVD*z)));
E=real(A.*exp(+1i*kl*z-1i*wl*t));
subplot(2,1,1)
graph1=plot(z,E);
title('Dispersion in Positive GVD Medium');
set(graph1,'LineWidth',1.5);
xlabel('Space (arb)');
ylabel('Amplitude (arb)');
axis([(-500+vg*t),(500+vg*t), -1, 1])


GVD=-1;
 z=-500:0.01:1800;
A=tp*(sqrt(tp^2-1i*GVD*z)).^(-1).*exp(-(t-z/vg).^2./(2*(tp^2-1i*GVD*z)));
E=real(A.*exp(+1i*kl*z-1i*wl*t));
subplot(2,1,2)
graph2=plot(z,E);
title('Dispersion in Negative GVD Medium');
set(graph2,'LineWidth',1.5);
xlabel('Space (arb)');
ylabel('Amplitude (arb)');
axis([(-500+vg*t),(500+vg*t), -1, 1])

frame = getframe(gcf);
writeVideo(v,frame);
pause(0.05)
end
close(v)


%%NORMAL DISPERSION 
v = VideoWriter('DISPERSION2.avi');
open(v);
tp=20;
GVD=-1;

vg=2
kl=0.5;
wl=100;
for t=-10:0.5:600
    
 z=-500:0.01:1800;
A=tp*(sqrt(tp^2-1i*GVD*z)).^(-1).*exp(-(t-z/vg).^2./(2*(tp^2-1i*GVD*z)));
E=real(A.*exp(+1i*kl*z-1i*wl*t));
graph1=plot(z,E);
title('Dispersion in Positive GVD Medium');
set(graph1,'LineWidth',1.5);
xlabel('Space (arb)');
ylabel('Amplitude (arb)');
axis([(-500+vg*t),(500+vg*t), -1, 1])

frame = getframe(gcf);
writeVideo(v,frame);
end
close(v)
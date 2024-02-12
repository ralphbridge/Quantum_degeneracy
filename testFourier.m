clear all
clc

Ts = 1/50;
t = 0:Ts:10-Ts;
%x = sin(2*pi*15*t) + sin(2*pi*20*t);
x=exp(-((t-5).^2)/(4*0.3^2));

% subplot(3,1,1)
% plot(t,x)
% xlabel('Time (seconds)')
% ylabel('Amplitude')

y = fft(x);
fs = 1/Ts;
f = (0:length(y)-1)*fs/length(y);

% subplot(3,1,2)
% plot(f,abs(y))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude')
% 
% n = length(x);
% fshift = (-n/2:n/2-1)*(fs/n);
% yshift = fftshift(y);
% subplot(3,1,3)
% plot(fshift,abs(yshift))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude')
% 
xf=ifft(y);
% figure
% subplot(3,1,1)
% plot(t,x)
% xlabel('Time (seconds)')
% ylabel('Amplitude')

yshift2=exp(-((f-5).^2)/(4*0.3^2));

subplot(2,1,1)
plot(f,yshift2)
xlabel('Frequency (Hz)')
ylabel('Magnitude')

xshift = ifftshift(yshift2);
subplot(2,1,2)
plot(t,abs(xshift))
xlabel('Time (seconds)')
ylabel('Magnitude')

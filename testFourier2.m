clear all
clc

n=1000;
f=linspace(1e13,10e14,n);
df=(max(f)-min(f))/n;

y=exp(-(abs(f-4e14)).^2/(4*(0.5e13^2)));
figure
subplot(2,1,1)
plot(f,y)
title('Test spectrum')
xlabel('Frequency (Hz)')
ylabel('Test spectrum')

t=(0:length(f)-1)/(n*df);
x=ifft(y);
subplot(2,1,2)
plot(t,abs(ifftshift(x)))
xlabel('Time (seconds)')
ylabel('Recovered x=ifft(y)')
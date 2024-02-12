clear all
clc

c=299792458;

S=readmatrix("spectrum.xlsx");
lambda=S(2:length(S),1)*1e-9;
I_lambda10fs=(S(2:length(S),2)-min(S(2:length(S),2)))/max(S(2:length(S),2));
I_lambda100fs_uv=(S(2:length(S),3)-min(S(2:length(S),3)))/max(S(2:length(S),3));
I_lambda100fs_ir=(S(2:length(S),4)-min(S(2:length(S),4)))/max(S(2:length(S),4));

I_f10fs=zeros(size(I_lambda10fs));
I_f100fs_uv=zeros(size(I_lambda100fs_uv));
I_f100fs_ir=zeros(size(I_lambda100fs_ir));

f=zeros(size(lambda));
for i=1:length(lambda)
    f(length(lambda)-i+1)=c/lambda(i);
    I_f10fs(length(lambda)-i+1)=lambda(i)^2*I_lambda10fs(i)/c;
    I_f100fs_uv(length(lambda)-i+1)=lambda(i)^2*I_lambda100fs_uv(i)/c;
    I_f100fs_ir(length(lambda)-i+1)=lambda(i)^2*I_lambda100fs_ir(i)/c;
end

t=linspace(1/max(f),1/min(f),length(f));

Iftest=exp(-(f-8e14).^2/(4*(4e13)^2));

spectrum=Iftest;
pulse=ifftshift(spectrum);

figure
subplot(2,1,1)
plot(f,spectrum)
xlabel('Frequency $Hz$','interpreter','latex')
ylabel('Amplitude')

subplot(2,1,2)
plot(t,pulse)
xlabel('Time $s$','interpreter','latex')
ylabel('Amplitude')

fmirror=zeros(2*length(f)+1,1);
Imirror=zeros(size(fmirror));

for i=1:length(f)
    fmirror(i)=-f(length(f)-i+1);
    fmirror(i+length(f)+1)=f(i);
    Imirror(i)=Iftest(length(f)-i+1);
    Imirror(i+length(f)+1)=Iftest(i);
end

figure
subplot(2,1,1)
plot(fmirror,Imirror)
xlabel('Frequency $Hz$','interpreter','latex')
ylabel('Amplitude')

subplot(2,1,2)
plot(t,abs(ifft(Imirror)))
xlabel('Time $s$','interpreter','latex')
ylabel('Amplitude')
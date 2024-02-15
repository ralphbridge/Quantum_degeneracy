clear all
clc

c=299792458;

S=readmatrix("spectrum.xlsx");

n=length(S)-1;

lambda=S(2:n+1,1)*1e-9;
Ilambda10fs=(S(1:n+1,2)-min(S(2:n+1,2)))/max(S(2:n+1,2));
Ilambda100fs_uv=(S(1:n+1,3)-min(S(2:n+1,3)))/max(S(2:n+1,3));
Ilambda100fs_ir=(S(2:n+1,4)-min(S(2:n+1,4)))/max(S(2:n+1,4));

f=zeros(n,1);
If10fs=zeros(n,1);
If100fs_uv=zeros(n,1);
If100fs_ir=zeros(n,1);

for i=1:n
    f(n-i+1)=c/lambda(i);
    If10fs(n-i+1)=lambda(i)^2*Ilambda10fs(i)/c;
    If100fs_uv(n-i+1)=lambda(i)^2*Ilambda100fs_uv(i)/c;
    If100fs_ir(n-i+1)=lambda(i)^2*Ilambda100fs_ir(i)/c;
end

df=(max(f)-min(f))/n;
t=(0:n-1)/(n*df);

Iftest=exp(-(f-8e14).^2/(4*(0.5e13)^2));

spectrum=If10fs;

for i=1:n
    if spectrum(i)<=0.01*max(spectrum)
        spectrum(i)=0;
    end
end

pulse=ifft(spectrum);

figure
subplot(2,1,1)
plot(f,spectrum)
xlabel('Frequency $Hz$','interpreter','latex')
ylabel('Amplitude')

It=abs(ifftshift(pulse));
subplot(2,1,2)
plot(t,It)
xlabel('Time $s$','interpreter','latex')
ylabel('Amplitude')

Iprofile=zeros(size(It));
tprofile=zeros(size(t));
j=1;
for i=1:n
    if It(i)>=0.005*max(It)
        Iprofile(j)=It(i);
        tprofile(j)=t(i);
        j=j+1;
    end
end

Iprofile=Iprofile(find(Iprofile,1,'first'):find(Iprofile,1,'last'));
tprofile=tprofile(find(tprofile,1,'first'):find(tprofile,1,'last'));
tprofile=tprofile-min(tprofile);

figure
plot(tprofile,Iprofile)
clear all
clc

c=299792458;

S=readmatrix("spectrum.xlsx");

n=length(S)-1;

lambda=S(2:n+1,1);
Ilambda10fs=(S(2:n+1,2)-min(S(2:n+1,2)))/max(S(2:n+1,2));
Ilambda100fs_uv=(S(2:n+1,3)-min(S(2:n+1,3)))/max(S(2:n+1,3));
Ilambda100fs_ir=(S(2:n+1,4)-min(S(2:n+1,4)))/max(S(2:n+1,4));

% figure
% plot(lambda,Ilambda10fs,'LineWidth',2)
% xlabel('Wavelength $\lambda\ nm$','interpreter','latex','FontSize',15)
% ylabel('Intensity as a function of wavelength','FontSize',15)
% grid on
% 
% hold on
% plot(lambda,0.75*exp(-((lambda-755)/(2*10)).^2),'LineWidth',2)
% plot(lambda,0.45*exp(-((lambda-815)/(2*20)).^2),'LineWidth',2)

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

spectruml=Ilambda10fs;
spectrum=If10fs;

% for i=1:n % Used to clean up the spectrum
%     if spectrum(i)<=0.01*max(spectrum)
%         spectrum(i)=0;
%     end
% end

pulse=ifft(spectrum);

figure
plot(lambda,spectruml,'linewidth',2)
xlabel('Wavelength $(nm)$','interpreter','latex','fontsize',20)
ylabel('Power spectral density (a.u.)','interpreter','latex',fontsize=20)
axis([600 1000 0 1])
grid on

It=abs(ifftshift(pulse));

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

% figure
% plot(tprofile,Iprofile)

FSAC=readmatrix("FSAC.xlsx");

m=length(FSAC)-1;

clear t

t=FSAC(2:m+1,1)-8.77e-5;
trace=FSAC(2:m+1,2);

figure
plot(t,trace,linewidth=2)
xlabel('Delay (a.u)','interpreter','latex','fontsize',20)
ylabel('Autocorrelation signal (a.u.)','interpreter','latex',fontsize=20)
axis([-1e-3 1e-3 0 1.1*max(trace)])
grid on
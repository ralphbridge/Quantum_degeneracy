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

% subplot(2,1,1)
% plot(lambda,I_lambda10fs)
% grid on
% 
% subplot(2,1,2)
% plot(f,I_f10fs)
% grid on
% 
% figure
% subplot(2,1,1)
% plot(lambda,I_lambda100fs_uv)
% grid on
% 
% subplot(2,1,2)
% plot(f,I_f100fs_uv)
% grid on
% 
% figure
% subplot(2,1,1)
% plot(lambda,I_lambda100fs_ir)
% grid on
% 
% subplot(2,1,2)
% plot(f,I_f100fs_ir)
% grid on

%w=2*pi*f;
% I_w10fs=zeros(size(I_f10fs));
% I_w100fs_uv=zeros(size(I_f100fs_uv));
% I_w100fs_ir=zeros(size(I_f100fs_ir));
% for i=1:length(w)
%     I_w10fs(i)=I_f10fs(i)/(2*pi);
%     I_w100fs_uv(i)=I_f100fs_uv(i)/(2*pi);
%     I_w100fs_ir(i)=I_f100fs_ir(i)/(2*pi);
% end

% subplot(2,1,1)
% plot(lambda,I_lambda10fs)
% grid on
% 
% subplot(2,1,2)
% plot(w,I_w10fs)
% grid on
% 
% figure
% subplot(2,1,1)
% plot(lambda,I_lambda100fs_uv)
% grid on
% 
% subplot(2,1,2)
% plot(w,I_w100fs_uv)
% grid on
% 
% figure
% subplot(2,1,1)
% plot(lambda,I_lambda100fs_ir)
% grid on
% 
% subplot(2,1,2)
% plot(w,I_w100fs_ir)
% grid on

% I_t10fs=real(ifft(I_w10fs));
% 
% figure
% subplot(2,1,1)
% plot(w,I_w10fs)
% xlabel('$\omega$ rad/s','Interpreter','latex')
% grid on
% subplot(2,1,2)
% plot(t,I_t10fs)
% xlabel('$t$ s')
% grid on
% title('10 fs laser')
% 
% I_t100fs_uv=real(ifft(I_w100fs_uv));
% 
% figure
% subplot(2,1,1)
% plot(w,I_w100fs_uv)
% xlabel('$\omega$ rad/s','Interpreter','latex')
% grid on
% subplot(2,1,2)
% plot(t,I_t100fs_uv)
% xlabel('$t$ s')
% grid on
% 
% I_t100fs_ir=real(ifft(I_w100fs_ir));
% 
% figure
% subplot(2,1,1)
% plot(w,I_w100fs_ir)
% xlabel('$\omega$ rad/s','Interpreter','latex')
% grid on
% subplot(2,1,2)
% plot(t,I_t100fs_ir)
% xlabel('$t$ s')
% grid on

w=linspace(-10,10,1001);

t=zeros(size(w));
for i=1:length(w)
    t(length(w)-i+1)=w(i);
end

figure
%I_w=exp(-(w).^2/(4*(2)^2));
I_w=zeros(size(w));
for i=1:length(w)
    if w(i)==5
        I_w(i)=1;
    end
end
subplot(3,1,1)
plot(w,I_w)

subplot(3,1,2)
plot(t,fft(I_w),'o')

subplot(3,1,3)
plot(t,abs(fft(I_w)),'o')

X = [0 1 3 1 0];
Y = fft(X);

figure
subplot(2,1,1)
plot(X)

subplot(2,1,2)
plot(abs(Y).^2)
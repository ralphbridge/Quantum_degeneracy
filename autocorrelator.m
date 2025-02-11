clear all
clc

eps0=8.85e-12;
c=299792458;
lambda0=756e-9; % laser central wavelength
w0=2*pi*c/lambda0;
E00=sqrt(2*1e14/(c*eps0));

GDD_laser=124e-30; % Check FSAC manual for correct values
TOD_laser=0;
GDD_fsac=230e-30;
TOD_fsac=345e-45;

n_bounces=14;

zd=2.06+0.05*(n_bounces-1); % Distance from laser pinhole to FSAC pinhole
zw=0e-3; % Window thickness (measure this again)

syms x

n_air=1+0.05792105/(238.0185-(1e-12)*(x/(2*pi*c))^2)+0.00167917/(57.362-(1e-12)*(x/(2*pi*c))^2);
n_bk7=(1+1.03961212*(1e12)*(2*pi*c/x)^2/((1e12)*(2*pi*c/x)^2-0.00600069867)+0.231792344*(1e12)*(2*pi*c/x)^2/((1e12)*(2*pi*c/x)^2-0.0200179144)+1.01046945*(1e12)*(2*pi*c/x)^2/((1e12)*(2*pi*c/x)^2-103.560653));

n0_air=double(subs(n_air,w0));
np0_air=double(subs(diff(n_air,x),w0));
npp0_air=double(subs(diff(diff(n_air,x),x),w0));
nppp0_air=double(subs(diff(diff(diff(n_air,x),x)),w0));

k0_air=n0_air*w0/c;
kp0_air=(n0_air+w0*np0_air)/c;
kpp0_air=(2*np0_air+w0*npp0_air)/c;
kppp0_air=(3*npp0_air+w0*nppp0_air)/c;

n0_bk7=double(subs(n_bk7,w0));
np0_bk7=double(subs(diff(n_bk7,x),w0));
npp0_bk7=double(subs(diff(diff(n_bk7,x),x),w0));
nppp0_bk7=double(subs(diff(diff(diff(n_bk7,x),x)),w0));

k0_bk7=n0_bk7*w0/c;
kp0_bk7=(n0_bk7+w0*np0_bk7)/c;
kpp0_bk7=(2*np0_bk7+w0*npp0_bk7)/c;
kppp0_bk7=(3*npp0_bk7+w0*nppp0_bk7)/c;

function res=InverseFourier(Fw,w,t)
    ft=zeros(size(w,1),1);
    for j=1:length(t)
        for i=1:length(w)
            ft(j)=ft(j)+Fw(i)*exp(1j*w(i)*t(j))+conj(Fw(i))*exp(-1j*w(i)*t(j));
        end
    end
    res=real(ft);
end

function res=E0shift(E0,j)
    res=zeros(length(E0));
    Eftemp=append(append(E0,E0),E0);
    for i=1:length(E0)
        idx=int((length(E0)-1)/2)+i+j;
        res(i)=Eftemp(idx);
    end
    clear Eftemp
end

function S=S_l(t,E0)
    S=zeros(length(t));
    j=0;
    for i=1:length(t)            
        S(j)=trapz((E0+E0shift(E0,j))^2,t(i));
        j=j+1;
    end
end

function S=S_q(t,E0)
    S=np.zeros(len(t));
    j=0;
    for i=1:length(t)            
        S(j)=trapz((E0+E0shift(E0,j))^4,t(i));
        j=j+1;
    end
end

function res=roundup(x)
    res=int32(ceil(x/100))*100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting measured spectra

S=readmatrix("spectrum.xlsx");
n=size(S,1);
lam=zeros(n,1);
Il=zeros(n,1);

for j=1:size(S,2)
    for i=1:size(S,1)
        %S(i,j)
        if j==1
            lam(i)=(S(i,j))*1e-9;
        elseif j==2
            Il(i)=S(i,j);
        end
    end
end

Il=(Il-min(Il))/max(Il); % Initial measured spectrum

f=zeros(n,1);
spectrumf=zeros(n,1);

for i=1:n
    f(n-i+1)=c/(lam(i));
    spectrumf(n-i+1)=(lam(i)^2)*Il(i)/c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Increasing time resolution (by increasing frequency range)

Ef0=sqrt(2*spectrumf/(c*eps0));
t0=linspace(0,50e-15,length(f));
Et0=InverseFourier(Ef0,2*pi*f,t0);

%%%%%%%%%%%% Trimming spectrum low and high frequencies and interpolating it to get equally spaced frequencies

NN=799; % Length of actual spectrum (obtained from matlab)
f_trim=zeros(NN,1);
spectrumf_trim=zeros(NN,1);

for i=1:NN
    f_trim(i)=f(i+266); % 266 is the position of the original spectrum where it starts being bigger than 0
    spectrumf_trim(i)=spectrumf(i+266); % New angular frequency spectrum (trimmed)
end

%spectrumf_interp_function=csapi(f_trim,spectrumf_trim);

N=3000;
f_interp=linspace(min(f_trim),max(f_trim),N)';
spectrumf_interp=zeros(N,1);
df=f_interp(2)-f_interp(1);

w=2*pi*f_interp;

for i=1:N
    if i==1
        spectrumf_interp(i)=spectrumf_trim(i);
    elseif i==N
        spectrumf_interp(i)=spectrumf_trim(NN);
    else
        spectrumf_interp(i)=csapi(f_trim,spectrumf_trim,f_interp(i));%spectrumf_interp_function(f_interp(i));
        if spectrumf_interp(i)<0
            spectrumf_interp(i)=0;
        end
    end
end

Ef0_interp=sqrt(2*spectrumf_interp/(c*eps0))*(1+0j); % Initial electric field spectrum (trimmed)

% Initial phases for original measured spectrum (up to third order)

for i=1:N
    % Dispersion phases introduced by the laser itself to the pulse
    Ef0_interp(i)=Ef0_interp(i)*exp(1j*GDD_laser*(w(i)-w0)^2/factorial(2));
    Ef0_interp(i)=Ef0_interp(i)*exp(1j*TOD_laser*(w(i)-w0)^3/factorial(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t0_interp=linspace(0,max(t0),length(f_interp));%linspace(0,50e-15,length(f_interp));
Et0_interp=InverseFourier(Ef0_interp,2*pi*f_interp,t0_interp);

% fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
% plt.subplot(2,1,1)
% line1,=plt.plot(t0*1e15,Et0)
% ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
% ax1.set_ylabel(r'$E_0\ V/m$', fontsize=16)
% ax1.grid(which='both')

% plt.subplot(2,1,2)
% line2,=ax2.plot((t0_interp)*1e15,Et0_interp)
% ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
% ax2.set_ylabel('$E_0^{interp}\ V/m$', fontsize=16)
% ax2.grid(which='both')
% plt.show()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting GDD data from Thorlabs chirped mirrors data

CM=readmatrix("UMxx-15FS_data.xlsx");
n_cm=size(CM,1);
lam_cm=zeros(n_cm,1);
GD_cm_lam=zeros(n_cm,1);

for j=1:size(CM,2)
    for i=1:n_cm
        if j==0
            lam_cm(i)=(CM(i,j))*1e-9;
        else
            GD_cm_lam(i)=(CM(i,j))*1e-15;
        end
    end
end

w_cm=zeros(n_cm,1);
GD_cm=zeros(n_cm,1);
for i=1:n_cm
    w_cm(n_cm-1-i)=2*pi*c/lam_cm(i);
    GD_cm(n_cm-1-i)=GD_cm_lam(i);
end

GD_cm_interp_function=csapi(w_cm,GD_cm);
GD_cm_interp=GD_cm_interp_function(w);
% GD_cm_interp=np.interp(w,w_cm,GD_cm)

GDD_cm=zeros(N,1);
TOD_cm=zeros(N,1);
for i=1:N % Higher accuracy order derivatives from Fornberg 1988
    if i==0
        m2=(GD_cm_interp(i+1)-GD_cm_interp(i))/(w(i+1)-w(i));
        GDD_cm(i)=m2;
    elseif i==N-1
        m1=(GD_cm_interp(i)-GD_cm_interp(i-1))/(w(i)-w(i-1));
        GDD_cm(i)=m1;
    elseif i==1||i==N-2
        GDD_cm(i)=(-GD_cm_interp(i-1)+GD_cm_interp(i+1))/(2*2*pi*df);
        TOD_cm(i)=(GD_cm_interp(i-1)-2*GD_cm_interp(i)+GD_cm_interp(i+1))/((2*pi*df)^2);
    else
        GDD_cm(i)=(GD_cm_interp(i-2)-8*GD_cm_interp(i-1)+8*GD_cm_interp(i+1)-GD_cm_interp(i+2))/(12*2*pi*df);
        TOD_cm(i)=(-GD_cm_interp(i-2)+16*GD_cm_interp(i-1)-30*GD_cm_interp(i)+16*GD_cm_interp(i+1)-GD_cm_interp(i+2))/(12*(2*pi*df)^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting GDD data from Thorlabs P01 mirrors data

P01=readmatrix("P01_data.xlsx");
n_p01=size(P01,1);
lam_p01=zeros(n_p01,1);
GDD_p01_data_lam=zeros(n_p01,1);

for j=1:size(P01,2)
    for i=1:n_p01
        if j==0
            lam_p01(i)=(P01(i,j))*1e-9;
        else
            GDD_p01_data_lam(i)=(P01(i,j))*1e-30;
        end
    end
end

w_p01=zeros(n_p01,1);
GDD_p01_data=zeros(n_p01,1);
for i=1:n_p01
    w_p01(n_p01-1-i)=2*pi*c/lam_p01(i);
    GDD_p01_data(n_p01-1-i)=GDD_p01_data_lam(i);
end

GDD_p01_interp_function=csapi(w_p01,GDD_p01_data);
GDD_p01_interp=GDD_p01_interp_function(w);

TOD_p01_data=zeros(n_p01,1);
for i=1:n_p01
    if i==0
        m2=(GDD_p01_data(i+1)-GDD_p01_data(i))/(w_p01(i+1)-w_p01(i));
        TOD_p01_data(i)=m2;
    elseif i==n_p01-1
        m1=(GDD_p01_data(i)-GDD_p01_data(i-1))/(w_p01(i)-w_p01(i-1));
        TOD_p01_data(i)=m1;
    else
        m1=(GDD_p01_data(i)-GDD_p01_data(i-1))/(w_p01(i)-w_p01(i-1));
        m2=(GDD_p01_data(i+1)-GDD_p01_data(i))/(w_p01(i+1)-w_p01(i));
        TOD_p01_data(i)=(m1+m2)/2;
    end
end

TOD_p01_interp_function=csapi(w_p01,GDD_p01_data);
TOD_p01_interp=TOD_p01_interp_function(w);

%print(GDD_p01_interp(2*pi*c/800e-9)*1e30) # Check this line, it does not interpolate correctly

%%%%%%%%%%%%%%%%%%%%%%%% Computing the phases due to each element in the layout

Ef=zeros(size(Ef0_interp,1));

for i=1:N
    % Dispersion phases introduced by air up to the third order
    % Ef(i)=Ef0_interp(i)*exp(1j*k0_air*zd);
    % Ef(i)=Ef(i)*exp(1j*kp0_air*(w(i)-w0)*zd);
    Ef(i)=Ef0_interp(i)*exp(1j*kpp0_air*(w(i)-w0)^2*zd/factorial(2));
    % Ef(i)=Ef(i)*exp(1j*kppp0_air*(w(i)-w0)^3*zd/factorial(3));
    
    % Dispersion phases introduced by BK7 Fused Silica glass window up to the third order
    % Ef(i)=Ef(i)*exp(1j*k0_bk7*zw);
    % Ef(i)=Ef(i)*exp(1j*kp0_bk7*(w(i)-w0)*zw);
    % Ef(i)=Ef(i)*exp(1j*kpp0_bk7*(w(i)-w0)^2*zw/factorial(2));
    % Ef(i)=Ef(i)*exp(1j*kppp0_bk7*(w(i)-w0)^3*zw/factorial(3));
    
    % Dispersion phases introduced by bounces off Chirped mirrors
    % Ef(i)=Ef(i)*exp(1j*GDD_cm(i)*(w(i)-w0)^2/factorial(2));
    % Ef(i)=Ef(i)*exp(1j*TOD_cm(i)*(w(i)-w0)^3/factorial(3));
    Ef(i)=Ef(i)*exp(1j*n_bounces*GD_cm_interp(i)*(w(i)-w0));
    
    % % Dispersion phases introduced by bounces off P01 Silver mirrors
    % Ef(i)=Ef(i)*exp(1j*GDD_p01(i)*(w(i)-w0)^2/factorial(2));
    % Ef(i)=Ef(i)*exp(1j*TOD_p01(i)*(w(i)-w0)^3/factorial(3));
    Ef(i)=Ef(i)*exp(1j*GDD_p01_interp(i)*(w(i)-w0)^2/factorial(2));
    
    % % Dispersion phases introduced by the FSAC up to the third order
    Ef(i)=Ef(i)*exp(1j*GDD_fsac*(w(i)-w0)^2/factorial(2));
    Ef(i)=Ef(i)*exp(1j*TOD_fsac*(w(i)-w0)^3/factorial(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fs=max(f_interp)-min(f_interp) % Sampling frequency (equals the bandwidth or maximum frequency)
t=linspace(-3e-12,3e-12,len(f_interp));
Et=InverseFourier(Ef,2*pi*f_interp,t);

% fig,(ax1,ax2)=subplots(2,1,tight_layout=True)
% subplot(2,1,1)
% line1,=plot(t0*1e15,Et0)
% ax1.set_xlim([-20,20])
% ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
% ax1.set_ylabel(r'$E_0\ V/m$', fontsize=16)
% % major_tick = np.arange(roundup(min(lam)*1e9), roundup(max(lam)*1e9),200)#[200, 400, 600, 800, 1000]
% % minor_tick = np.arange(roundup(min(lam)*1e9)+100, roundup(max(lam)*1e9),200)#[300, 500, 700, 900]
% % ax1.set_xticks(major_tick) # Grid
% % ax1.set_xticks(minor_tick, minor=True)
% ax1.grid(which='both')

% plt.subplot(2,1,2)
% line2,=ax2.plot(t*1e15,Et)
% ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
% ax2.set_ylabel('$E\ V/m$', fontsize=16)
% % ax2.set_xlim([-100,100])
% % major_tick = np.arange(-100,100,20)
% % minor_tick = np.arange(-100,100,10)
% % ax2.set_xticks(major_tick)
% % ax2.set_xticks(minor_tick, minor=True)
% ax2.grid(which='both')

% % fig.savefig("Frequency_Time.pdf",bbox_inches='tight')
% plt.show()

%%%%%%%%%%%%%%

Slinear=S_l(t,Et);
Squad=S_q(t,Et);

% subplot(3,1,1)

% line1,=ax1.plot(t*1e15,Et,lw=1)
% ax1.plot(t*1e15,Et,'.',alpha=0.01)
% ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
% ax1.set_ylabel(r'$E\ V/m$', fontsize=16)
% ax1.set_xlim([-50,50])
% ax1.grid()

% line2,=ax2.plot(t*1e15 ,Slinear,lw=1)
% ax2.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
% ax2.set_ylabel(r'$S_{linear}\ W/m^2$', fontsize=16)
% ax2.set_xlim([-20,20])
% ax2.grid()

plot((t-t(argmax(Squad)))*1e15,Squad)
% ax3.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
% ax3.set_ylabel(r'$S_{quadratic}\ W/m^2$', fontsize=16)
% ax3.set_xlim([-250,250])
grid on
% plt.savefig('FieldTrace_'+str(n_bounces)+'b_far.pdf',bbox_inches='tight')

% Check first, second and third order phases added due to air and BK7

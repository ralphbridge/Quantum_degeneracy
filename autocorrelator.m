% Converted from Python -> MATLAB
clear; clc;

%% Constants
eps0 = 8.85e-12;
c = 299792458;
lambda0 = 756e-9; % laser central wavelength
w0 = 2*pi*c/lambda0;
E00 = sqrt(2*1e14/(c*eps0));

GDD_laser = 124e-30; % Check FSAC manual for correct values
TOD_laser = 0;
GDD_fsac = 230e-30;
TOD_fsac = 345e-45;

n_bounces = 8;

zd = 2.06 + 0.05*(n_bounces-1); % Distance from laser pinhole to FSAC pinhole
zw = 0e-3; % Measure this again

%% Symbolic refractive index derivatives (follow Python sympy use)
syms x real

n_air = 1 + 0.05792105./(238.0185 - (1e-12)*(x/(2*pi*c)).^2) + 0.00167917./(57.362 - (1e-12)*(x/(2*pi*c)).^2);
n_bk7 = sqrt( 1 ...
    + 1.03961212*(1e12)*(2*pi*c/x).^2 ./((1e12)*(2*pi*c/x).^2 - 0.00600069867) ...
    + 0.231792344*(1e12)*(2*pi*c/x).^2 ./((1e12)*(2*pi*c/x).^2 - 0.0200179144) ...
    + 1.01046945*(1e12)*(2*pi*c/x).^2 ./((1e12)*(2*pi*c/x).^2 - 103.560653) );

% evaluate function and derivatives at w0
n0_air = double(subs(n_air, x, w0));
np0_air = double(subs(diff(n_air, x), x, w0));
npp0_air = double(subs(diff(n_air, x, 2), x, w0));
nppp0_air = double(subs(diff(n_air, x, 3), x, w0));

k0_air = double(n0_air * w0 / c);
kp0_air = double((n0_air + w0*np0_air)/c);
kpp0_air = double((2*np0_air + w0*npp0_air)/c);
kppp0_air = double((3*nppp0_air + w0*nppp0_air)/c); % matches Python structure though last term has duplication

n0_bk7 = double(subs(n_bk7, x, w0));
np0_bk7 = double(subs(diff(n_bk7, x), x, w0));
npp0_bk7 = double(subs(diff(n_bk7, x, 2), x, w0));
nppp0_bk7 = double(subs(diff(n_bk7, x, 3), x, w0));

k0_bk7 = double(n0_bk7 * w0 / c);
kp0_bk7 = double((n0_bk7 + w0*np0_bk7)/c);
kpp0_bk7 = double((2*np0_bk7 + w0*npp0_bk7)/c);
kppp0_bk7 = double((3*nppp0_bk7 + w0*nppp0_bk7)/c);

%% -- Helper functions (as nested functions) --
% Note: MATLAB nested functions capture workspace variables. If you prefer,
% move them to separate local functions in the same file.

    function ft = InverseFourier(Fw, w, t)
        % Implements: ft[j] = sum_i Fw[i]*exp(1j*w[i]*t[j]) + conj(Fw[i])*exp(-1j*w[i]*t[j])
        Nt = numel(t);
        Nw = numel(w);
        ft = zeros(Nt,1);
        for j = 1:Nt
            val = 0 + 0i;
            for i = 1:Nw
                val = val + Fw(i)*exp(1i*w(i)*t(j)) + conj(Fw(i))*exp(-1i*w(i)*t(j));
            end
            ft(j) = real(val);
        end
    end

    function Ef = E0shift(E0, jshift)
        % replicate behavior of Python function E0shift
        N = numel(E0);
        Ef = zeros(N,1);
        Eftemp = [E0; E0; E0]; % vertical concat
        center_idx = floor((N-1)/2) + 1; % convert to 1-based
        for ii = 1:N
            idx = center_idx + ii - 1 + jshift;
            Ef(ii) = Eftemp(idx);
        end
    end

    function S = S_l(tvec, E0)
        S = zeros(numel(tvec),1);
        for jj = 1:numel(tvec)
            Esum = E0 + E0shift(E0, jj-1);
            S(jj) = trapz(tvec, (Esum).^2 );
        end
    end

    function S = S_q(tvec, E0)
        S = zeros(numel(tvec),1);
        for jj = 1:numel(tvec)
            Esum = E0 + E0shift(E0, jj-1);
            S(jj) = trapz(tvec, (Esum).^4 );
        end
    end

    function out = roundup(x)
        out = ceil(x/100) * 100;
    end

%% #################### Getting measured spectra
Smat = readmatrix('spectrum.xlsx'); % assumes two columns: wavelength (nm), intensity
% If header exists, readmatrix should skip automatically; adjust if needed.

% ensure numeric
Smat = double(Smat);
[nrows, ncols] = size(Smat);
lam = zeros(nrows,1);
Il = zeros(nrows,1);

for j = 1:ncols
    for i = 1:nrows
        if j == 1
            lam(i) = Smat(i,j) * 1e-9; % convert nm -> m
        elseif j == 2
            Il(i) = Smat(i,j);
        end
    end
end

Il = (Il - min(Il)) / max(Il); % Initial measured spectrum (normalized)

% build frequency array reversed to match Python ordering
n = nrows;
f = zeros(n,1);
spectrumf = zeros(n,1);
for i = 1:n
    f(n - i + 1) = c / lam(i);
    spectrumf(n - i + 1) = (lam(i)^2) * Il(i) / c;
end

%% Increasing time resolution (by increasing frequency range)
Ef0 = sqrt(2 * spectrumf ./ (c * eps0));
t0 = linspace(0,50e-15,numel(f));
Et0 = InverseFourier(Ef0, 2*pi*f, t0);

%% Trimming spectrum low and high frequencies and interpolating it to get equally spaced frequencies
NN = 799; % Length of actual spectrum (obtained from matlab)
f_trim = zeros(NN,1);
spectrumf_trim = zeros(NN,1);

startIndex = 267; % Python used i+266 since zero-based; matlab is 1-based -> 266+1=267
for i = 1:NN
    f_trim(i) = f(i + startIndex - 1);
    spectrumf_trim(i) = spectrumf(i + startIndex - 1);
end

% cubic spline interpolation function
spectrumf_interp_function = @(xq) ppval(spline(f_trim, spectrumf_trim), xq);

N = 3500;
f_interp = linspace(min(f_trim), max(f_trim), N)';
spectrumf_interp = zeros(N,1);
df = f_interp(2) - f_interp(1);
w = 2*pi*f_interp;

for i = 1:N
    if i == 1
        spectrumf_interp(i) = spectrumf_trim(1);
    elseif i == N
        spectrumf_interp(i) = spectrumf_trim(end);
    else
        val = spectrumf_interp_function(f_interp(i));
        if val < 0
            val = 0;
        end
        spectrumf_interp(i) = val;
    end
end

Ef0_interp = sqrt(2 * spectrumf_interp ./ (c * eps0)) .* (1 + 0i); % complex

% Initial phases for original measured spectrum (up to third order)
for i = 1:N
    Ef0_interp(i) = Ef0_interp(i) * exp(1i * GDD_laser * (w(i)-w0)^2 / factorial(2));
    Ef0_interp(i) = Ef0_interp(i) * exp(1i * TOD_laser * (w(i)-w0)^3 / factorial(3));
end

t0_interp = linspace(0,50e-15,numel(f_interp));
Et0_interp = InverseFourier(Ef0_interp, 2*pi*f_interp, t0_interp);

%% Getting GDD data from Thorlabs chirped mirrors data
CM = readmatrix('UMxx-15FS_data.xlsx');
CM = double(CM);
[n_cm, m_cm] = size(CM);
lam_cm = zeros(n_cm,1);
GD_cm_lam = zeros(n_cm,1);

for j = 1:m_cm
    for i = 1:n_cm
        if j == 1
            lam_cm(i) = CM(i,j) * 1e-9;
        else
            GD_cm_lam(i) = CM(i,j) * 1e-15;
        end
    end
end

w_cm = zeros(n_cm,1);
GD_cm = zeros(n_cm,1);
for i = 1:n_cm
    w_cm(n_cm - i + 1) = 2*pi*c / lam_cm(i);
    GD_cm(n_cm - i + 1) = GD_cm_lam(i);
end

% cubic spline interpolation over angular frequency
GD_cm_interp = ppval(spline(w_cm, GD_cm), w);

GDD_cm = zeros(N,1);
TOD_cm = zeros(N,1);

% Finite-difference formulas (ported from Python)
for i = 1:N
    if i == 1
        m2 = (GD_cm_interp(i+1) - GD_cm_interp(i)) / (w(i+1) - w(i));
        GDD_cm(i) = m2;
    elseif i == N
        m1 = (GD_cm_interp(i) - GD_cm_interp(i-1)) / (w(i) - w(i-1));
        GDD_cm(i) = m1;
    elseif i == 2 || i == N-1
        GDD_cm(i) = (-GD_cm_interp(i-1) + GD_cm_interp(i+1)) / (2*2*pi*df);
        TOD_cm(i) = (GD_cm_interp(i-1) - 2*GD_cm_interp(i) + GD_cm_interp(i+1)) / ((2*pi*df)^2);
    else
        GDD_cm(i) = (GD_cm_interp(i-2) - 8*GD_cm_interp(i-1) + 8*GD_cm_interp(i+1) - GD_cm_interp(i+2)) / (12 * 2*pi*df);
        TOD_cm(i) = (-GD_cm_interp(i-2) + 16*GD_cm_interp(i-1) - 30*GD_cm_interp(i) + 16*GD_cm_interp(i+1) - GD_cm_interp(i+2)) / (12 * (2*pi*df)^2);
    end
end

%% Getting GDD data from Thorlabs P01 mirrors data
P01 = readmatrix('P01_data.xlsx');
P01 = double(P01);
[n_p01, m_p01] = size(P01);
lam_p01 = zeros(n_p01,1);
GDD_p01_data_lam = zeros(n_p01,1);

for j = 1:m_p01
    for i = 1:n_p01
        if j == 1
            lam_p01(i) = P01(i,j) * 1e-9;
        else
            GDD_p01_data_lam(i) = P01(i,j) * 1e-30;
        end
    end
end

w_p01 = zeros(n_p01,1);
GDD_p01_data = zeros(n_p01,1);
for i = 1:n_p01
    w_p01(n_p01 - i + 1) = 2*pi*c / lam_p01(i);
    GDD_p01_data(n_p01 - i + 1) = GDD_p01_data_lam(i);
end

GDD_p01_interp = ppval(spline(w_p01, GDD_p01_data), w);

% approximate TOD_p01_data (slope averaging)
TOD_p01_data = zeros(n_p01,1);
for i = 1:n_p01
    if i == 1
        m2 = (GDD_p01_data(i+1) - GDD_p01_data(i)) / (w_p01(i+1) - w_p01(i));
        TOD_p01_data(i) = m2;
    elseif i == n_p01
        m1 = (GDD_p01_data(i) - GDD_p01_data(i-1)) / (w_p01(i) - w_p01(i-1));
        TOD_p01_data(i) = m1;
    else
        m1 = (GDD_p01_data(i) - GDD_p01_data(i-1)) / (w_p01(i) - w_p01(i-1));
        m2 = (GDD_p01_data(i+1) - GDD_p01_data(i)) / (w_p01(i+1) - w_p01(i));
        TOD_p01_data(i) = (m1 + m2) / 2;
    end
end

TOD_p01_interp = ppval(spline(w_p01, GDD_p01_data), w); % as in Python (note: this uses same GDD data)

%% Computing the phases due to each element in the layout
Ef = zeros(numel(Ef0_interp),1,'like',1+1i);

for i = 1:N
    % Start with initial interpolated spectrum
    Ef(i) = Ef0_interp(i);
    % Dispersion phases introduced by air up to the second order (as in your uncommented code)
    Ef(i) = Ef(i) * exp(1i * kpp0_air * (w(i)-w0)^2 * zd / factorial(2));
    
    % Dispersion phases introduced by BK7 Fused Silica window (commented out in Python)
    % Ef(i) = Ef(i) * exp(1i * kpp0_bk7 * (w(i)-w0)^2 * zw / factorial(2));
    
    % Dispersion phases introduced by n_bounces bounces off Chirped mirrors: applied as n_bounces*GD_cm_interp*(w-w0)
    Ef(i) = Ef(i) * exp(1i * n_bounces * GD_cm_interp(i) * (w(i)-w0));
    
    % Dispersion from P01
    Ef(i) = Ef(i) * exp(1i * GDD_p01_interp(i) * (w(i)-w0)^2 / factorial(2));
    
    % FSAC terms
    Ef(i) = Ef(i) * exp(1i * GDD_fsac * (w(i)-w0)^2 / factorial(2));
    Ef(i) = Ef(i) * exp(1i * TOD_fsac * (w(i)-w0)^3 / factorial(3));
end

%% Inverse Fourier to time domain
t = linspace(-3e-12, 3e-12, numel(f_interp));
Et = InverseFourier(Ef, 2*pi*f_interp, t);

%% Compute quadratic S (and optionally linear S)
% Slinear = S_l(t, Et);  % commented out as in Python
Squad = S_q(t, Et);

%% Interpolating for smoother plot
% Suppose your existing time-domain vectors are:
% t : original time array (1xN)
% y : original data (1xN)

% Choose how much denser you want the smooth curve:
Nfactor = 10;  % e.g. 10x more points

t=(t - t(find(Squad==max(Squad),1)));

% Create a dense time vector covering the same range:
t_dense = linspace(min(t), max(t), Nfactor * numel(t));

% Interpolate the data using a smooth method ('spline' is usually best):
Squad_dense = interp1(t, Squad, t_dense, 'spline');


%% Plotting
figure('Name','Squadratic');
plot(t * 1e15, Squad, 'LineWidth', 1);
hold on;
plot(t_dense*1e15, Squad_dense, '-', 'LineWidth', 1.5, 'DisplayName', 'Smoothed interpolation');
xlabel('Time delay \tau (fs)', 'FontSize', 12);
ylabel('S_{quadratic} (W/m^2)', 'FontSize', 12);
xlim([-250 250]);
grid on;
% Save figure
print(sprintf('FieldTrace_%db_far.pdf', n_bounces), '-dpdf','-bestfit');

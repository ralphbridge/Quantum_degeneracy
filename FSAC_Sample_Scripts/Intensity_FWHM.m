% Intensity FWHM
%
% Reads in interferometric autocorrelation data and calculates the FWHM
%  of the pulse assuming a Gaussian envelope.
%
% Central wavelength should be set manually below.
%
% Marshall Scott (mscott@thorlabs.com)
% 20170417 - Initial version
%

lam = 0.8;  % [um] Central wavelength
trace_fname = 'waveplate_8b_min_2.csv';
offset_fname = 'waveplate_8b_min_2_background.csv';

% Load interferometric autocorrelation data
M = csvread(trace_fname, 18);
t = M(:,4);  % [s] Time axis data
inf_trace = M(:,5);  % [AU] Interferometric autocorrelation trace

% Load offset data (same settings, laser blocked)
% If no offset data is available, comment this out.
M = csvread(offset_fname, 18);
offset = mean(M(:,5));  % mean offset
inf_trace = inf_trace - offset;  % Subtract offset from trace data

% Shift time axis so that t=0 at signal max
[inf_mx, mx_ind] = max(inf_trace);
t = t - t(mx_ind);

% Plot offset-corrected AC trace.
% format_plot(t*1000, inf_trace, 'Time (ms)', 'Signal (AU)', ...
%            'Interferometric trace');

% Compute the Fourier transform of the autocorrelation trace
dt = t(2) - t(1);  % [s] Compute the time step
f = linspace(-1/2/dt, 1/2/dt, length(t));
inf_fft = ifftshift(fft(fftshift(inf_trace)))*dt;
inf_fft_m = abs(inf_fft);

% Plot FFT of the interferometric AC trace.
format_plot(f, abs(inf_fft), 'Frequency (Hz)', 'Signal (AU)', ...
            'Interferometric trace FFT');

% Find the fundamental frequency
peaks_index = [0; (inf_fft_m(2:end) - inf_fft_m(1:end-1) > 0)] .* ...
              [(inf_fft_m(1:end-1) - inf_fft_m(2:end) > 0); 0];
[~, I] = sort(inf_fft_m .* peaks_index, 'descend');
f_1 = (abs(f(I(2))) + abs(f(I(3)))) / 2;  % [Hz] fundamental frequency
    
% Convert the time axis to optical delay
T = 1 / f_1;  % [s] fringe period
delay = t * f_1 * lam / 0.3;  % [um] optical delay

% Plot offset-corrected AC trace vs delay.
format_plot(delay, inf_trace, 'Delay (fs)', 'Signal (AU)', ...
            'Interferometric trace');

% Filter out the fringes
inf_fft_f = inf_fft .* (abs(f) < f_1/2)';

% Plot filtered FFT of the interferometric AC trace.
%format_plot(f, abs(inf_fft_f), 'Frequency (Hz)', 'Signal (AU)', ...
%            'Filtered interferometric trace FFT');
        
% Inverse FFT to get the intensity autocorrelation with background
int_bg = fftshift(ifft(ifftshift(inf_fft_f)))/dt;
int_bg = abs(int_bg);  % remove phase information

% Remove the background to get the background free intensity trace.
% use average of 10 points from each end.
bg = (sum(int_bg(1:5)) + sum(int_bg(end-4:end)))/10;
int_trace = int_bg - bg;

% Plot offset-corrected AC trace vs delay.
format_plot(delay, int_trace, 'Delay (fs)', 'Signal (AU)', ...
            'Intensity autocorrelation trace');

% Calculate the FWHM of the intensity plot.

tau_ac = fwhm(delay, int_trace);  % [fs] BG free intensity AC FWHM
tau_gauss = tau_ac / sqrt(2);  % Pulse width assuming Gaussian shape
tau_sec = tau_ac / 1.543;  % Pulse width assuming sec^2 shape
disp(' ')
disp(['bg-free Intensity AC FWHM: ', num2str(tau_ac), ' fs'])
disp(['Gaussian pulse FWHM: ', num2str(tau_gauss), ' fs'])
disp(['Sech^2 pulse FWHM: ', num2str(tau_sec), ' fs'])
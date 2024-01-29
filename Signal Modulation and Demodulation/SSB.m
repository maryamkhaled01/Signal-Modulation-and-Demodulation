%% 1. Original Message

% Read the Sound File
[voice,fs]=audioread('eric.wav'); % fs = 48kHz
sec = 8;
voice = voice(1:sec*fs,1); 
time=linspace(0,sec,sec*fs);

% Plot Original Message in Time Domain
figure('Name','Original Signal','NumberTitle','on');
subplot(2,1,1);
plot(time,voice)
title('Message in Time Domain');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Compute the spectrum
voice_spectrum=fftshift(fft(voice));
freq=linspace(-fs/2,fs/2,length(voice_spectrum));

% Plot Original Message Spectrum
subplot(2,1,2);
plot(freq,abs(voice_spectrum))
title('Message in Frequency Domain')
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Playing Sound
sound (voice,fs);
pause(sec);
%% 2 & 3. Ideal Filter

cutoff_frequency = 4000;
filter_order = 50;

% Signal Filtering
fltr=designfilt('lowpassfir','FilterOrder',filter_order,'CutoffFrequency',cutoff_frequency, 'SampleRate',fs);
filtered_voice=filter(fltr,voice);

% Plot Filtered Message in Time Domain
figure('Name','Filtered Message','NumberTitle','on');
subplot(2,1,1);
plot(time,filtered_voice)
title('Filtered Message in Time Domain');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot Filtered Message Spectrum
subplot(2,1,2);
plot(linspace(-fs/2,fs/2,length(filtered_voice)),abs(fftshift(fft(filtered_voice))));
title('Filtered Message in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Playing Sound
sound (filtered_voice,fs);
pause(sec);
%% 4.DSB-SC

% Upsample
fc = 100000; % Carier Frequency
fs_new = 5*fc ; % Sampling Frequency (new fs =500k)
resampled_voice=resample(filtered_voice,fs_new,fs); % Resmaple by 500k/48k
time=linspace(0,sec,sec*fs_new);

% Message Modulation
Ac = 1; 
carrier = Ac .* cos(2*pi*fc*time'); % Carrier Signal
suprsd_carrier = carrier.*resampled_voice; % DSB-SC

% Plot DSB-SC in Time Domain
figure('Name','DSB-SC','NumberTitle','on');
subplot(2,1,1);
plot(time,suprsd_carrier);
title('DSB-SC in Time Domain');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot DSB-SC Spectrum
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(suprsd_carrier)),abs(fftshift(fft(suprsd_carrier))))
title('DSB-SC in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Playing Sound
% sound(suprsd_carrier,fs_new)
% pause(sec);
%% 5. SSB-SC

% Design an ideal bandpass filter
f_passband = [fc - 4000 , fc]; % Passband frequencies
filter_order = 10000; % Filter order
filter_coeffs_mod = fir1(filter_order, f_passband / (fs_new/2), 'bandpass', hamming(filter_order + 1));

% Apply the filter to the DSB-SC signal
LSB = filter(filter_coeffs_mod, 1, suprsd_carrier);

% Plot LSB in Time Domain
figure('Name','LSB','NumberTitle','on');
subplot(2,1,1);
plot(time',LSB)
title('LSB in Time Domain');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot LSB Spectrum
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(LSB)),abs(fftshift(fft(LSB))))
title('LSB in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
%% 6. Coherent Detection

% Signal Demodulation
t = linspace(0, length(LSB)/fs_new, length(LSB));
demodulated_signal = LSB .* carrier;

% Plot Demoodulated LSB (before LPF) in Time Domain
figure('Name','Demodulated LSB (Before LPF)','NumberTitle','on');
subplot(2,1,1);
plot(t,demodulated_signal) 
title('Demodulated LSB in Time Domain (Before LPF)');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot Demoodulated LSB (before LPF) in Frequency Domain
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(demodulated_signal)),abs(fftshift(fft(demodulated_signal))))
title('Demodulated LSB in Frequency Domain (Before LPF)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Design an ideal LPF
f_cutoff = fs; % Cutoff frequency in Hz
filter_order = 1000; % Filter order 
filter_coeffs = fir1(filter_order, f_cutoff / (fs_new/2), 'low');

% Apply LPF on the demodulated signal to remove interferences
filtered_demodulated_signal = filter(filter_coeffs, 1, demodulated_signal);

filtered_t = linspace(0, length(filtered_demodulated_signal)/fs_new, length(filtered_demodulated_signal));
freq_axis_filtered = linspace(-fs_new/2, fs_new/2, length(filtered_demodulated_signal));

% Plot Demoodulated LSB (After LPF) in Time Domain
figure('Name','Demodulated LSB (After LPF)','NumberTitle','on');
subplot(2,1,1);
plot(filtered_t, filtered_demodulated_signal);
title('Demodulated LSB in Frequency Domain (After LPF)');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot Demoodulated LSB (After LPF) in Frequency Domain
subplot(2,1,2);
plot(freq_axis_filtered, abs(fftshift(fft(filtered_demodulated_signal))));
title('Demodulated LSB in Time Domain (After LPF)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Playing Sound
r_demodulated_signal=resample(filtered_demodulated_signal,fs,fs_new);
sound(r_demodulated_signal,fs)
pause(sec);
%% 7. Butterworth Filter

% Buttterworth Filter to get LSB
[b,a]= butter(4,[((fc-4000)/(fs_new/2)) (fc/(fs_new/2))],'bandpass');
butter_LSB = filter(b, a, suprsd_carrier);

% Plot Butter LSB in Time Domain
figure('Name',' LSB (Butterworth)','NumberTitle','on');
subplot(2,1,1);
plot(time,butter_LSB)
title(' LSB (Butterworth) in Time Domain');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot Butter LSB in Frequency Domain
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(butter_LSB)),abs(fftshift(fft(butter_LSB))))
title(' LSB (Butterworth) in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Coherent Detection 
butter_t = linspace(0, length(butter_LSB)/fs_new, length(butter_LSB));
butter_demodulated_signal = butter_LSB .* carrier;

% Plot Demoodulated butter_LSB (before LPF) in Time Domain
figure('Name','Demodulated LSB (Butterworth) (Before LPF)','NumberTitle','on');
subplot(2,1,1);
plot(butter_t,butter_demodulated_signal) 
title('Demodulated LSB (Butterworth) in Time Domain (Before LPF)');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot Demoodulated butter_LSB (before LPF) in Frequency Domain
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(butter_demodulated_signal)),abs(fftshift(fft(butter_demodulated_signal))))
title('Demodulated LSB (Butterworth) in Frequency Domain (Before LPF)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Apply LPF (in no.6) on the demodulated signal to remove interferences
filtered_butter_demodulated_signal = filter(filter_coeffs, 1, butter_demodulated_signal);

filtered_t = linspace(0, length(filtered_butter_demodulated_signal)/fs_new, length(filtered_butter_demodulated_signal));
freq_axis_filtered = linspace(-fs_new/2, fs_new/2, length(filtered_butter_demodulated_signal));

% Plot Demoodulated LSB (After LPF) in Time Domain
figure('Name','Demodulated LSB (Butterworth) (After LPF)','NumberTitle','on');
subplot(2,1,1);
plot(filtered_t, filtered_butter_demodulated_signal);
title('Demodulated LSB (Butterworth) in Frequency Domain (After LPF)');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot Demoodulated LSB (After LPF) in Frequency Domain
subplot(2,1,2);
plot(freq_axis_filtered, abs(fftshift(fft(filtered_butter_demodulated_signal))));
title('Demodulated LSB (Butterworth) in Time Domain (After LPF)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Playing Sound
r_butter_demodulated_signal=resample(filtered_butter_demodulated_signal,fs,fs_new);
sound(r_butter_demodulated_signal,fs)
pause(sec);
%% 8. Ideal Filter Case with Noise

% SNR values to simulate
snr_values = [0, 10, 30];

% Plot demodulated waveform and spectrum for each SNR
for snr_value = snr_values
    % Add noise to the recieved LSB
    noisy_LSB = awgn(LSB, snr_value);
    
    % Time vector of noisy signal
    t_noise = linspace(0, length(noisy_LSB)/fs_new, length(noisy_LSB)); 
    local_carrier = Ac .* cos(2*pi*fc*t_noise'); % Carrier Signal
    
    % Coherent detection
    noisy_demodulated_signal = noisy_LSB .* local_carrier;
    
    % Apply LPF (in no.6) on the demodulated signal to remove interferences
    filtered_noisy_demodulated_signal = filter(filter_coeffs, 1, noisy_demodulated_signal);

    % Plot received demodulated noisy signal in Time Domain
    figure('Name', ['Demodulated Noisy Signal (SNR = ' num2str(snr_value) ' dB)'], 'NumberTitle', 'on');
    subplot(2, 1, 1);
    plot(t_noise, filtered_noisy_demodulated_signal);
    title(['Demodulated Noisy Signal in Time Domain (SNR = ' num2str(snr_value) ' dB)']);
    xlabel('Time (sec)');
    ylabel('Magnitude');
    grid on;

    % Plot received demodulated noisy signal in Frequency Domain
    subplot(2, 1, 2);
    freq_axis = linspace(-fs_new/2, fs_new/2, length(filtered_noisy_demodulated_signal));
    plot(freq_axis, abs(fftshift(fft(filtered_noisy_demodulated_signal))));
    title(['Demodulated Noisy Signal in Frequency Domain (SNR = ' num2str(snr_value) ' dB)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on; 
    
    % Play back the demodulated signal
    r_noisy_demodulated_signal=resample(filtered_noisy_demodulated_signal,fs,fs_new);
    sound(r_noisy_demodulated_signal, fs);
    pause(sec);
    
end
%% 9. SSB-TC

% DC_bias
dc_bias = 2 .* max(resampled_voice);
trnsm_carrier=(dc_bias + resampled_voice) .* carrier; % DSB-TC

% Apply the bandpass filter (in no. 5) to the DSB-TC signal
LSB_tc = filter(filter_coeffs_mod, 1, trnsm_carrier); % SSB-TC

% Plot LSB in Time Domain
figure('Name','LSB-TC','NumberTitle','on');
subplot(2,1,1);
plot(time',LSB_tc)
title('LSB-TC in Time Domain');
xlabel('Time (sec)');
ylabel('Magnitude');
grid on;

% Plot LSB Spectrum
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(LSB_tc)),abs(fftshift(fft(LSB_tc))))
title('LSB-TC in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Envelope Detector
envelope = abs(hilbert(LSB_tc));

% Plot the Demodulated Signal in Time Domain
figure('Name', 'Envelope (SSB-TC)', 'NumberTitle', 'on');
subplot(2,1,1);
plot(time, envelope);
title('Envelope (SSB-TC) in Time Domain');

% Plot the Demodulated Signal in Frequency Domain
subplot(2,1,2);
plot(linspace(-fs_new/2,fs_new/2,length(envelope)),abs(fftshift(fft(envelope))))
title('Envelop (SSB-TC) in Frequency Domain');

% Playing Sound
envelope = resample(envelope,fs,fs_new);
sound(envelope,fs)
pause(sec);
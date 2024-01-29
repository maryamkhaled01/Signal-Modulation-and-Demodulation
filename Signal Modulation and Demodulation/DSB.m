
%% PART 1
clear all;
%% 1

filename = 'eric.wav';
[audio, fs] = audioread(filename);
sound(audio,fs);

len = length(audio);
audio_freq = fftshift(fft(audio));
f_axis = fs/2*linspace(-1,1,len);

% Plot the spectrum
figure;
plot(f_axis, abs(audio_freq)/len);
xlabel('f(Hz)');
ylabel('Magnitude');
title('unfiltered signal in frequency domain');


%% 2 and 3

% Filter as 4kHz
BW = 4000;
audio_freq(f_axis >= BW | f_axis <= -BW) = 0;
filtered_signal = ifft(ifftshift(audio_freq));
len = length(filtered_signal);
audio_freq = fftshift(fft(filtered_signal));
f_axis = fs/2*linspace(-1,1,len);

% Plot the filtered spectrum
figure;
plot(f_axis, abs(audio_freq)/len);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Filtered Spectrum');

t1 = linspace(0,len/fs,len); 
t1=t1';
% Plot the filtered waveform
figure;
plot(t1, filtered_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Signal in the Time Domain');


%% 4 sound
sound(abs(filtered_signal),fs);

%% 5

fc = 100000;
m = 0.5;
Am=max(filtered_signal);
Ac = Am/m;

filtered_signal  = resample(filtered_signal,5*fc,fs);
fs=5*fc;
t1=linspace(0,length(filtered_signal)/fs,length(filtered_signal));
t1=t1';
carrier = Ac.*cos(2*pi*fc*t1);

% DSB-TC
dsbtc = carrier.*(1+m*filtered_signal/Am);

% DSB-SC
dsbsc = carrier .* filtered_signal;

spectrum_dsbtc = fftshift(fft(dsbtc));
spectrum_dsbsc = fftshift(fft(dsbsc));

len = length(dsbtc);
f_axis=fs/2*linspace(-1,1,len);

% Plot the spectrum of the modulated signals
figure;
plot(f_axis, abs(spectrum_dsbtc)/len);
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');
title('DSB-TC Modulated Signal Spectrum');

len = length(dsbsc);
f_axis=fs/2*linspace(-1,1,len);

figure;
plot(f_axis, abs(spectrum_dsbsc)/len);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-SC Modulated Signal Spectrum');

%% 6 Envelope

DSBSC_envelope = abs(hilbert(dsbsc));
DSBTC_envelope = abs(hilbert(dsbtc));
%% Plot
figure; 
plot(t1,dsbsc);
hold on;
plot(t1,-DSBSC_envelope,'k-',t1, DSBSC_envelope,'-k','Linewidth',1.5); % phase reversal occurs
hold off;
title('DSBSC (in blue) in time domain with envelope detector (in black)');
ylim([-2 2]);
xlim([2 2.5]);

figure; 
plot(t1,dsbtc);
hold on;
plot(t1,-DSBTC_envelope,'k',t1,DSBTC_envelope,'-k','Linewidth',1.5); % no phase reversal occurs
title('DSBTC (in blue) in time domain with envelope detector (in black)');
hold off;
ylim([-5 5]);
xlim([3 3.5]);




%% 7 Resample and sound for DSB-SC

DSBSC_envelope = resample(abs(DSBSC_envelope), fs/5, fs);
sound(abs(DSBSC_envelope), fs/5);
% Observation: not detected very well, distorted sound
% Envelope detector can only be used for DSB-TC


%% Resample and Sound for DSB-TC

DSBTC_envelope = resample(abs(DSBTC_envelope), fs/5, fs);
sound(abs(DSBTC_envelope), fs/5);
% Observation: more accurately detected, less distortion.



%% 8
SNR = [0,10,30];
for i=1:length(SNR)
    
    % generate signal+noise
    dsbsc_noise = awgn(dsbsc, SNR(i));
    
    % demodulate using coherent detector
    demodulatedsignal = dsbsc_noise.*cos(2*pi*fc*t1);
    
    clear dsbsc_noise;
    
    % fourier transform
    demodulatedsignal_in_FD = fftshift(fft(demodulatedsignal));
    
    % LPF at Fm
    demodulatedsignal_in_FD(f_axis >= BW | f_axis <= -BW) = 0;
    
    % inverse fourier transform to get demodulated signal in time domain
    demodulatedsignal = ifft(ifftshift(demodulatedsignal_in_FD));

    % plot demodulated signal in time domain
    figure; plot(t1, demodulatedsignal); title([num2str(SNR(i)),' SNR demodulated signal in time domain']);

    % fourier transform
    len = length(demodulatedsignal);
    spectrum_demodulated = fftshift(fft(demodulatedsignal));
    f_axis = fs/2*linspace(-1,1,len);
    
    % plot demodulated signal in frequency domain
    figure; plot(f_axis, abs(spectrum_demodulated) / len); title([num2str(SNR(i)),' SNR demodulated signal in frequency domain']);

    % resample to sound the demodulated signal
    demodulatedsignal = resample(abs(demodulatedsignal), fs/5, fs);
    sound(abs(demodulatedsignal), fs/5);
    pause(10);
end
clear demodulatedsignal;
clear demodulatedsignal_in_FD;
%% 9 with frequency error

fc = 100100;
for i=1:length(SNR)
    
    % generate signal+noise
    dsbsc_noise = awgn(dsbsc, SNR(i));
    
    % demodulate using coherent detector
    demodulatedsignal = dsbsc_noise.*cos(2*pi*fc*t1);
    
    clear dsbsc_noise;
    
    % fourier transform
    demodulatedsignal_in_FD = fftshift(fft(demodulatedsignal));
    
    % LPF at Fm
    demodulatedsignal_in_FD(f_axis >= BW | f_axis <= -BW) = 0;
    
    % inverse fourier transform to get demodulated signal in time domain
    demodulatedsignal = ifft(ifftshift(demodulatedsignal_in_FD));

    % plot demodulated signal in time domain 
    figure; plot(t1, demodulatedsignal); title([num2str(SNR(i)),' SNR demodulated signal with frequency error in time domain']);

    % fourier transform
    len = length(demodulatedsignal);
    spectrum_demodulated = fftshift(fft(demodulatedsignal));
    f_axis = fs/2*linspace(-1,1,len);
    
    % plot demodulated signal in frequency domain
    figure; plot(f_axis, abs(spectrum_demodulated) / len); title([num2str(SNR(i)),' SNR demodulated signal with frequency error in frequency domain']);
    % resample to sound the demodulated signal
    demodulatedsignal = resample(demodulatedsignal, fs/5,fs);
    sound(abs(demodulatedsignal), fs/5);
    
    pause(10);
    
end

clear demodulatedsignal;
clear demodulatedsignal_in_FD;

%% 10 with phase error 
fc = 100000;
for i=1:length(SNR)
    
    % generate signal+noise
    dsbsc_noise = awgn(dsbsc, SNR(i));
    
    % demodulate using coherent detector
    demodulatedsignal = dsbsc_noise.*cos(2*pi*fc*t1 + pi/9);
    
    clear dsbsc_noise;
    
    % fourier transform
    demodulatedsignal_in_FD = fftshift(fft(demodulatedsignal));
    
    % LPF at Fm
    demodulatedsignal_in_FD(f_axis >= BW | f_axis <= -BW) = 0;
    
    % inverse fourier transform to get demodulated signal in time domain
    demodulatedsignal = ifft(ifftshift(demodulatedsignal_in_FD));

    % plot demodulated signal in time domain 
    figure; plot(t1, demodulatedsignal); title([num2str(SNR(i)),' SNR demodulated signal with phase error in time domain']);

    % fourier transform
    len = length(demodulatedsignal);
    spectrum_demodulated = fftshift(fft(demodulatedsignal));
    f_axis = fs/2*linspace(-1,1,len);
    
    % plot demodulated signal in frequency domain
    figure; plot(f_axis, abs(spectrum_demodulated) / len); title([num2str(SNR(i)),' SNR demodulated signal with phase error in frequency domain']);

    % resample to sound the demodulated signal
    demodulatedsignal = resample(demodulatedsignal, fs/5,fs);
    sound(abs(demodulatedsignal), fs/5);
    pause(10);
    
end

clear demodulatedsignal;
clear demodulatedsignal_in_FD;





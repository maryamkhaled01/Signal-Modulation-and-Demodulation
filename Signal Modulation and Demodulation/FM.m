clear;
clear all;
close all;

%Read audio file
[S, Fs] = audioread('eric.wav');
%sound(abs(S),Fs);

%Find the spectrum
%Fourier transform
L = length(S);
F = fftshift(fft(S));
f = Fs/2*linspace(-1,1,L);

%Plotting the spectrum
figure; plot(f,abs(F)/L); title('The original signal before filteration'); %abs(F)/L -> spectrum magnitude

%Filtering and plot
cutoff_frequency = 4000;
F(f>=cutoff_frequency|f<=-cutoff_frequency) = 0;
figure; plot(f,abs(F)/L); xlim([-5000,5000]); title('Spectrum of filtered signal');

%In time domain
X = ifft(ifftshift(F));
%sound(abs(X),Fs);

%calculate time
tStart = 0;
tEnd = tStart + length(X) / Fs;
t=linspace(tStart,tEnd,length(X));
t=t';
figure; plot(t,X); title('Time of filtered signal');
%sound(abs(X),Fs);

%constants
fc=100000;
Ac=1;
omega_c=2*pi*fc;
%Get maximum deviation of the integrated(X) signal from its mean.
max_dev=max(abs(cumsum(X)));
%normalization of the deviation
norm=2*pi*max_dev./Fs;
k_fm=0.2/norm;

%resampling
X=resample(X,5*fc,Fs);
Fs=5*fc;

tStart = 0;
tEnd = tStart + length(X) / Fs;
t=linspace(tStart,tEnd,length(X));
t=t';

%FM modulation
X=Ac*cos(omega_c*t + 2*pi*k_fm*cumsum(X)./Fs);

%Fourier transform
L = length(X);
F = fftshift(fft(X));
f = Fs/2*linspace(-1,1,L);

%Plotting the spectrum
figure; plot(f,abs(F)/L); title('FM modulation spectrum');

%descriminator
dy=diff(X); %calculating difference between X samples
dy=[0;dy]; %Making sure to let the length be same as the original X

% envelope detector
envelopeFM = abs(hilbert(dy)) -  mean(abs(hilbert(dy)));

%plotting for FM demod. signal in time
figure; plot(t,envelopeFM); title('FM demodulated signal in Time');
ylim([-0.00015 0.00015]);

% resample to hear the signal
envelopeFM = resample(envelopeFM, Fs/5, Fs);
%sound(100.*abs(envelopeFM), Fs/5);
sound(1000.*abs(envelopeFM), Fs/5);





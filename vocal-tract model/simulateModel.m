% Using simulateBMP in 2D models of speech synthesis
% First section runs a simulation to obtain the impulse response
% Second section convolves a glottal pulse train (Rosenberg type) with the
% impulse response to obtain the synthesized vowel
% ------------------- Jorge Petrosino (june, 2022)
clear; clc

imageFileName='TubeModel_o.bmp';

scale=1e-2;
duration=20e-3;
simSource.type = 'impulse';
simSource.amplitude = 100;      
simSource.fCut=15000;
simMedium.c0 = 344; %Default value
simMedium.density = 1.2; %Default value
simMedium.CFL=0.5;
%simMedium.WallMask=true;
recordVideo=false;

[sensorData, t, dt, pVariable, lx, ly] = ...
       simulateBMP(imageFileName, scale, duration, ...
       simMedium, simSource,recordVideo);
%%
% Logical flags
writeAudioFile=true;
flagGlottal=true; %true uses Rosenberg pulses, false uses pulse train

% Audio parameters
fskWave=round(1/dt); % simulation sample rate
fs48k=48000; % audio sample rate
dt48k=1/fs48k; 
t48k=0:dt48k:t(end)+dt-dt48k; % time vector of audio
f48k=(0:length(t48k)-1)/length(t48k)*fs48k; % freq vector of audio

% This windowing is tu fade out the impulse response
% The surfaces of the simultacion are ideal relflection surfaces. This
% means the impulse response duration is longer than expected.
h=sensorData.*exp(-t/5e-3); %windowing the response (exponential fade out)


% resample impulse response (from fskWave to fs48k)
h48k=interp1(t,h,t48k);

% Fourier Transform of impulse response (Transference function H48k)
H48k=fft(h48k)*2/length(h48k);

% Generate glottal pulse train
durationGlottal=0.5;
tglottal=0:dt48k:durationGlottal-dt48k;
fglottal=(0:length(tglottal)-1)/length(tglottal)*fs48k;
f0=100;
if flagGlottal
    q=4; 
    % Uses Rosenberg pulses
    glottalTrain=glottalSource(f0,durationGlottal,fs48k,q);
else
    nt=round(tglottal/dt48k);
    T=1/f0; nT=round(T/dt48k);
    % Uses an ideal pulse train
    nsamplesDutyCycle=2;
    glottalTrain=double(mod(nt,nT)<=(nsamplesDutyCycle-1));
    %Using nsamplesDutyCucle=2 gives a freq response droping down at
    %fNyquist=fs48k/2
end

% Convolving glottal pulses whith impulse response
audio=conv(h48k,glottalTrain);
maxAudio=max(audio);
audio=audio/maxAudio; % Normalize audio levels
tConv=0:dt48k:(length(audio)-1)*dt48k;
Audio=fft(audio)*2/length(audio);
fConv=(0:length(audio)-1)/length(audio)*fs48k;
% Fourier transform of glottalTrain
GlottalTrain=fft(glottalTrain)*2/length(glottalTrain);

% Graphs -----------------------------------------------------
figure(1)

subplot(3,2,1)
maxh48k=max(h48k);
plot(t48k*1e3,h48k)
axis([0 inf -1.1*maxh48k 1.1*maxh48k])
grid on; grid minor
title('impulse response (h)')
xlabel('t [ms]')

subplot(3,2,2)
maxH=max(abs(H48k));
semilogx(f48k,20*log10(abs(H48k/maxH)))
%axis([100 fs48k/2 maxH-30 maxH+10])
axis([100 fs48k/2 -30 10])
grid on; grid minor
xticks([200 1000 10000])
xlabel('f [Hz]')
ylabel('dB [rel]')

subplot(3,2,3)
plot(tglottal*1e3,glottalTrain)
axis([0 50 -1.1 1.1])
grid on; grid minor
title('glottal pulse (g)')
xlabel('t [ms]')

subplot(3,2,4)
maxG=max(abs(GlottalTrain));
semilogx(fglottal,20*log10(abs(GlottalTrain/maxG)+1e-8))
axis([100 fs48k/2 -30 10])
grid on; grid minor
xticks([200 1000 10000])
xlabel('f [Hz]')
ylabel('dB [rel]')

maxaudio=max(abs(audio));
audio=audio/maxaudio*0.8; % Normalize the audio
subplot(3,2,5)
plot(tConv*1e3,audio)
axis([0 50 -1.1*maxaudio 1.1*maxaudio])
grid on; grid minor
title('conv(g,h)')
xlabel('t [ms]')

maxA=max(abs(Audio));
subplot(3,2,6)
semilogx(fConv,20*log10(abs(Audio/maxA)+1e-6))
axis([100 fs48k/2 -30 10])
grid on
xticks([200 1000 10000])
xlabel('f [Hz]')
ylabel('dB [rel]')

% Play audio samples
% 1: glottal pulses (for reference), 2: sound result of processing
maxGlotis=max(glottalTrain);
silence500ms=zeros(1,24000);
sound([glottalTrain/maxGlotis, silence500ms, audio]*0.55,fs48k)

% Saving audio to wav file
if writeAudioFile
    audiowrite([imageFileName(1:end-4) '.wav'],audio,fs48k);
end


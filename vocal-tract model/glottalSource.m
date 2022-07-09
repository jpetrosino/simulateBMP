function [glottalPulseTrain]=glottalSource(f0,duration,fs,q)
% function [glottalPulseTrain]=glottalSource(f0,duration,fs,q);
%                            coded by Jorge Petrosino
% Generate audio signal with a train of glottal pulses, after
% Rosenberg, A. E. (1971). 
% Effect of glottal pulse shape on the quality of natural vowels. 
% The Journal of the Acoustical Society of America, 49(2B), 583-590.
%
% f0 is the fundamental frequency of the pulse train
% duration desired of audio in seconds
% fs is the sampling frequency of the audio
% q is a parameter related to the sharpness of each pulse. Greater values
% of q means more sharp waveform. Recommended value 8 (suggested 1<=q<=32)
dt=1/fs;
t=0:dt:duration-dt;
%fglotis=(0:length(t)-1)/length(t)*fs48k;
nt=round(t/dt);
T=1/f0; nT=round(T/dt);
T1=0.25*T/q; nT1=round(T1/dt);
T2=0.1*T/q; nT2=round(T2/dt);

% Glottal pulse 
glottalPulseTrain= ...
   and(0<=mod(nt,nT),mod(nt,nT)<=mod(nT1,nT)).* ...
     0.5.*(1-cos(2*pi*t/(2*T1))) ...
   +and(mod(nT1,nT)<mod(nt,nT),mod(nt,nT)<=mod(nT1+nT2,nT)).* ...
     abs(cos( 2*pi*(t-T1)/(4*T2) ));

end
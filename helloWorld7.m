% Ejemplo de versión mínima (utilizando valores por defecto)
%% 1-Setup parameters
imageFileName='HelloWorld7.bmp';
scale=1e-2; % 1 cm
duration = 10e-3; % 10 ms
simSource.type='nCycles';
simSource.n=0.5;
simSource.f0=8000;
simMedium.c0=344; % 344 m/s
simMedium.alpha=0.8;
recordVideo=false;

%% 2-Call simulateImage256
[sensor_data, t, dt, pVariable, lx, ly] = ...
          simulateBMP(imageFileName, scale, duration, ...
          simMedium, simSource,recordVideo);

%% 3-Plot the results
figure(1)
tVariable=(0:length(pVariable)-1)*dt;
plot(t,sensor_data)
title('Sensor'); xlabel('t [s]'); ylabel('p [Pa]'); grid on

% Ejemplo de versión mínima (utilizando valores por defecto)
%% 1-Setup parameters
imageFileName='HelloWorld2.bmp';
scale=1e-2; % 1 cm
duration = 10e-3; % 10 ms
%simSource.fs=48000; dt=1/simSource.fs; t=(0:dt:duration-dt);
simSource.fs=48000; dt=1/simSource.fs; t=(0:dt:duration-dt);
simSource.type='audio';
simSource.p=10*sin(2*pi*1000*t);
simMedium.c0=344; % 344 m/s
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

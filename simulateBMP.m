function [sensorData, t, dt, pVariable, lx, ly] = ...
          simulateBMP(imageFileName, scale, duration, ...
          simMedium, simSource,recordVideo)

%%     simulateBMP  v3.0 2020 (upgrade from simulateImage256)
%
% [sensorData, t, dt, pVariable, lx, ly] = ...
%          simulateBMP(imageFileName, scale, duration, ...
%          simMedium, simSource,recordVideo);
%
% This function runs time-domain simulations of acoustic wave propagation
% based on the information contained in a bitmap image file.                                     
%
% Designed to work with the k-Wave toolbox, available for free at: 
% http://www.k-wave.org/
%
%--------------------------------------------------------------------------
% Input arguments
%--------------------------------------------------------------------------
% - Name of the 24-bits or 256-colour BMP image file to be used.
%   Recomended size 2^n x 2^m (fast computation), for example 256x512
%   imageFilaName = 'name.bmp' 
%
% - Desired side measurement of the minimum square of the grid, in metres.
%   scale = 1e-2; %This scale ensures fNyquist = 17 kHz
%
% - Desired duration of the simulation, in seconds.
%   duration = 2e-3;
%
% . Medium properties struct
%   -  Sound speed value, in metres per second.
%      simMedium.c0 = 344; %Default value
%   -  Air density, in kg/m3.
%      simMedium.density = 1.2; %Default value
%   -  Sabine Absorption coefficient.
%      simMedium.alpha = 0.5; %Defauult value (1 means total absorption)
%
% - Source struct. Explained in the following section.
%
% - Option to record the simulation to a video file, true or false.
%   recordVideo = true;
%
%
%--------------------------------------------------------------------------
% Source properties
%--------------------------------------------------------------------------
% The source argument is a structure arrays that admit 5 different formats. 
% These can be set by assigning certain values to the source struct.
% 
% 1) Impulse:
%    simSource.type = 'impulse';
%    simSource.amplitude = 4;        (desired amplitude in pascals)
%    simSource.mode = 'dirichlet';   (mandatory pressure on source point)
%
% 2) Sine wave (n cycles): 
%    simSource.type = 'nCycles';     
%    simSource.amplitude = 4;        (desired amplitude in pascals)
%    simSource.n=2;                  (desired number of cicles)
%    simSource.f0 = 1000; [Hz]       (desired frequency)
%    simSource.mode = 'additive'; (pressure aded at source point)
%
% 3) White noise:
%     simSource.type = 'whiteNoise'; 
%     simSource.amplitude = 2; 
%     simSource.mode = 'additive';
%     simSource.duration = 2e-3; [s] (desired noise duration)
%
% 4) Free form (any mathematical expression accepted by MATLAB):
%     simSource.type = '4*sin(2*pi*1000*t+pi/4)'; (example)
%     simSource.mode = 'additive':
%
% 5) Free variable (e.g. audio = any one-dimension MATLAB variable)
%     simSource.type = 'audio'
%     simSource.mode='dirichlet';
%     simSource.fs = 44100; % sample rate of audio variable
%     simSource.p=audioSamples % must be the MATLAB variable
%     sensorData is returned using simSource.fs when using 'audio'
%
%     Sample frequency inide the simulation is determined by fs=1/dt
%     (scale=1e-2 m, determines fs= 114667 Hz, and fNyquist=17 kHz)
%
%--------------------------------------------------------------------------
% Creation of sources, sensors and refletive surfaces
%--------------------------------------------------------------------------
% RED pixels are treated as SOURCES.
% GREEN pixels are treated as SENSORS.
% BLACK pixels are treated as perfectly reflective SURFACES. 
% BROWN pixels are partially reflectant surfaces.
%
% Colour values are as follows (standard palette of Paint):
% red256 = 79; red = [237    28    36];
% green256 = 113; green =[34   177    76];
% black256 = 0; black = [0 0 0];
% brown256 = 1; brown = [136     0    21];
%
% These are the standard colour codes used by MS Paint, the most commonly
% available image editing software.
%
% The rest of the colours are ignored by the simulation, and can be used to 
% add notes and other useful information to the simulation's animation,
% without it affecting the results.
%
%--------------------------------------------------------------------------
% Reference paper
%--------------------------------------------------------------------------
% "MATLAB-based simulation software as teaching aid for physical acoustics"
% Jorge Petrosino, Lucas Landini, Georgina Lizaso, Ian Kuri, Ianina Canalis
% Universidad Nacional de Lanus, Argentina
% 23rd International Congress on Acoustics, 2019.
% http://pub.dega-akustik.de/ICA2019/data/articles/001342.pdf
% https://youtu.be/aKQ2G7RMh64
% 
% Sample simulations and complementary functions available at:
% https://github.com/GLizaso/Teaching_aid_for_physical_acoustics


%% Image loading and format checking

image = imread(imageFileName);
format = imfinfo(imageFileName);

if not(or(format.BitDepth == 8, format.BitDepth == 24)) 
    disp('Unsupported image file. 24-bits or 256-colour BMP files must be used');
end 

if not(strcmp(format.Format,'bmp')) 
    disp('Unsupported image file. 24-bits or 256-colour BMP files must be used');
end 

%% Grid settings

%CFL = 0.3; Se define más abajo
dx = scale; dy = scale;
if format.BitDepth == 8
    Nx = length(image(:,1));  Ny = length(image(1,:)); 
end
if format.BitDepth == 24
    Nx = length(image(:,1,1));  Ny = length(image(1,:,1)); 
end
lx = Nx*dx; ly = Ny*dy;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% Medium properties
%   -  Sound speed value, in metres per second.
if not(isfield(simMedium,'c0')); simMedium.c0=344; end %INUTIL POR AHORA YA QUE EL STRUCT DEBE EXISTIR
if not(isfield(simMedium,'density')); simMedium.density=1.2; end % Deafult value
if not(isfield(simMedium,'alpha')); simMedium.alpha=0.5; end
if not(isfield(simMedium,'speedRatio')); simMedium.speedRatio=1.1; end
if not(isfield(simMedium,'sign')); simMedium.sign=1; end
if not(isfield(simMedium,'CFL')); simMedium.CFL=0.3; end
if not(isfield(simMedium,'showAbsorptionMask')); simMedium.showAbsorptionMask=true; end
if not(isfield(simMedium,'showWallMask')); simMedium.showWallMask=true; end
if not(isfield(simMedium,'showSourceMask')); simMedium.showSourceMask=true; end
if not(isfield(simMedium,'showSensorMask')); simMedium.showSensorMask=true; end
if not(isfield(simMedium,'showLegendMask')); simMedium.showLegendMask=true; end

airSpeed = simMedium.c0; %Default value
%   -  Air density, in kg/m3.
airDensity=simMedium.density; %Default value
%   -  Sabine Absorption coefficient.
alphaSabine=simMedium.alpha; %Defauult value (1 means total absorption)

%%%%%%%%% REVISAR
%medium.density = 1.24;    % Air density is the default. Can be modified.
%medium.sound_speed = c0;  % c0 is a user-defined input argument

%%%%%REVISAR
%airDensity=1.24; airSpeed=c0; 
% Para luego agregar wallDensity; wallSoundSpeed;

%% Absorptive material

%%%%%%REVISAR brown = 1; % Brown pixels must have this value to work as sensors

brown256=1; brown=[136 0 21];

if format.BitDepth==8
    absorptiveMask = (image == brown256);
else
    absorptiveMask = (squeeze(image(:,:,1)) == brown(1)&...
                    squeeze(image(:,:,2)) == brown(2)&...
                    squeeze(image(:,:,3)) == brown(3));
end
if sum(sum(absorptiveMask)==0); AbsorptionExist=false; else; AbsorptionExist=true; end
% Next mask is used to show absorptive material in the simulation
pointMask=zeros(Nx,Ny); 
pointMask(1:4:end,1:4:end)=1; pointMask(3:4:end,3:4:end)=1;

%%%% REVISAR . QUIZÁS DEBA CAMBIAR wallSpeed según alfa
wallSpeed=simMedium.speedRatio*airSpeed;
% Positive reflection coefficient => simMedium.sign=1; 
% Negative reflection coefficient => simMedium.sign=-1; 
reflCoef=simMedium.sign*sqrt(1-alphaSabine);
wallDensity=airSpeed/wallSpeed*airDensity*(1+reflCoef)/(1- reflCoef);
% Negative reflection coefficient
%wallDensity=airSpeed/wallSpeed*airDensity*(1-reflCoef)/(1+ reflCoef);

medium.density=absorptiveMask*wallDensity + (~absorptiveMask)*airDensity;
medium.sound_speed=absorptiveMask*wallSpeed + (~absorptiveMask)*airSpeed;
z0Air=airDensity*airSpeed;
z0Wall=wallDensity*wallSpeed;

%{
airAlphaCoeff=0.00;
wallAlphaCoeff=0; %0.00005;
medium.alpha_coeff = absorptiveMask*wallAlphaCoeff + (~absorptiveMask)*airAlphaCoeff;
medium.alpha_power = 1.5;
%}

%% Time array creation

kgrid.t_array = makeTime(kgrid, medium.sound_speed, simMedium.CFL , duration);
dt = kgrid.dt;            % Two of the function's  
t = kgrid.t_array;        % output arguments


%dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1000/dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%round(1/dt) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Source properties

%%%%% REVISAR red = 79; % Red pixels in the image must have this value to work as sources

red256=79; red=[237 28 36]; % Red pixels in the image must have this value to work as sources

if format.BitDepth==8
    source.p_mask = (image == red256);
else
    source.p_mask = (squeeze(image(:,:,1)) == red(1)&...
                    squeeze(image(:,:,2)) == red(2)&...
                    squeeze(image(:,:,3)) == red(3));
end
nSources=sum(sum(source.p_mask));
% Set default values if neccesary
simSource=checkDefault(simSource,t(end));
source.p_mode = simSource.mode;   % User defined input arguments

% Verify the user-defined source type
switch simSource.type
    case 'impulse'
        pVariable = simSource.amplitude*impulseSource(dt, t(end)); 
    case 'nCycles'
        pVariable = simSource.amplitude * nCyclesSource(...
            simSource.f0, simSource.n, dt, t(end));
    case 'whiteNoise'
        pVariable = simSource.amplitude * whiteNoiseSource(...
            simSource.duration, dt, t(end));
    case 'audio'
        if length(simSource.p(1,:))==nSources
            simSource.p=simSource.p';
        end
        %disp(size(simSource.p))
        % set simSource.p row vector
      
        clear pVariable
        tLengthAudio=length(simSource.p)/simSource.fs;
        tAudio=0:1/simSource.fs:tLengthAudio-1/simSource.fs;
        tkWave=0:dt:tLengthAudio - dt;
        if length(simSource.p(:,1))==1 %Only one source function
            pVariable=interp1(tAudio,simSource.p,tkWave);
 
        else  % Multiple source functions
            disp(['1/dt = ' num2str(1/dt)])
            disp(['tLengthAudio = ' num2str(tLengthAudio)])
            disp([size(tAudio);size(simSource.p);size(tkWave)])
            for i=1:nSources 
                pVariable(i,:)=interp1(tAudio,simSource.p(i,:),tkWave);
            end
        end

                                %{
                                %% En revisión (proceso a descartar con resample)       
                                disp(['simSource.fs =' num2str(simSource.fs)])
                                disp(['round(1/dt) =' num2str(round(1/dt))])
                                if length(simSource.p(:,1))==1 %Only one source function

                                    % Some simSource.fs conversion are not allowed by resample()
                                    pVariable=resample(simSource.p,round(1/dt)...
                                            ,simSource.fs);
                                else  
                                    for i=1:nSources
                                        pVariable(i,:)=resample(simSource.p(i,:),round(1/dt)...
                                            ,simSource.fs);
                                    end
                                end
                                disp(size(pVariable))
                                %%
                                %}
    otherwise
        eval(['pVariable =' simSource.type ';']);
end

% Assign the equation to the sources 
%States cut freq for adequate simulation
%REVIRAS XXXX 0.15/DT VÁLIDO PARA HOMOGÉNEO
cMin=min(simMedium.c0,simMedium.c0*simMedium.speedRatio);
fNyquist=cMin/2/dx;
if not(isfield(simSource,'fCut')); simSource.fCut=fNyquist*75/86 ;end 
% For cMin=344 and dx=0.01; fNyquist=17200 Hz; fCut=15kHz
if isfield(simSource,'fCorte'); ...
        simSource.fCut=simSource.fCorte;end % Se mantiene por compatibilidad
if AbsorptionExist 
    %normfCut
    %simSource.fCut
    disp('--- Absorption surfaces detected ---');
    disp(['  z0 air  = ' num2str(z0Air)])
    disp(['  z0 wall = ' num2str(z0Wall)])
end
normfCut=simSource.fCut; 
%REVISAR: si se utiliza exactamente la fNyquist como fCut, entonces en el
%ruido blanco o en el impulso existirá señal por encima de dicha frecuencia
%que provocrá alias y quizás inestabilidad en las soluciones.
if not(isfield(simSource,'order')); simSource.order=15;end
%{
if normfCut<20000
%}    
        [b,a]=butter(simSource.order,normfCut*dt*2); % fNyquist
        % disp('max. frequency < 20 kHz');
        disp(['fNyquist = ' num2str(normfCut) 'Hz'])
%{
else


        normfCut=20000;
        [b,a]=butter(simSource.order,20000*dt*2); % fNyquist=20 kHz
        disp('max. frequency > 20 kHz');
        disp('fNyquist = 20 KHz')

end
%}
filteredVar=0*pVariable; 
if length(filteredVar(:,1))==1 %same signal all sources
    filteredVar=filter(b,a,pVariable);
else
    for i=1:nSources
        filteredVar(i,:)=filter(b,a,pVariable(i,:));
    end
end
pVariable=filteredVar;
source.p=pVariable;

%{
% Clear unneeded fields from the source struct
[simSource.f0, simSource.amplitude, simSource.n] = deal(1);
simSource = rmfield(simSource, {'amplitude'; 'mode' ; 'type'; 'f0'; 'n'});
%}

%% Perfectly reflective surfaces

%%%%%%REVISAR black = 0; % Black pixels must have this value to work as surfaces

black256=0; black=[0 0 0]; % Black pixels must have this value to work as surfaces

if format.BitDepth==8
    source.u_mask = (image == black256);
else
    source.u_mask = (squeeze(image(:,:,1)) == black(1)&...
                    squeeze(image(:,:,2)) == black(2)&...
                    squeeze(image(:,:,3)) == black(3));
end

source.u_mask(1,1) = 1; % Ensuring at least one black pixel
source.ux = 0*kgrid.t_array;
source.uy = 0*kgrid.t_array;
source.u_mode = 'dirichlet';

%% Sensors 

%%%%%%REVISAR green = 113; % Green pixels must have this value to work as sensors

green256=113; green=[34 177 76];% Green pixels must have this value to work as sensors

if format.BitDepth==8
    sensor.mask = (image == green256);
else
    sensor.mask = (squeeze(image(:,:,1)) == green(1)&...
                    squeeze(image(:,:,2)) == green(2)&...
                    squeeze(image(:,:,3)) == green(3));
end


%% Simulation implementation
%{
if format.BitDepth==8
    fullMask = (image < 255);
else
    fullMask = (squeeze(image(:,:,1)) == 255&...
                    squeeze(image(:,:,2)) == 255&...
                    squeeze(image(:,:,3)) == 255);
end
%}
gray256=164; gray=[127 127 127];

if format.BitDepth==8
    legendMask = (image == gray256);
else
    legendMask = (squeeze(image(:,:,1)) == gray(1)&...
                    squeeze(image(:,:,2)) == gray(2)&...
                    squeeze(image(:,:,3)) == gray(3));
end

fullMask=sensor.mask.*simMedium.showSensorMask ...
    +absorptiveMask.*pointMask.*simMedium.showAbsorptionMask ...
    +source.u_mask.*simMedium.showWallMask ...
    +legendMask.*simMedium.showLegendMask ...
    +source.p_mask.*simMedium.showSourceMask ;

% fullMask = image < 255; % This argument shows the image while simulating
%{1
%inputArgs = {...
%'DisplayMask',fullMask, 'RecordMovie',recordVideo, 'LogScale',false};
inputArgs = {...
'DisplayMask',fullMask, 'RecordMovie',recordVideo, 'LogScale',false};

% Run simulation
sensorData = kspaceFirstOrder2D( ...
kgrid, medium, source, sensor, inputArgs{:});
%}

%[b,a]=butter(simSource.order,normfCut*1.2*dt*2); % fNyquist=20 kHz
%nSensors=length(sensorData(:,1));
%{
for i=1:nSensors
    sensorData(i,:)=filter(b,a,sensorData(i,:));
end
%}
%%%% BORRAR LO SIGUIENTE, SÓLO CHEQUEO Y SALTEO LA SIMULACIÓN
%{
sensorData=0;
figure(1);imagesc(source.p_mask)
figure(2);imagesc(source.u_mask)
figure(3);imagesc(sensor.mask)
figure(4);imagesc(absorptive)
figure(5);imagesc(fullMask)
%}

end


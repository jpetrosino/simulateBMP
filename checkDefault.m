function [simSource] = checkDefault(simSource,duration)
%[simSource] = checkDefault(simSource,duration)
%  Check fields of simSource structure to set necesary defaults
% variable 'duration' is total time length of simulation.
% Is used to set whiteNoiseDuration to 1/5 

if not(isfield(simSource,'mode'))
    simSource.mode='additive';
else
    if not(strcmp(simSource.mode,'dirichlet'))
        simSource.mode='additive';
        disp('additive mode set')    
    end
end
if not(isfield(simSource,'amplitude')); simSource.amplitude=10;end
switch simSource.type
    case 'audio'
        if not(isfield(simSource,'fs')); simSource.fs=44100; end
    case 'nCycles'
        if not(isfield(simSource,'f0')); simSource.f0=1000; end
        if not(isfield(simSource,'n')); simSource.n=3; end
    case 'whiteNoise'
        if not(isfield(simSource,'duration')); simSource.duration=duration/5; end
    otherwise
end


function [source] = nCyclesSource(f0,n,dt,duration)
%%                             nCyclesSource                             
%
% [source] = nCyclesSource(f0,n,dt,duration)
%
% Generates a sine wave of f0 frequency, n cycles and the specified
% duration at the source location.
% Completes the vector with zeros after all the cycles were emitted.
%
% If the requested cycle quantity is greater than the duration, it outputs
% a warning.

t0 = 1/f0; 
tCycles = (0:dt:n*t0-dt);
source = zeros(1,round(duration/dt)+1);
source(1:length(tCycles+1)) = sin(2*pi*f0*tCycles);

if tCycles(end) > duration
    disp(' ')
    disp('Source cycles exceed simulation time');
    disp(' ')
end
end


function [source] = whiteNoiseSource(noiseDuration, dt, duration)
%%                             whiteNoiseSource                           
%
% [source] = whiteNoiseSource(noiseDuration, dt, duration)
%
% Generates white noise of the specified length at the source location, 
% with an amplitude value of 1.

%source = zeros(1,round(duration/dt)+1);
source = zeros(1,round(duration/dt));
source(1,1:round(noiseDuration/dt))=2*(rand(1,round(noiseDuration/dt))-0.5);

end

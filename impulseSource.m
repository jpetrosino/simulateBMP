function [ source ] = impulseSource(dt,duration)
%%                             impulseSource 
%
% [source] = impulseSource(dt,duration)
%
% Generates an impulse at the source location

t = 0:dt:duration;
source = zeros(1,round(duration/dt)+1);
source(1) = 1;

end


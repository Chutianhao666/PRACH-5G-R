function [fd,rangeM,propagationDelay] = railwayDopplerModel( ...
    x0,v,d,h,fc,t)
% railwayDopplerModel  Straight-track railway LOS geometry.
% x0: initial longitudinal coordinate (m)
% v: signed velocity (m/s)
% d: lateral BS-track distance (m)
% h: BS/UE antenna height difference (m)
% fc: carrier frequency (Hz)
% t: time vector (s)
% fd: instantaneous Doppler frequency (Hz)
% rangeM: LOS distance (m)
% propagationDelay: rangeM/c (s)

c = 299792458;
assert(isscalar(x0) && isscalar(v) && isscalar(d) && ...
    isscalar(h) && isscalar(fc) && d >= 0 && fc > 0, ...
    'Invalid railway geometry inputs.');
x = x0 + v.*t;
rangeM = sqrt(x.^2 + d.^2 + h.^2);
assert(all(rangeM(:) > 0),'BS and UE must not occupy identical coordinates.');
fd = -(fc/c) .* (v.*x./rangeM);
propagationDelay = rangeM/c;
end

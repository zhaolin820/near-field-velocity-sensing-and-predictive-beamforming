function [a] = array_response(r, theta, M, d, lambda)
%near-field array response vector
%   [a] = array_response(r, theta, M, d, lambda)
%Inputs:
%   r: distance of the target
%   theta: direction of the target
%   M: number of antennas at the BS
%   d: antenna spacing at the BS
%   lambda: signal wavelength
%Outputs:
%   a: array response vector
%Date: 29/12/2023
%Author: Zhaolin Wang

delta_m = (-(M-1)/2 : (M-1)/2)' * d;

% propagation distance
r_m = sqrt(r^2 + delta_m.^2 - 2*r*delta_m*cos(theta));

% array response vector
a = exp( -1i * 2 * pi /lambda * r_m ) ./ r_m;

end
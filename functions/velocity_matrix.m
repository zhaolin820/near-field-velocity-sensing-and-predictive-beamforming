function [d] = velocity_matrix(r, theta, M, d, lambda, vr, vt, Ts, N)
%calculating the near-field doppler-frequency vector at time index n
%   [d] = velocity_matrix(r, theta, M, d, lambda, vr, vt, Ts, N)
%Inputs:
%   r: distance of the target
%   theta: direction of the target
%   M: number of antennas at the BS
%   d: antenna spacing at the BS
%   lambda: signal wavelength
%   vr: radial velocity of the target
%   vt: transverse velocity of the target
%   T_s: symbol period
%   N: length of one CPI
%Outputs:
%   d: doppler-frequenct matrix for one CPI
%Date: 29/12/2023
%Author: Zhaolin Wang

delta_m = (-(M-1)/2 : (M-1)/2)' * d;

r_m = sqrt(r^2 + delta_m.^2 - 2*r*delta_m*cos(theta));

vr_m = (r - delta_m*cos(theta))./r_m * vr;
vt_m = delta_m*sin(theta)./r_m * vt;

v_m = vr_m + vt_m;

N = 1:N;
% Doppler frequency matrix for one CPI
d = exp( -1i * 2 * pi /lambda * v_m*N*Ts );

end



function [X] = match_filter(r, theta, vr, vt, M, N, A, d, lambda, Ts, s)
%matched filter for a specific location and a specific velocity
%   [X] = match_filter(r, theta, vr, vt, M, N, A, d, lambda, Ts, s)
%Inputs:
%   r: distance of the target
%   theta: direction of the target
%   vr: radial velocity of the target
%   vt: transverse velocity of the target
%   M: number of antennas at the BS
%   N: length of one CPI
%   A: array response matrix
%   d: antenna spacing at the BS
%   lambda: signal wavelength
%   T_s: symbol period
%   s: transmit signal
%Outputs:
%   X: matched filter
%Date: 29/12/2023
%Author: Zhaolin Wang


X = zeros(M, N);
for n = 1:N
    s_n = s(:,n);
    [d_n] = velocity_vector(r, theta, M, d, lambda, vr, vt, Ts,n);
    D_n = d_n*d_n.';
    
    X(:,n) = (A.*D_n)*s_n;
end


end


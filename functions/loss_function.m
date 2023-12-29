function [f, g] = loss_function(eta, r, theta, Y, M, N, A, d, lambda, Ts, s)
%Loss function and gradient for maximum likelihood estimation
%   loss_function(eta, r, theta, Y, M, N, A, d, lambda, Ts, s)
%Inputs:
%   eta: vector of velocities
%   r: distance of the target
%   theta: direction of the target
%   Y: signal received over one CPI
%   M: number of antennas at the BS
%   N: length of one CPI
%   A: array response matrix
%   d: antenna spacing at the BS
%   lambda: signal wavelength
%   T_s: symbol period
%   s: transmit signal
%Outputs:
%   f: loss function
%   g: gradient
%Date: 29/12/2023
%Author: Zhaolin Wang

vr = eta(1); vt = eta(2);

% objective function
X = match_filter(r, theta, vr, vt, M, N, A, d, lambda, Ts,s);
f = -2*real(trace(X*Y')) + trace(X*X');

% gradient of objective function
if nargout >= 2 
    g = zeros(2,1);
    delta_m = (-(M-1)/2 : (M-1)/2)' * d;
    r_m = sqrt(r^2 + delta_m.^2 - 2*r*delta_m*cos(theta));
    a = beamfocusing(r, theta, M, d, lambda);
    A = a*a.';
    
    g_X = -conj(Y) + conj(X);
    
    X_vr = zeros(M, N); X_vt = zeros(M, N);
    for n = 1:N
        dn = velocity_vector(r, theta, M, d, lambda, vr, vt, Ts,n);
        qm = (r - delta_m*cos(theta))./r_m;
        pm = delta_m*sin(theta)./r_m;
    
        d_vr = -1i*2*pi/lambda*n*Ts*qm .* dn;
        d_vt = -1i*2*pi/lambda*n*Ts*pm .* dn;
        
        
        X_vr(:,n) = ( 2*A .* ( d_vr*dn.') ) * s(:,n);
        X_vt(:,n) = ( 2*A .* ( d_vt*dn.') ) * s(:,n);
    end

    g(1) = 2*real(trace(g_X.' * X_vr ));
    g(2) = 2*real(trace(g_X.' * X_vt ));
end

function [f, g] = loss_function(Y, r, theta, v, M, N, d, lambda, Ts, s)
%Loss function and gradient for maximum likelihood estimation
%   loss_function(eta, r, theta, Y, M, N, A, d, lambda, Ts, s)
%Inputs:
%   Y: signal received over one CPI
%   r: distance of the target
%   theta: direction of the target
%   v: vector of velocities
%   M: number of antennas at the BS
%   N: length of one CPI
%   d: antenna spacing at the BS
%   lambda: signal wavelength
%   Ts: symbol period
%   s: transmit signal
%Outputs:
%   f: loss function
%   g: gradient
%Date: 03/06/2024
%Author: Zhaolin Wang


vr = v(1); vt = v(2);
delta_m = (-(M-1)/2 : (M-1)/2)' * d;
r_m = sqrt(r^2 + delta_m.^2 - 2*r*delta_m*cos(theta));
X = match_filter(r, theta, vr, vt, M, N, d, lambda, Ts, s);

f = -abs(trace(X'*Y))^2/real(trace(X*X'));   

if nargout > 1 % gradient required
    n = 1:N;
    [d_n] = velocity_vector(r, theta, M, d, lambda, vr, vt, Ts,n);
    qm = (r - delta_m*cos(theta))./r_m;
    pm = delta_m*sin(theta)./r_m;
    
    d_vr = -1i*2*pi/lambda*Ts*qm*n .* d_n;
    d_vt = -1i*2*pi/lambda*Ts*pm*n .* d_n;
    a = array_response(r, theta, M, d, lambda);
    
    H = a.*d_n;
    H_vr = a.*d_vr;
    H_vt = a.*d_vt;
    
    X_vr = zeros(M, N); X_vt = zeros(M, N);
    for n = 1:N
        
        hn = H(:,n);
        hn_vr = H_vr(:,n);
        hn_vt = H_vt(:,n);
        X_vr(:,n) = ( hn_vr*hn.' + hn*hn_vr.' ) * s(:,n);
        X_vt(:,n) = ( hn_vt*hn.' + hn*hn_vt.' ) * s(:,n);
    end
    
    
    X_norm = norm(X,"fro")^2;
    Theta = trace(Y*X')*X_norm;
    Omega = abs(trace(Y*X'))^2;
    g_X = (Theta*Y' - Omega*X')/X_norm^2;
    
    g = zeros(2,1);
    g(1) = -2*real(trace(g_X*X_vr));
    g(2) = -2*real(trace(g_X*X_vt));
end

end

clc
clear all
close all



% Estimating the radial and trasnverse velocities based on maximum
% likelihood estimation and quasi-Newton method

addpath("functions\");
%% system parameter
Delta = 2e-3; % coherence time (s);
c = 3e8; % speed of light (m/s)
f = 28e9; % carrier frequency (Hz)
lambda = c/f; % signal wavelength (m)
M = 512; % number of antennas
d = lambda/2; % antenna spacing (m)
B = 100e3; % system bandwidth (Hz)
Ts = 1/B; % symbol duration (s)
D = M*d; % array apearture (m)
N = floor(Delta/Ts); % coherence period interval
Pt = 10^(10/10); % transmit power (mW)
N0 = 10^(-174/10)*B; % noise power (mW)
SNR = Pt*lambda^2/((4*pi)^3*N0); % signal-to-noise ratio

%% target parameter
r = 20; % distance (m)
theta = 90/180*pi; % direction (rad)
vr = 10; % radial velocity (m/s)
vt = 8; % transverse velocity (m/s)

%% signal model
a = beamfocusing(r, theta, M, d, lambda);
w = conj(a)/sqrt(M); % beamformer
s = sqrt(SNR)*w/sqrt(2) *(randn(1,N) + 1i * randn(1,N)); % transmit signal

% receive signal
A = a*a.';
Y = zeros(M, N);
for n = 1:N
    s_n = s(:,n);
    d_n = velocity_vector(r, theta, M, d, lambda, vr, vt, Ts,n);
    D_n = d_n*d_n.';
    
    z_n = 1/sqrt(2) * (randn(M,1) + 1i * randn(M,1));
    Y(:,n) = (A.*D_n)*s_n + z_n;
end

%% velocity sensing, maximum likelihood estimation

% initialization
vr_init = 9;
vt_init = 9;
eta_init = [vr_init; vt_init];

% quasi-Newton method
options = optimoptions(@fminunc, 'Display','iter-detailed', 'HessianApproximation', 'lbfgs', 'SpecifyObjectiveGradient',true);
eta = fminunc(@(eta)(loss_function(eta, r, theta, Y, M, N, A, d, lambda, Ts,s)),eta_init,options);

disp(['Ground truth: radial velocity - ' num2str(vr) ' m/s, transverse velocity - ' num2str(vt) ' m/s']);
disp(['Estimated results: radial velocity - ' num2str(eta(1)) ' m/s, transverse velocity - ' num2str(eta(2)) ' m/s']);





clc
clear all
close all

% Estimating the radial and trasnverse velocities based on maximum
% likelihood estimation and quasi-Newton method

addpath("functions/");
%% system parameter
Delta = 2e-3; % coherence time (s);
c = 3e8; % speed of light (m/s)
f = 28e9; % carrier frequency (Hz)
lambda = c/f; % signal wavelength (m)
M = 256; % number of antennas
d = lambda/2; % antenna spacing (m)
B = 100e3; % system bandwidth (Hz)
Ts = 1/B; % symbol duration (s)
D = M*d; % array apearture (m)
N = floor(Delta/Ts); % coherence period interval
Pt = 10^(-10/10); % transmit power (mW)
N0 = 10^(-174/10)*B; % noise power (mW)
SNR = Pt*lambda^2/((4*pi)^3*N0); % signal-to-noise ratio


%% target parameter
r = 20; % distance (m)
theta = 60/180*pi; % direction (rad)
vr = 10; % radial velocity (m/s)
vt = 8; % transverse velocity (m/s)

%% signal model
a = array_response(r, theta, M, d, lambda);
w = conj(a); % beamformer
w = sqrt(SNR)*w./norm(w);
s = w/sqrt(2) *(randn(1,N) + 1i * randn(1,N));
Z = 1/sqrt(2) * (randn(M,N) + 1i * randn(M,N));


[X] = match_filter(r, theta, vr, vt, M, N, d, lambda, Ts,s);
beta = randn(1) + 1i*randn(1); beta = beta ./ abs(beta);
Y = beta*X + Z;

%% MLE estimation
v0 = [8; 6];
options = optimoptions(@fminunc, 'Display','iter-detailed', HessianApproximation="lbfgs", SpecifyObjectiveGradient=true);
v = fminunc(@(eta)(loss_function(Y, r, theta, eta, M, N, d, lambda, Ts, s)),v0,options);

disp(['Ground truth: radial velocity - ' num2str(vr) ' m/s, transverse velocity - ' num2str(vt) ' m/s']);
disp(['Estimated results: radial velocity - ' num2str(v(1)) ' m/s, transverse velocity - ' num2str(v(2)) ' m/s']);





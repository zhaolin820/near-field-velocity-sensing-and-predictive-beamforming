clc
clear all
close all

% plotting Fig. 2

addpath("functions\");
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
r = 10; % distance (m)
theta = 90/180*pi; % direction (rad)
vr = 10; % radial velocity (m/s)
vt = 8; % transverse velocity (m/s)

%% signal model
a = array_response(r, theta, M, d, lambda);
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

%% plot Fig. 2
vt_all = -20:0.5:20;
t_velocity = zeros(length(vt_all),1);
for i = 1:length(vt_all)
    t_velocity(i) = -loss_function(Y, r, theta, [vr, vt_all(i)], M, N, d, lambda, Ts, s);
end
t_velocity = t_velocity ./ max(t_velocity);
figure; box on; hold on;
plot(vt_all, t_velocity, 'LineWidth', 1.5);
y = 0:0.2:1;
plot(vt*ones(1, length(y)), y, '--m', 'LineWidth', 1);
set(gca,'linewidth',1);
xlabel('Transverse velocity (m/s)', 'Interpreter','latex');
ylabel('Normalized $g(\mathbf{Y}, \mbox{\boldmath $\eta$}, \mbox{\boldmath $v$})$', 'Interpreter','latex');

vr_all = -20:0.5:20;
r_velocity = zeros(length(vr_all),1);
for i = 1:length(vr_all)
    r_velocity(i) = -loss_function(Y, r, theta, [vr_all(i), vt], M, N, d, lambda, Ts, s);
end
r_velocity = r_velocity ./ max(r_velocity);
figure; box on; hold on;
plot(vr_all, r_velocity, 'LineWidth', 1.5);
y = 0:0.2:1;
plot(vr*ones(1, length(y)), y, '--m', 'LineWidth', 1);
set(gca,'linewidth',1);
xlabel('Radial velocity (m/s)', 'Interpreter','latex');
ylabel('Normalized $g(\mathbf{Y}, \mbox{\boldmath $\eta$}, \mbox{\boldmath $v$})$', 'Interpreter','latex');

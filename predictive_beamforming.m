clc
clear all
close all


% Predictive beamforming based on the near-field veloity sensing.
% Runing this script can take a long time.
% plot Fig. 3, Fig. 4, and Fig. 5.

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

%% user motion model
t_total = 20e-3; % total time (s)
s = 20; % speed of user 20 m/s
R = 2*s/pi;
rho = 10 + R;
s_rad = s/R;
t = 0:Delta:t_total; % time 
alpha = s_rad*t;

r = sqrt( R^2 + rho^2 - 2*R*rho*sin(alpha) );
theta = acos(R./r.*cos(alpha));

r_diff = [r(2:end), 10000] - r; 
vr = r_diff/Delta;

theta_diff = [theta(2:end), 10000] - theta;
vt = theta_diff/Delta .* r;

%% predictive beamforming
r_predicted_all = zeros(1, length(t)); theta_predicted_all = zeros(1, length(t));
vr_estimated_all = zeros(1, length(t)); vt_estimated_all = zeros(1, length(t));
R_all = zeros(length(t), 1);
R_optimal_all = zeros(1, length(t));
options = optimoptions(@fminunc, 'Display','iter-detailed', HessianApproximation="lbfgs", SpecifyObjectiveGradient=true);
r_predicted = r(2); theta_predicted = theta(2); vr_estimated = vr(1); vt_estimated = vt(1);

hw = waitbar(0,'Running Predictive Beamforming ...');
for l = 2:(length(t)-1)
    disp(['******************** CPI - ' num2str(l) ' ********************']);
    r_predicted_all(l) = r_predicted;
    theta_predicted_all(l) = theta_predicted;

    % predicted beamforming
    D_predicted = velocity_vector(r_predicted, theta_predicted, M, d, lambda, vr_estimated, vt_estimated, Ts, 1:N);
    a_predicted = array_response(r_predicted, theta_predicted, M, d, lambda);
    w_predicted = diag(conj(a_predicted))*conj(D_predicted);
    w_predicted = w_predicted*sqrt(SNR)./vecnorm(w_predicted);
    
    % optimal beamforming
    D_optimal= velocity_vector(r(l), theta(l), M, d, lambda, vr(l), vt(l), Ts, 1:N);
    a_optimal = array_response(r(l), theta(l), M, d, lambda);
    w_optimal = diag(conj(a_optimal))*conj(D_optimal);
    w_optimal = w_optimal*sqrt(SNR)./vecnorm(w_optimal);
    
    % transmit signal
    s = 1/sqrt(2) * w_predicted .* (randn(1,N) + 1i * randn(1,N));


    % echo signal and communication rate
    R = 0; R_optimal = 0;
    Y = zeros(M, N);
    a = array_response(r(l), theta(l), M, d, lambda); A = a*a.';  
    D = velocity_vector(r(l), theta(l), M, d, lambda, vr(l), vt(l), Ts, 1:N);
    for n = 1:N
        d_n = D(:,n);
        
        % echo signal
        s_n = s(:,n);
        D_n = d_n*d_n.'; 
        z_n = 1/sqrt(2) * (randn(M,1) + 1i * randn(M,1));
        Y(:,n) = (A.*D_n)*s_n + z_n;

        % communication rate with predicted beamforming
        h_n = a.'*diag(d_n);
        R = R + 1/N*log2(1 + abs(h_n*w_predicted(:,n))^2 );

        % communication rate with optimal beamforming
        R_optimal = R_optimal + 1/N*log2(1 + abs(h_n*w_optimal(:,n))^2 );
    end

    % velocity sensing
    v0 = [vr_estimated; vt_estimated];
    v = fminunc(@(eta)(loss_function(Y, r_predicted, theta_predicted, eta, M, N, d, lambda, Ts, s)),v0,options);
    
    vr_estimated = v(1); 
    vt_estimated = v(2);
    
    % user location prediction
    theta_predicted = theta_predicted + vt_estimated/r_predicted * Delta;
    r_predicted = r_predicted + vr_estimated*Delta;

    

    vr_estimated_all(l) = vr_estimated;
    vt_estimated_all(l) = vt_estimated;
    R_all(l) = R; R_optimal_all(l) = R_optimal;

    waitbar((l-1)/(length(t)-1),hw);
end
close(hw);

figure; box on; hold on;
plot(t(2:end-1)*1e3, vr(2:end-1), 'LineWidth', 1.5);
plot(t(2:end-1)*1e3, vr_estimated_all(2:end-1), '--', 'LineWidth', 1.5);
ylabel('Radial velocity (m/s)', 'Interpreter','latex');
xlabel('Time (ms)', 'Interpreter','latex');
set(gca,'linewidth',1);
legend('Radial velocity, ground-truth', 'Radial velocity, estimated', 'Interpreter','latex');

figure; box on; hold on;
plot(t(2:end-1)*1e3, vt(2:end-1), 'LineWidth', 1.5);
plot(t(2:end-1)*1e3, vt_estimated_all(2:end-1), '--', 'LineWidth', 1.5);
ylabel('Transverse velocity (m/s)', 'Interpreter','latex');
xlabel('Time (ms)', 'Interpreter','latex');
set(gca,'linewidth',1);
legend('Transverse velocity, ground-truth', 'Transverse velocity, estimated', 'Interpreter','latex');

figure; box on; hold on
[X_real, Y_real] = pol2cart(theta, r);
plot(X_real,Y_real, 'LineWidth', 1.5);
[X_predicted, Y_predicted] = pol2cart(theta_predicted_all, r_predicted_all);
plot(X_predicted(2:end-1),Y_predicted(2:end-1), '--', 'LineWidth', 1.5);
set(gca,'linewidth',1);
xlabel('$x$-axis (m)', 'Interpreter','latex');
ylabel('$y$-axis (m)', 'Interpreter','latex');
legend('Ground-truth trajectory', 'Predicted trajectory', 'Interpreter','latex');

figure; box on; hold on;
plot(t(2:end-1)*1e3, R_optimal_all(2:end-1), 'LineWidth', 1.5);
plot(t(2:end-1)*1e3, R_all(2:end-1), '--', 'LineWidth', 1.5);
ylabel('Achievable rate (bit/s/Hz)', 'Interpreter','latex');
xlabel('Time (ms)', 'Interpreter','latex');
set(gca,'linewidth',1);
legend('Optimal', 'Proposed', 'Interpreter','latex');


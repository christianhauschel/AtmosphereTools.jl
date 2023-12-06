clear all;
close all; clc; 

seed = [23341 23342 23343 23344];

h = 20000.0;
V = 10.0;

DCM = [1 0 0; 0 1 0; 0 0 1];

fs = 1000;

scale_length_med = 533.4;
dt = 1/fs;
wingspan = 1.0;
V_wind_at6m = 15; % 
angle_wind_at6m = 0.0; % degrees clowckwise from north

noise_seeds = [23341 23342 23343 23344];

t = 0:dt:5;
t = t';
n_t = length(t);

h_vec = ones(n_t, 1) * h;
V_vec = ones(n_t, 1) * h;
signal_h = timeseries(h_vec, t);
signal_V = timeseries(V, t);

DCM_vec = zeros(3, 3, n_t);
for i = 1:n_t
    DCM_vec(:,:,i) = DCM(:,:);
end
signal_DCM = timeseries(DCM_vec, t);


res = sim("model_dryden.slx");


V_wind = res.V_wind.Data;
omega_wind = res.omega_wind.Data;
noise = res.noise.Data;

figure()
plot(t, V_wind)
xlabel("t [s]")
ylabel("V [m/s]")

% figure()
% plot(t, omega_wind)
% xlabel("t [s]")
% ylabel("\omega [rad/s]")


table_Vwind = table(t, V_wind(:,1), V_wind(:,2), V_wind(:,3), VariableNames=["t", "u", "v", "w"]);
writetable(table_Vwind, "V.csv")

table_omega = table(t, omega_wind(:,1), omega_wind(:,2), omega_wind(:,3), VariableNames=["t", "p", "q", "r"]);
writetable(table_omega, "omega.csv")

table_noise = table(t, noise(:,1), noise(:,2), noise(:,3), noise(:,4), VariableNames=["t", "u", "v", "w", "r"]);
writetable(table_noise, "noise.csv")

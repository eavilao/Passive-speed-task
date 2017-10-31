%% translation
v_mu = 13.75; % mean forward speed cm/s
d = 1500; % total distance to simulate
nt = round(d/v_mu); % number of time points to simulate d @ v_mu
nsim = 1000; % number of simulations

figure; hold on;

v_sig_vis = 4.9; % visual uncertainty in decoding linear speed
v_sig_ves = 3.1; % vestibular uncertainty in decoding linear speed
for i=1:nsim
    plot(v_mu*(1:nt),cumsum(normrnd(0,v_sig_vis,[1 nt])),'Color',[0.5 1 0.5]);
    plot(v_mu*(1:nt),cumsum(normrnd(0,v_sig_ves,[1 nt])),'r');
end

%% rotation
w_mu = 37.5; % mean rotation speed deg/s
R = 10; % number of rotations to simulate
nt = round(R*360/w_mu); % number of time points to simulate R @ w_mu
nsim = 1000; % number of simulations

figure; hold on;

w_sig_vis = 13.2; % visual uncertainty in decoding angular speed
w_sig_ves = 7.8; % vestibular uncertainty in decoding angular speed
for i=1:nsim
    plot(w_mu*(1:nt)/360,cumsum(normrnd(0,w_sig_vis,[1 nt])/360),'Color',[0.5 1 0.5]);
    plot(w_mu*(1:nt)/360,cumsum(normrnd(0,w_sig_ves,[1 nt]))/360,'r');
end
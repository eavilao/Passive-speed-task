% code to perform multiple linear regression usign least squares for the
% linear and angular speed protocols. % load experiments.m

units = experiments(1).singleunits;
% bi-sensory neurons
resp_ves_exc = experiments(1).populations(1).ves.exc.indx;
resp_vis_exc = experiments(1).populations(1).vis.exc.indx;
ves_vis_exc = resp_ves_exc & resp_vis_exc;

resp_ves_sup = experiments(1).populations(1).ves.sup.indx;
resp_vis_sup = experiments(1).populations(1).vis.sup.indx;
ves_vis_sup = resp_ves_sup & resp_vis_sup;

units = units(ves_vis_exc | ves_vis_sup); 

% linearspeed
cond = {'ves','vis'};
for cellNum = 1:length(units)
    for i=1:length(cond)
        pst = mean(units(cellNum).(cond{i}).rate_pst)';
        speed_pst(:,i,cellNum) = pst(1:251); % time
        speed_mu(:,i,cellNum) = cell2mat({units(cellNum).(cond{i}).rate_avg.mu})'; % 4 speeds
    end
end
cond = {'com'};
for cellNum = 1:length(units)
    pst =  mean(units(cellNum).(cond{1}).rate_pst)';
    speed_pst_com(:,cellNum) = pst(1:251); % time
    speed_mu_com(:,cellNum) = cell2mat({units(cellNum).(cond{1}).rate_avg.mu})'; % 4 speeds
end

%% weighted linear sum model

% over time
for cellNum = 1:length(units)
[reg_coeff_pst(:,cellNum),conf_int_pst(:,:,cellNum)] = regress(speed_pst_com(:,cellNum),speed_pst(:,:,cellNum),0.32);
end

% over avg
for cellNum = 1:length(units)
[reg_coeff_avg(:,cellNum),conf_int_avg(:,:,cellNum)] = regress(speed_mu_com(:,cellNum),speed_mu(:,:,cellNum),0.32);
end

% plot
figure; hold on;
for cellNum = 1:length(units)
    plot(reg_coeff_pst(1,cellNum),reg_coeff_pst(2,cellNum), '.r', 'MarkerSize', 25)
    errorbar(reg_coeff_pst(1,cellNum),reg_coeff_pst(2,cellNum),...
        reg_coeff_pst(1,cellNum)-conf_int_pst(1,1,cellNum),reg_coeff_pst(1,cellNum)-conf_int_pst(1,2,cellNum),...
        reg_coeff_pst(2,cellNum)-conf_int_pst(2,1,cellNum),reg_coeff_pst(2,cellNum)-conf_int_pst(2,2,cellNum), 'k')
end
plot(0:3,0:3, '-r'); axis([-1.5 1.5 -1.5 1.5])
set(gca, 'TickDir','out'); vline(0);hline(0)
box off
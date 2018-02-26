% compute singleunits distance correlation (on psth) between singleunits/multiunits and lfps of the same channel. 

%% Linspeed
%%
sua_units_lin = experiments(1).singleunits; mua_units_lin = experiments(1).multiunits;
sua_nunits_lin = length(sua_units_lin); mua_nunits_lin = length(mua_units_lin);
lfp_units_lin = experiments(1).lfps; lfp_nunits_lin = length(lfp_units_lin);

cond={'ves','vis','com'}; nconds = length(cond);

%% pairwise correlations in psth vs distance singleunits

% extract multiunits
t_samp_sua = sua_units_lin(1).ves.time>0.5 & sua_units_lin(1).ves.time<1.2; 
for i=1:sua_nunits_lin
    sua_session_id(i) = sua_units_lin(i).session_id;
    sua_channel_no(i) = sua_units_lin(i).channel_no;
    for j=1:nconds
        resp_sua_lin_conds(i,j,:) = mean(sua_units_lin(i).(cond{j}).rate_pst(:,t_samp_sua)); % just motion period
    end
%     resp_lfp_lin(i,:) = mean(squeeze(resp_lfp_lin_conds(i,:,:))); 
end

dists = 1:15; ndists = length(dists);
sua_sua_corr_lin = cell(nconds,ndists);

for l=1:nconds
    for k=1:length(dists)
        for i=1:sua_nunits_lin
            for j=1:sua_nunits_lin
                if (sua_session_id(i)==sua_session_id(j)) && (abs(sua_channel_no(i) - sua_channel_no(j))==dists(k))
                    sua_sua_corr_lin{l,k} = [sua_sua_corr_lin{l,k} ; corr(squeeze(resp_sua_lin_conds(i,l,:)),squeeze(resp_sua_lin_conds(j,l,:)))];
                end
            end
        end
        corr_sua_sua_pairs_mu_lin(l,k) = mean(sua_sua_corr_lin{l,k});
        corr_sua_sua_pairs_std_lin(l,k) = std(sua_sua_corr_lin{l,k});
    end
end

distance = 0:0.1:1.5;
figure; hold on;
plot(distance, corr_sua_sua_pairs_mu_lin(1,:)', 'r', 'LineWidth', 2); 
plot(distance, corr_sua_sua_pairs_mu_lin(2,:)', 'g','LineWidth', 2); 
 box off;
set(gca, 'TickDir', 'out', 'ylim',([-0.05 0.2]), 'FontSize',18); 
xlabel('Distance (mm)'); ylabel('Corr coeff.');


%% pairwise correlations in psth vs distance lfps

%extract lfps
t_samp = lfp_units_lin(1).ves.time>0.5 & lfp_units_lin(1).ves.time<1.2; 

for i=1:lfp_nunits_lin
    lfp_session_id(i) = lfp_units_lin(i).session_id;
    lfp_channel_no(i) = lfp_units_lin(i).channel_no;
    for j=1:nconds
        resp_lfp_lin_conds(i,j,:) = mean(lfp_units_lin(i).(cond{j}).wave_pst(:,t_samp)); % just motion period
    end
%     resp_lfp_lin(i,:) = mean(squeeze(resp_lfp_lin_conds(i,:,:))); 
end


dists = 1:15; ndists = length(dists);
lfp_lfp_corr_lin = cell(nconds,ndists);

for l=1:nconds
    for k=1:length(dists)
        for i=1:lfp_nunits_lin
            for j=1:lfp_nunits_lin
                if (lfp_session_id(i)==lfp_session_id(j)) && (abs(lfp_channel_no(i) - lfp_channel_no(j))==dists(k))
                    lfp_lfp_corr_lin{l,k} = [lfp_lfp_corr_lin{l,k} ; corr(squeeze(resp_lfp_lin_conds(i,l,:)),squeeze(resp_lfp_lin_conds(j,l,:)))];
                end
            end
        end
        corr_lfp_lfp_pairs_mu_lin(l,k) = mean(lfp_lfp_corr_lin{l,k});
        corr_lfp_lfp_pairs_std_lin(l,k) = std(lfp_lfp_corr_lin{l,k});
    end
end

distance = 0:0.1:1.5;
figure; hold on;
plot(distance, corr_lfp_lfp_pairs_mu_lin(1,:)', 'r', 'LineWidth', 2); 
plot(distance, corr_lfp_lfp_pairs_mu_lin(2,:)', 'g','LineWidth', 2); 
 box off;
set(gca, 'TickDir', 'out', 'ylim',([0.7 1]), 'ytick', [0.7 0.8 0.9 1], 'FontSize',18)
xlabel('Distance (mm)'); ylabel('Corr coeff.');


%% pairwise correlations in psth vs distance MUA

% extract multiunits
t_samp_mua = mua_units_lin(1).ves.time>0.5 & mua_units_lin(1).ves.time<1.2; 
for i=1:mua_nunits_lin
    mua_session_id(i) = mua_units_lin(i).session_id;
    mua_channel_no(i) = mua_units_lin(i).channel_no;
    for j=1:nconds
        resp_mua_lin_conds(i,j,:) = mean(mua_units_lin(i).(cond{j}).rate_pst(:,t_samp_mua)); % just motion period
    end
%     resp_lfp_lin(i,:) = mean(squeeze(resp_lfp_lin_conds(i,:,:))); 
end

dists = 1:15; ndists = length(dists);
mua_mua_corr_lin = cell(nconds,ndists);

for l=1:nconds
    for k=1:length(dists)
        for i=1:mua_nunits_lin
            for j=1:mua_nunits_lin
                if (mua_session_id(i)==mua_session_id(j)) && (abs(mua_channel_no(i) - mua_channel_no(j))==dists(k))
                    mua_mua_corr_lin{l,k} = [mua_mua_corr_lin{l,k} ; corr(squeeze(resp_mua_lin_conds(i,l,:)),squeeze(resp_mua_lin_conds(j,l,:)))];
                end
            end
        end
        corr_mua_mua_pairs_mu_lin(l,k) = mean(mua_mua_corr_lin{l,k});
        corr_mua_mua_pairs_std_lin(l,k) = std(mua_mua_corr_lin{l,k});
    end
end

distance = 0:0.1:1.5;
figure; hold on;
plot(distance, corr_mua_mua_pairs_mu_lin(1,:)', 'r', 'LineWidth', 2); 
plot(distance, corr_mua_mua_pairs_mu_lin(2,:)', 'g','LineWidth', 2); 
 box off;
set(gca, 'TickDir', 'out', 'ylim',([0.1 0.6]), 'FontSize',18)
xlabel('Distance (mm)'); ylabel('Corr coeff.');


%% Angular speed
%%
sua_units_ang = experiments(2).singleunits; mua_units_ang = experiments(2).multiunits;
sua_nunits_ang = length(sua_units_ang); mua_nunits_ang = length(mua_units_ang);
lfp_units_ang = experiments(2).lfps; lfp_nunits_ang = length(lfp_units_ang);

cond={'ves','vis'}; nconds = length(cond);


%% pairwise correlations in psth vs distance singleunits

% extract singleunits angular speed
t_samp_sua = sua_units_ang(1).ves.time>0.5 & sua_units_ang(1).ves.time<1.2; 
for i=1:sua_nunits_ang
    sua_session_id(i) = sua_units_ang(i).session_id;
    sua_channel_no(i) = sua_units_ang(i).channel_no;
    for j=1:nconds
        resp_sua_ang_conds(i,j,:) = mean(sua_units_ang(i).(cond{j}).rate_pst(:,t_samp_sua)); % just motion period
    end
%     resp_lfp_ang(i,:) = mean(squeeze(resp_lfp_ang_conds(i,:,:))); 
end

dists = 1:15; ndists = length(dists);
sua_sua_corr_ang = cell(nconds,ndists);

for l=1:nconds
    for k=1:length(dists)
        for i=1:sua_nunits_ang
            for j=1:sua_nunits_ang
                if (sua_session_id(i)==sua_session_id(j)) && (abs(sua_channel_no(i) - sua_channel_no(j))==dists(k))
                    sua_sua_corr_ang{l,k} = [sua_sua_corr_ang{l,k} ; corr(squeeze(resp_sua_ang_conds(i,l,:)),squeeze(resp_sua_ang_conds(j,l,:)))];
                end
            end
        end
        corr_sua_sua_pairs_mu_ang(l,k) = mean(sua_sua_corr_ang{l,k});
        corr_sua_sua_pairs_std_ang(l,k) = std(sua_sua_corr_ang{l,k});
    end
end

distance = 0:0.1:1.5;
figure; hold on;
plot(distance, corr_sua_sua_pairs_mu_ang(1,:)', 'r', 'LineWidth', 2); 
plot(distance, corr_sua_sua_pairs_mu_ang(2,:)', 'g','LineWidth', 2); 
 box off;
set(gca, 'TickDir', 'out', 'ylim',([-0.05 0.23]), 'ytick', [0 0.1 0.2], 'FontSize',18); 
xlabel('Distance (mm)'); ylabel('Corr coeff.');


%% pairwise correlations in psth vs distance lfps

%extract lfps
t_samp = lfp_units_ang(1).ves.time>0.5 & lfp_units_ang(1).ves.time<1.2; 

for i=1:lfp_nunits_ang
    lfp_session_id(i) = lfp_units_ang(i).session_id;
    lfp_channel_no(i) = lfp_units_ang(i).channel_no;
    for j=1:nconds
        resp_lfp_ang_conds(i,j,:) = mean(lfp_units_ang(i).(cond{j}).wave_pst(:,t_samp)); % just motion period
    end
%     resp_lfp_ang(i,:) = mean(squeeze(resp_lfp_ang_conds(i,:,:))); 
end


dists = 1:15; ndists = length(dists);
lfp_lfp_corr_ang = cell(nconds,ndists);

for l=1:nconds
    for k=1:length(dists)
        for i=1:lfp_nunits_ang
            for j=1:lfp_nunits_ang
                if (lfp_session_id(i)==lfp_session_id(j)) && (abs(lfp_channel_no(i) - lfp_channel_no(j))==dists(k))
                    lfp_lfp_corr_ang{l,k} = [lfp_lfp_corr_ang{l,k} ; corr(squeeze(resp_lfp_ang_conds(i,l,:)),squeeze(resp_lfp_ang_conds(j,l,:)))];
                end
            end
        end
        corr_lfp_lfp_pairs_mu_ang(l,k) = mean(lfp_lfp_corr_ang{l,k});
        corr_lfp_lfp_pairs_std_ang(l,k) = std(lfp_lfp_corr_ang{l,k});
    end
end

distance = 0.1:0.1:1.5;
figure; hold on;
plot(distance, corr_lfp_lfp_pairs_mu_ang(1,:)', 'r', 'LineWidth', 2); 
plot(distance, corr_lfp_lfp_pairs_mu_ang(2,:)', 'g','LineWidth', 2); 
 box off;
set(gca, 'TickDir', 'out', 'ylim',([0.75 1]), 'ytick', [0.8 0.9 1], 'FontSize',18)
xlabel('Distance (mm)'); ylabel('Corr coeff.');


%% pairwise correlations in psth vs distance MUA

% extract multiunits
t_samp_mua = mua_units_ang(1).ves.time>0.5 & mua_units_ang(1).ves.time<1.2; 
for i=1:mua_nunits_ang
    mua_session_id(i) = mua_units_ang(i).session_id;
    mua_channel_no(i) = mua_units_ang(i).channel_no;
    for j=1:nconds
        resp_mua_ang_conds(i,j,:) = mean(mua_units_ang(i).(cond{j}).rate_pst(:,t_samp_mua)); % just motion period
    end
%     resp_lfp_ang(i,:) = mean(squeeze(resp_lfp_ang_conds(i,:,:))); 
end

dists = 1:15; ndists = length(dists);
mua_mua_corr_ang = cell(nconds,ndists);

for l=1:nconds
    for k=1:length(dists)
        for i=1:mua_nunits_ang
            for j=1:mua_nunits_ang
                if (mua_session_id(i)==mua_session_id(j)) && (abs(mua_channel_no(i) - mua_channel_no(j))==dists(k))
                    mua_mua_corr_ang{l,k} = [mua_mua_corr_ang{l,k} ; corr(squeeze(resp_mua_ang_conds(i,l,:)),squeeze(resp_mua_ang_conds(j,l,:)))];
                end
            end
        end
        corr_mua_mua_pairs_mu_ang(l,k) = mean(mua_mua_corr_ang{l,k});
        corr_mua_mua_pairs_std_ang(l,k) = std(mua_mua_corr_ang{l,k});
    end
end

distance = 0.1:0.1:1.5;
figure; hold on;
plot(distance, corr_mua_mua_pairs_mu_ang(1,:)', 'r', 'LineWidth', 2); 
plot(distance, corr_mua_mua_pairs_mu_ang(2,:)', 'g','LineWidth', 2); 
 box off;
set(gca, 'TickDir', 'out', 'ylim',([0.15 0.6]), 'FontSize',18)
xlabel('Distance (mm)'); ylabel('Corr coeff.');
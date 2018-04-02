% compute singleunits distance correlation (on psth) between singleunits/multiunits and lfps of the same channel. 

sua_units_lin = experiments(1).singleunits; mua_units_lin = experiments(1).multiunits;
sua_nunits_lin = length(sua_units_lin); mua_nunits_lin = length(mua_units_lin);
cond={'ves','vis','com'}; nconds = length(cond);
% extract singleunits
for i=1:sua_nunits_lin
    sua_session_id(i) = sua_units_lin(i).session_id;
    sua_channel_no(i) = sua_units_lin(i).channel_no;
    for j=1:nconds
        resp_sua_lin(i,j,:) = mean(sua_units_lin(i).(cond{j}).rate_pst(:,[91:161])); % just motion period
    end
end
%extract multiunits
for i=1:mua_nunits_lin
    mua_session_id(i) = mua_units_lin(i).session_id;
    mua_channel_no(i) = mua_units_lin(i).channel_no;
    for j=1:nconds
        resp_mua_lin(i,j,:) = mean(mua_units_lin(i).(cond{j}).rate_pst(:,[91:161])); % just motion period
    end
end


%% pairwise correlations in psth vs distance
dists = 0:15; ndists = length(dists);
sua_mua_corr_lin = cell(nconds,ndists);
for l=1:nconds
    for k=1:length(dists)
        for i=1:sua_nunits_lin
            for j=1:mua_nunits_lin
                if sua_session_id(i)==mua_session_id(j) && abs(sua_channel_no(i) - mua_channel_no(j))==dists(k)
                    sua_mua_corr_lin{l,k} = [sua_mua_corr_lin{l,k} ; corr(squeeze(resp_sua_lin(i,l,:)),squeeze(resp_mua_lin(j,l,:)))];
                end
            end
        end
        corr_sua_mua_pairs_mu_lin(l,k) = mean(sua_mua_corr_lin{l,k});
        corr_sua_mua_pairs_std_lin(l,k) = std(sua_mua_corr_lin{l,k});
        
    end
end

%% Compute for angular speed
sua_units_ang = experiments(2).singleunits; mua_units_ang = experiments(2).multiunits;
sua_nunits_ang = length(sua_units_ang); mua_nunits_ang = length(mua_units_ang);
cond={'ves','vis'}; nconds = length(cond);
% extract singleunits
for i=1:sua_nunits_ang
    sua_session_id(i) = sua_units_ang(i).session_id;
    sua_channel_no(i) = sua_units_ang(i).channel_no;
    for j=1:nconds
        resp_max_sua_ang(i,j,:) = mean(sua_units_ang(i).(cond{j}).rate_pst(:,[91:161]));
    end
end
%extract multiunits
for i=1:mua_nunits_ang
    mua_session_id(i) = mua_units_ang(i).session_id;
    mua_channel_no(i) = mua_units_ang(i).channel_no;
    for j=1:nconds
        resp_max_mua_ang(i,j,:) = mean(mua_units_ang(i).(cond{j}).rate_pst(:,[91:161]));
    end
end


%% pairwise correlations in psth vs distance
dists = 0:15; ndists = length(dists);
sua_mua_corr_ang = cell(nconds,ndists);
for l=1:nconds
    for k=1:length(dists)
        for i=1:sua_nunits_ang
            for j=1:mua_nunits_ang
                if sua_session_id(i)==mua_session_id(j) && abs(sua_channel_no(i) - mua_channel_no(j))==dists(k)
                    sua_mua_corr_ang{l,k} = [sua_mua_corr_ang{l,k} ; corr(squeeze(resp_max_sua_ang(i,l,:)),squeeze(resp_max_mua_ang(j,l,:)))];
                end
            end
        end
        corr_sua_mua_pairs_mu_ang(l,k) = mean(sua_mua_corr_ang{l,k});
        corr_sua_mua_pairs_std_ang(l,k) = std(sua_mua_corr_ang{l,k});
        
    end
end

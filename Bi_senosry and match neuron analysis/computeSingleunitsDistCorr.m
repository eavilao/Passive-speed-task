% compute singleuntis distance correlation between waveforms

sua_units = experiments(1).singleunits; 
nunits = length(sua_units);
cond={'ves','vis','com'}; nconds = length(cond);
% extract singleunits
for i=1:nunits
    session_id(i) = units(i).session_id;
    channel_no(i) = units(i).channel_no;
    for j=1:nconds
         resp_max(i,j) = max(mean(units(i).(cond{j}).rate_pst)); % to compute the max reponse correlation
    end
end


%% pairwise correlations in peak response vs distance
dists = 1:15; ndists = length(dists);
su_pairs = cell(nconds,ndists);
for l=1:nconds
    for k=1:length(dists)
        for i=1:nunits
            for j=1:nunits
                if session_id(i)==session_id(j) && abs(channel_no(i) - channel_no(j))==dists(k)
                    su_pairs{l,k} = [su_pairs{l,k} ; [resp_max(i,l) resp_max(j,l)]];
                end
            end
        end
        corr_respmax_mu(l,k) = corr(su_pairs{l,k}(:,1),su_pairs{l,k}(:,2));
        corr_respmax_std(l,k) = sqrt(((1 - corr_respmax_mu(k))^2)/size(su_pairs{l,k},1));
    end
end

%% plot
interDist =[100:100:1500];
figure; hold on;
shadedErrorBar(interDist,corr_respmax_mu(1,:),corr_respmax_std(1,:),'r');
shadedErrorBar(interDist,corr_respmax_mu(2,:),corr_respmax_std(2,:),'g');
shadedErrorBar(interDist,corr_respmax_mu(3,:),corr_respmax_std(3,:),'b');
box off;
set(gca, 'TickDir', 'out', 'ylim',([-0.5 0.5]), 'ytick', [-0.5 0 0.5], 'FontSize',18)

% compute bi-sensory neurons psth correlation 
% TODO integrate with main code. 

units = experiments(1).singleunits;
nunits = length(units);
cond={'ves','vis','com'}; nconds = length(cond);

% bi-sensory neurons
resp_ves_exc = experiments(1).populations(1).ves.exc.indx;
resp_vis_exc = experiments(1).populations(1).vis.exc.indx;
ves_vis_exc = resp_ves_exc & resp_vis_exc;

resp_ves_sup = experiments(1).populations(1).ves.sup.indx;
resp_vis_sup = experiments(1).populations(1).vis.sup.indx;
ves_vis_sup = resp_ves_sup & resp_vis_sup;

units = units(ves_vis_exc | ves_vis_sup); 
t_samp =  units(1).ves.time>0 & units(1).ves.time<1.5;

for i=1:length(units)
    session_id(i) = units(i).session_id;
    channel_no(i) = units(i).channel_no;
    for j=1:nconds
        resp_psth(i,j,:) = mean(units(i).(cond{j}).rate_pst(:,t_samp)); %just motion period %mean(units(i).(cond{j}).rate_pst(:,[1:251])); 
    end
end

%% correlations in psth vs condition
for k=1:length(units)
    corr_ves_vis(k) = corr(squeeze(resp_psth(k,1,:)),squeeze(resp_psth(k,2,:)));
    std_ves_vis(k) = sqrt(((1 - corr_ves_vis(k))^2)/nunits);
    
    corr_ves_com(k) = corr(squeeze(resp_psth(k,1,:)),squeeze(resp_psth(k,3,:)));
    std_ves_com(k) = sqrt(((1 - corr_ves_com(k))^2)/nunits);
    
    corr_vis_com(k) = corr(squeeze(resp_psth(k,2,:)),squeeze(resp_psth(k,3,:)));
    std_ves_com(k) = sqrt(((1 - corr_vis_com(k))^2)/nunits);
end

for i=1:length(units)
    if find(units(i).ves.stats.flags.exc == 1);
        indx_all(:,i) = 1;
    else
        indx_all(:,i) = 0;
    end
end
indx_all = logical(indx_all);

%% plot
figure; hold on;
scatter(corr_vis_com(indx_all),corr_ves_com(indx_all),'ok', 'filled');
scatter(corr_vis_com(~indx_all),corr_ves_com(~indx_all),'ok');
plot(-1:1,-1:1,'k');
box off; axis equal;
set(gca, 'TickDir', 'out', 'ylim',([0 1]),'xlim',([0 1]), 'ytick', [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'xtick', [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], 'FontSize',18)
xlabel('Corr (r_{vis},r_{com})', 'FontSize',18); ylabel('Corr (r_{ves},r_{com})','FontSize',20); 

% x axis histogram
figure; hold on;
hist(corr_vis_com(indx_all | ~indx_all),25);
xlim ([0 1])
set(gca,'ylim', [0 7], 'TickDir', 'out','ytick', [0 7],'xtick',[0 0.5 1], 'Fontsize', 18);
xlabel('Corr (r_{vis},r_{com})','Fontsize', 18);
vline(mean(corr_vis_com(indx_all | ~indx_all)));

% y axis histogram
figure; hold on;
hist(corr_ves_com(indx_all | ~indx_all),25);
xlim ([0 1])
vline(mean(corr_ves_com(indx_all | ~indx_all)));
set(gca,'ylim', [0 6], 'TickDir', 'out','ytick', [0 6],'xtick',[0 0.5 1], 'Xdir','reverse','Fontsize', 18); xlabel('Corr (r_{ves},r_{com})', 'Fontsize', 18);

% diagonal histogram
figure; hold on; 
hist(corr_vis_com(indx_all | ~indx_all)-corr_ves_com(indx_all | ~indx_all),25);
%xlim ([0 1])
set(gca,'ylim', [0 6],'xlim',[-0.707 0.707] ,'TickDir', 'out','ytick', [0 6],'Fontsize', 18); xlabel('diagonal', 'Fontsize', 18); % 1/sqrt(2) 0.702 
vline(mean(corr_vis_com(indx_all | ~indx_all)-corr_ves_com(indx_all | ~indx_all)));
vline(0);

%[p,h] = signrank(corr_ves_com,corr_vis_com); % paired data
[p,h] = ttest(corr_ves_com,corr_vis_com); % paired data




%% 




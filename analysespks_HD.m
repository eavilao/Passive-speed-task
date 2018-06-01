function results = analysespks_HD(data,modality,prs)

switch modality
    case 'ves'
        speeds = data.stim.speed(data.stim.modality == 1);        
        spks = data.spks(data.stim.modality == 1);
        choice = data.stim.choice(data.stim.modality == 1);
    case 'vis'
        speeds = data.stim.speed(data.stim.modality == 2);
        spks = data.spks(data.stim.modality == 2);
        choice = data.stim.choice(data.stim.modality == 2);
    case 'com'
        speeds = data.stim.speed(data.stim.modality == 3);
        spks = data.spks(data.stim.modality == 3);
        choice = data.stim.choice(data.stim.modality == 3);
end

%% trial-averaged responses 
results.time = data.stim.time;
results.tstim = mean(data.stim.tstim) - (prs.nspk(1)-prs.nspk(2));
results.choice = choice;
results.headings = speeds; 

all_speeds = unique(speeds);
correctTrials = choice==0; incorrectTrials =  choice==5;
for i=1:length(all_speeds)
    results.stim(i) = all_speeds(i);
    spk = spks(speeds == all_speeds(i)); ntrls = length(spk);
    spk_corr = spks(speeds == all_speeds(i) & correctTrials);
    spk_err = spks(speeds == all_speeds(i) & incorrectTrials);
    
    results.nspk{i,:}=[spk.nspk];
    results.rate_pst(i,:) = hist(cell2mat({spk.tspk}'),data.stim.time)/(ntrls*prs.dt);
    results.rate_pst(i,:) = smooth_pst(results.rate_pst(i,:),prs.dt,prs.tsmooth);
    results.rate_avg(i).mu = mean([spk.nspk]/results.tstim); 
    results.rate_avg(i).sig = std([spk.nspk]/results.tstim)/sqrt(ntrls);
    
    % correct choices
    results.rate_pst_correct(i,:) = hist(cell2mat({spk_corr.tspk}'),data.stim.time)/(ntrls*prs.dt);
    results.rate_pst_correct(i,:) = smooth_pst(results.rate_pst_correct(i,:),prs.dt,prs.tsmooth);
    results.rate_avg_correct(i).mu = mean([spk_corr.nspk]/results.tstim);  
    results.rate_avg_correct(i).sig = std([spk_corr.nspk]/results.tstim)/sqrt(ntrls);
    
    % incorrect choices
    results.rate_pst_err(i,:) = hist(cell2mat({spk_err.tspk}'),data.stim.time)/(ntrls*prs.dt);
    results.rate_pst_err(i,:) = smooth_pst(results.rate_pst_err(i,:),prs.dt,prs.tsmooth);
    results.rate_avg_err(i).mu = mean([spk_err.nspk]/results.tstim);
    results.rate_avg_err(i).sig = std([spk_err.nspk]/results.tstim)/sqrt(ntrls);
    
    % compute CP
    if length(spk_corr)>=3 & length(spk_err)>=3
        results.stats.CP(i) =  newROC([spk_corr.nspk],[spk_err.nspk],1); %compute choice prob
    else
        results.stats.CP(i) = NaN;
    end
end
%% statistical tests
if ~strcmp(modality,'null')
    alpha = prs.alpha; %1-(1-prs.alpha)^(1/length(all_speeds));
    for i=1:length(all_speeds) % not speeds but headings
        spk = spks(speeds == all_speeds(i)); ntrls = length(spk(1:16));
        spk = spk(1:16); 
        nspk(i,:) = [spk.nspk]; 
        tspk = {spk.tspk}; t_prestim = prs.tbeg_acc - prs.tstim_on;
        for j=1:length(spk)
            nspk_prestim(i,j) = length(tspk{j}(tspk{j}>prs.tstim_on & tspk{j}<prs.tbeg_acc));
        end
%         [~,pval_exc(i)] = ttest(nspk(i,:),nspk_0,'tail','right');
%         [~,pval_sup(i)] = ttest(nspk(i,:),nspk_0,'tail','left');
        r(i) = mean(nspk(i,:))/results.tstim;
        results.rate_prestim(i).mu = mean(nspk_prestim(i,:))/t_prestim;
        results.rate_prestim(i).sig = (std(nspk_prestim(i,:))/t_prestim)/sqrt(ntrls);
        sse(i) = sum(((nspk(i,:) - mean((nspk(i,:))))/results.tstim).^2);
    end
    % discrimination index
    rmax = max(r); rmin = min(r); sse = sum(sse);
    N = (length(all_speeds) + 1)*length(spk); M = length(all_speeds) + 1;
    results.stats.discindx = (rmax - rmin)/((rmax - rmin) + 2*sqrt(sse/(N - M)));
    % excitatory-ness in the middle of the trial
%     results.stats.pvals.exc_mid = min(pval_exc); results.stats.flags.exc_mid = (results.stats.pvals.exc_mid <= alpha);
    % suppressive-ness in the middle of the trial
%     results.stats.pvals.sup_mid = min(pval_sup); results.stats.flags.sup_mid = (results.stats.pvals.sup_mid <= alpha);
    % preferred stimulus
    [~,ind_pref]=max([results.rate_avg.mu]);
    results.stats.stim.pref = results.stim(ind_pref);
    % anti-preferred stimulus
    [~,ind_antipref]=min([results.rate_avg.mu]);
    results.stats.stim.antipref = results.stim(ind_antipref);
    % event times
    results.stats.events = compute_eventtimes(spks,prs);    % response latency, rise time, peak time
    % excitatory-ness
    results.stats.flags.exc = strcmp(results.stats.events.rise.type_on,'exc') | ...
        strcmp(results.stats.events.rise.type_off,'exc');
    % suppressive-ness
    results.stats.flags.sup = strcmp(results.stats.events.rise.type_on,'sup') | ...
        strcmp(results.stats.events.rise.type_off,'sup');
end

% compute choice correlation
results.npsk_choice_heading = [([spks.nspk])' choice' speeds'];
[rho,pval] = partialcorr(results.npsk_choice_heading); % choice controlling for heading

results.parCorr.rho_spks_choice = rho(2);
results.parCorr.pval_spks_choice = pval(2);

results.parCorr.rho_spks_heading = rho(3);
results.parCorr.pval_spks_heading = pval(3);

% correlation
[rho2,pval2]=corrcoef(results.npsk_choice_heading);
results.corr.rho_spks_choice = rho2(2); 
results.corr.pval_spks_choice = pval2(2); 

results.corr.rho_spks_heading = rho2(3); 
results.corr.pval_spks_heading = pval2(3); 

%% com-ves and com-vis correlations
if strcmp(modality,'com')
    r_ves = [data.ves.rate_avg.mu];
    r_vis = [data.vis.rate_avg.mu];
    r_com = [results.rate_avg.mu];
    [results.stats.rves_avg,results.stats.pves_avg] = corr(r_ves',r_com');
    [results.stats.rvis_avg,results.stats.pvis_avg] = corr(r_vis',r_com');
    
    t_indx = (results.time>prs.tstim_on & results.time<prs.tstim_off);
    r_ves = mean(data.ves.rate_pst(:,t_indx));
    r_vis = mean(data.vis.rate_pst(:,t_indx));
    r_com = mean(results.rate_pst(:,t_indx));
    [results.stats.rves_pst,results.stats.pves_pst] = corr(r_ves',r_com');
    [results.stats.rvis_pst,results.stats.pvis_pst] = corr(r_vis',r_com');
end



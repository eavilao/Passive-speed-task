function results = analysespks_1DAzi_Stefania(data,modality,prs)

switch modality
    case 'null'
        speeds = data.stim.speed(data.stim.modality == -1);
        spks = data.spks(data.stim.modality == -1);
    case 'ves_nofix'
        speeds = data.stim.speed(data.stim.modality == 0);
        spks = data.spks(data.stim.modality == 0);
    case 'ves'
        data.stim.speed = data.stim.speed;
        data.stim.modality = data.stim.modality;
        speeds = data.stim.speed(data.stim.modality == 1);        
        spks = data.spks(data.stim.modality == 1);
    case 'vis'
        data.stim.speed = data.stim.speed;
        data.stim.modality = data.stim.modality;
        speeds = data.stim.speed(data.stim.modality == 2);
        spks = data.spks(data.stim.modality == 2);
    case 'com'
        data.stim.speed = data.stim.speed;
        data.stim.modality = data.stim.modality;
        speeds = data.stim.speed(data.stim.modality == 3);
        spks = data.spks(data.stim.modality == 3);
end

%% trial-averaged responses
results.time = data.stim.time;
results.tstim = mean(data.stim.tstim) - (prs.nspk(1)-prs.nspk(2));
all_speeds = unique(speeds);
for i=1:length(all_speeds)
    results.stim(i) = all_speeds(i);
    spk = spks(speeds == all_speeds(i)); ntrls = length(spk);
    results.rate_pst(i,:) = hist(cell2mat({spk.tspk}'),data.stim.time)/(ntrls*prs.dt);
    results.rate_pst(i,:) = smooth_pst(results.rate_pst(i,:),prs.dt,prs.tsmooth);
    results.rate_avg(i).mu = mean([spk.nspk]/results.tstim);
    results.rate_avg(i).sig = std([spk.nspk]/results.tstim)/sqrt(ntrls);
end

%% statistical tests
if ~strcmp(modality,'null')
    nspk_0 = [data.spks(data.stim.modality==-1).nspk];
    r_0 = mean(nspk_0)/results.tstim;
    sse_0 = sum(((nspk_0 - mean(nspk_0))/results.tstim).^2);
    alpha = prs.alpha; %1-(1-prs.alpha)^(1/length(all_speeds));
    for i=1:length(all_speeds)
        spk = spks(speeds == all_speeds(i)); ntrls = length(spk);
        nspk(i,:) = [spk.nspk]; 
        tspk = {spk.tspk}; t_prestim = prs.tbeg_acc - prs.tstim_on;
        for j=1:length(spk)
            nspk_prestim(i,j) = length(tspk{j}(tspk{j}>prs.tstim_on & tspk{j}<prs.tbeg_acc));
        end
        [~,pval_exc(i)] = ttest(nspk(i,:),nspk_0,'tail','right');
        [~,pval_sup(i)] = ttest(nspk(i,:),nspk_0,'tail','left');
        r(i) = mean(nspk(i,:))/results.tstim;
        results.rate_prestim(i).mu = mean(nspk_prestim(i,:))/t_prestim;
        results.rate_prestim(i).sig = (std(nspk_prestim(i,:))/t_prestim)/sqrt(ntrls);
        sse(i) = sum(((nspk(i,:) - mean((nspk(i,:))))/results.tstim).^2);
    end
    % discrimination index
    rmax = max([r_0 r]); rmin = min([r_0 r]); sse = sse_0 + sum(sse);
    N = (length(all_speeds) + 1)*length(spk); M = length(all_speeds) + 1;
    results.stats.discindx = (rmax - rmin)/((rmax - rmin) + 2*sqrt(sse/(N - M)));
    % excitatory-ness in the middle of the trial
    results.stats.pvals.exc_mid = min(pval_exc); results.stats.flags.exc_mid = (results.stats.pvals.exc_mid <= alpha);
    % suppressive-ness in the middle of the trial
    results.stats.pvals.sup_mid = min(pval_sup); results.stats.flags.sup_mid = (results.stats.pvals.sup_mid <= alpha);
    % stimulus tuning
    results.stats.pvals.tuning = anova1([nspk_0 ; nspk]',[],'off'); results.stats.flags.tuning = (results.stats.pvals.tuning <= prs.alpha);
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
    % preferred direction for angular rotation
    if unique(sign(results.stim))==2
        [hr,pr] = ttest([spks(speeds<0).nspk],[spks(speeds>0).nspk],'tail','left');
        [hl,pl] = ttest([spks(speeds<0).nspk],[spks(speeds>0).nspk],'tail','right');
        results.stats.direction.pref = hr-hl;
        results.stats.direction.pval = min(pr,pl);
    end
end

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
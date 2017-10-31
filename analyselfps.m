function results = analyselfps(data,modality,prs)

switch modality
    case 'null'
        speeds = data.stim.speed(data.stim.modality == -1);
        waves = data.wave(data.stim.modality == -1);
    case 'ves_nofix'
        speeds = data.stim.speed(data.stim.modality == 0);
        waves = data.wave(data.stim.modality == 0);
    case 'ves'
        speeds = data.stim.speed(data.stim.modality == 1);
        waves = data.wave(data.stim.modality == 1);
    case 'vis'
        speeds = data.stim.speed(data.stim.modality == 2);
        waves = data.wave(data.stim.modality == 2);
    case 'com'
        speeds = data.stim.speed(data.stim.modality == 3);
        waves = data.wave(data.stim.modality == 3);
end

%% trial-averaged responses
results.time = waves(1).t;
all_speeds = unique(speeds);
for i=1:length(all_speeds)
    results.stim(i) = all_speeds(i);
    wave = waves(speeds == all_speeds(i)); ntrls = length(wave); Pww = [];
    for j=1:length(wave),[Pww(:,j),results.freq] = pmtm(wave(j).v(results.time>prs.tbeg_acc & ...
            results.time<prs.tbeg_dec+0.1),2,[],200); end
    results.wave_psd(i,:) = mean(Pww,2);
    results.wave_pst(i,:) = mean(cell2mat({wave.v}),2);
    results.wave_pst(i,:) = smooth_pst(results.wave_pst(i,:),results.time(2)-results.time(1),prs.tsmooth);
end

%% statistics
if ~strcmp(modality,'null')
    wave_pst = mean(results.wave_pst);
    SEP_on = min(wave_pst(results.time > prs.tstim_on & results.time < (prs.tstim_on + 0.4)));
    SEP_acc = min(wave_pst(results.time > prs.tbeg_acc & results.time < (prs.tbeg_acc + 0.5)));
    SEP_dec = min(wave_pst(results.time > prs.tbeg_dec & results.time < (prs.tbeg_dec + 0.5)));
    SEP_off = min(wave_pst(results.time > prs.tstim_off & results.time < (prs.tstim_off + 0.4)));
    lfp_dc = mean(wave_pst(results.time>prs.tstim_on & results.time<prs.tbeg_acc));
    % response latency to stimulus onset
    results.stats.SEP(1).t_peak = results.time(wave_pst==SEP_on); 
    results.stats.SEP(1).type_peak = {'stim_on'};
    results.stats.SEP(1).amplitude = wave_pst(wave_pst==SEP_on) - lfp_dc;
    % response latency to accelaration
    results.stats.SEP(2).t_peak = results.time(wave_pst==SEP_acc); 
    results.stats.SEP(2).type_peak = {'acc'};
    results.stats.SEP(2).amplitude = wave_pst(wave_pst==SEP_acc) - lfp_dc;
    % response latency to deccelaration
    results.stats.SEP(3).t_peak = results.time(wave_pst==SEP_dec); 
    results.stats.SEP(3).type_peak = {'dec'};
    results.stats.SEP(3).amplitude = wave_pst(wave_pst==SEP_dec) - lfp_dc;
    % response latency to stimulus offset
    results.stats.SEP(4).t_peak = results.time(wave_pst==SEP_off); 
    results.stats.SEP(4).type_peak = {'stim_off'};
    results.stats.SEP(4).amplitude = wave_pst(wave_pst==SEP_off) - lfp_dc;
    % preferred stimulus
    [~,ind_pref] = min(results.wave_pst(:,wave_pst==SEP_acc)); % most negative is preferred
    results.stats.stim.pref = results.stim(ind_pref);
    % anti-preferred stimulus
    [~,ind_antipref] = max(results.wave_pst(:,wave_pst==SEP_acc)); % least negative is antipreferred
    results.stats.stim.antipref = results.stim(ind_antipref);
    % event times
    results.stats.events = compute_eventtimes(waves,prs);    % response latency, rise time, peak time
    % excitatory-ness
    results.stats.flags.exc = strcmp(results.stats.events.rise.type_on,'exc');
    % suppressive-ness
    results.stats.flags.sup = strcmp(results.stats.events.rise.type_on,'sup');
end

%% com-ves and com-vis correlations
if strcmp(modality,'com')
    r_ves = mean(data.ves.wave_pst);
    r_vis = mean(data.vis.wave_pst);
    r_com = mean(results.wave_pst);
    [results.stats.rves_pst,results.stats.pves_pst] = corr(r_ves',r_com');
    [results.stats.rvis_pst,results.stats.pvis_pst] = corr(r_vis',r_com');
end
function results = analysecorr(stim,nspk,modality)

speeds = stim.speed;
modalities = stim.modality;

switch modality
    case 'ves_nofix'
        indx = modalities==0;
    case 'ves'
        indx = modalities==1;
    case 'vis'
        indx = modalities==2;
    case 'com'
        indx = modalities==3;
end

speeds = speeds(indx);
nspk = nspk(:,indx);

%% correlations
all_speeds = unique(speeds);
results.stim = all_speeds;
for i=1:length(all_speeds)
    indx = (speeds == all_speeds(i));
    mu_nspk(:,i) = mean(nspk(:,indx),2);
    mu_nspk2 = repmat(mean(nspk(:,indx),2),[1 sum(indx)]);
    sig_nspk = repmat(std(nspk(:,indx),1,2),[1 sum(indx)]);
    z_nspk(:,indx) = (nspk(:,indx) - mu_nspk2)./sig_nspk;                   % z-scored spike counts
    results.noise.r_vs_stim(i,:,:) = corr(nspk(:,indx)');
end
[results.signal.r,results.signal.p]=corr(mu_nspk');
[results.noise.r,results.noise.p]=corr(z_nspk');

%% slopes
drds1 = [];
if any(all_speeds<0)
    % negative velocities
    dr_temp = diff(mu_nspk(:,all_speeds<0)')';
    ndr = size(dr_temp,2); dr = [];
    for i=2:ndr
        dr(:,i) = (dr_temp(:,i-1)+dr_temp(:,i))/2;
    end
    dr(:,1) = dr_temp(:,1); dr(:,end+1) = dr_temp(:,end);
    drds1 = dr/unique(diff(all_speeds(all_speeds<0)));
end

% positive velocities
dr_temp = diff(mu_nspk(:,all_speeds>0)')';
ndr = size(dr_temp,2); dr = [];
for i=2:ndr
    dr(:,i) = (dr_temp(:,i-1)+dr_temp(:,i))/2;
end
dr(:,1) = dr_temp(:,1); dr(:,end+1) = dr_temp(:,end);
drds2 = dr/unique(diff(all_speeds(all_speeds>0)));
results.signal.dfds = [drds1 drds2];


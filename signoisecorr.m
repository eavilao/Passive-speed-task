function [r_sig,p_sig,r_noise,p_noise] = signoisecorr(stim,nspk1,nspk2,modality)

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
nspk10 = nspk1(modalities==-1); nspk1 = nspk1(indx);
nspk20 = nspk2(modalities==-1); nspk2 = nspk2(indx);


%% correlations

all_speeds = unique(speeds);
for i=1:length(all_speeds)
    indx = (speeds == all_speeds(i));
    mu_nspk1(i) = mean(nspk1(indx)); sig_nspk1(i) = std(nspk1(indx));
    mu_nspk2(i) = mean(nspk2(indx)); sig_nspk2(i) = std(nspk2(indx));
    z_nspk1(indx) = (nspk1(indx) - mu_nspk1(i))./sig_nspk1(i);              % z-scored spike counts
    z_nspk2(indx) = (nspk2(indx) - mu_nspk2(i))./sig_nspk2(i);              % z-scored spike counts
end
[r_sig,p_sig] = corr([mean(nspk10) mu_nspk1]',[mean(nspk20) mu_nspk2]');
[r_noise,p_noise] = corr(z_nspk1',z_nspk2');

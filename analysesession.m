function corr = analysesession(data)

%% average data from sessions (second-order statistics)
thisdata = data{2};
nsessions = length(thisdata);

% distance correlations
if isfield(thisdata(1).corr_distance,'r_nspk')
    % single/multiunits
    dx = [];
    r_nspk = []; r_tspk = [];
    rves_sig = []; rves_noise = [];
    rvis_sig = []; rvis_noise = [];
    for i=1:nsessions
        dx = [dx ; thisdata(i).corr_distance.dx(:)];
        r_nspk = [r_nspk ; thisdata(i).corr_distance.r_nspk(:)];
        r_tspk = [r_tspk ; thisdata(i).corr_distance.r_tspk(:)];
        rves_sig = [rves_sig ; thisdata(i).corr_distance.ves.r_sig(:)];
        rves_noise = [rves_noise ; thisdata(i).corr_distance.ves.r_noise(:)];
        rvis_sig = [rvis_sig ; thisdata(i).corr_distance.vis.r_sig(:)];
        rvis_noise = [rvis_noise ; thisdata(i).corr_distance.vis.r_noise(:)];
    end
    groupsize = 1;
    [~,r_nspk]=groupmean(dx,r_nspk,groupsize);
    [~,r_tspk]=groupmean(dx,r_tspk,groupsize);
    [~,rves_sig]=groupmean(dx,rves_sig,groupsize);
    [~,rves_noise]=groupmean(dx,rves_noise,groupsize);
    [~,rvis_sig]=groupmean(dx,rvis_sig,groupsize);
    [dx2,rvis_noise]=groupmean(dx,rvis_noise,groupsize);
    corr.dx = dx2;
    corr.r_nspk = r_nspk; corr.r_tspk = r_tspk;
    corr.ves.r_sig = rves_sig; corr.ves.r_noise = rves_noise;
    corr.vis.r_sig = rvis_sig; corr.vis.r_noise = rvis_noise;
else
    % lfps
    dx = [];
    r = []; rves = []; rvis = [];
    for i=1:nsessions
        dx = [dx ; thisdata(i).corr_distance.dx(:)];
        r = [r ; thisdata(i).corr_distance.r_mu(:)];
        rves = [rves ; thisdata(i).corr_distance.ves.r_mu(:)];
        rvis = [rvis ; thisdata(i).corr_distance.vis.r_mu(:)];
    end
    [~,r]=groupmean(dx,r,1);
    [~,rves]=groupmean(dx,rves,1);
    [dx2,rvis]=groupmean(dx,rvis,1);
    corr.dx = dx2;
    corr.r = r;
    corr.ves.r = rves; corr.vis.r = rvis;
end
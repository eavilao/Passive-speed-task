function corr2 = corr_ch2dist(corr,interchdist,exp_name)

x = corr.x;

if strcmp(exp_name,'HD')
    
    if isfield(corr,'r_nspk')
        r_nspk = corr.r_nspk;
        r_tspk = corr.r_tspk;
        rves_sig = corr.ves.r_sig;
        rvis_sig = corr.vis.r_sig;
        rves_noise = corr.ves.r_noise;
        rvis_noise = corr.vis.r_noise;
    else
        r = corr.r;
        r_ves = corr.ves.r;
        r_vis = corr.vis.r;
    end
    
else
    if isfield(corr,'r_nspk')
        r_nspk = corr.r_nspk;
        r_tspk = corr.r_tspk;
        rves_sig = corr.ves.r_sig;
        rvis_sig = corr.vis.r_sig;
        rves_noise = corr.ves.r_noise;
        rvis_noise = corr.vis.r_noise;
    else
        r = corr.r;
        r_null = corr.null.r;
        r_ves = corr.ves.r;
        r_vis = corr.vis.r;
    end
    
end




if ~isfield(corr,'t')
    if isfield(corr,'r_nspk')
        if length(unique(x))==max(x)% multiunit
            corr2.dx = (x-x(1))*interchdist;
            for i=1:length(x)
                corr2.r_nspk(i) = mean(diag(r_nspk,i-1));
                corr2.r_tspk(i) = mean(diag(r_tspk,i-1));
                corr2.ves.r_sig(i) = mean(diag(rves_sig,i-1));
                corr2.vis.r_sig(i) = mean(diag(rvis_sig,i-1));
                corr2.ves.r_noise(i) = mean(diag(rves_noise,i-1));
                corr2.vis.r_noise(i) = mean(diag(rvis_noise,i-1));
            end
        else% singleunit
            corr2.dx = [];
            corr2.r_nspk = []; corr2.r_tspk = [];
            corr2.ves.r_sig = []; corr2.vis.r_sig = [];
            corr2.ves.r_noise = []; corr2.vis.r_noise = [];
            for i=1:length(x)
                for j=i+1:length(x)
                    corr2.dx = [corr2.dx ; (x(j)-x(i))*interchdist];
                    corr2.r_nspk = [corr2.r_nspk ; r_nspk(i,j)];
                    corr2.r_tspk = [corr2.r_tspk ; r_tspk(i,j)];
                    corr2.ves.r_sig = [corr2.ves.r_sig ; rves_sig(i,j)];
                    corr2.vis.r_sig = [corr2.vis.r_sig ; rvis_sig(i,j)];
                    corr2.ves.r_noise = [corr2.ves.r_noise ; rves_noise(i,j)];
                    corr2.vis.r_noise = [corr2.vis.r_noise ; rvis_noise(i,j)];
                end
            end
        end
    else % lfp
        corr2.dx = (x-x(1))*interchdist;
        if strcmp(exp_name,'HD')
            for i=1:length(x)
                corr2.r_mu(i) = mean(diag(r,i-1));
                corr2.ves.r_mu(i) = mean(diag(r_ves,i-1));
                corr2.vis.r_mu(i) = mean(diag(r_vis,i-1));
            end
        else
            for i=1:length(x)
                corr2.r_mu(i) = mean(diag(r,i-1));
                corr2.null.r_mu(i) = mean(diag(r_null,i-1));
                corr2.ves.r_mu(i) = mean(diag(r_ves,i-1));
                corr2.vis.r_mu(i) = mean(diag(r_vis,i-1));
            end
        end
    end
else
    corr.dx = []; corr.r_mu = [];
    for i=1:length(x)
        for j=i+1:length(x)
            corr.dx = [corr.dx ; (x(j)-x(i))*interchdist];
            corr.r_mu = [corr.r_mu ; squeeze(r(i,j,:))'];
        end
    end
end
    
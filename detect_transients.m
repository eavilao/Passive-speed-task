function latencies = detect_transients(spks,prs)

tstim_on = prs.tstim_on; tstim_off = prs.tstim_off;
twin = prs.twin; dt = prs.dt;
t = linspace(tstim_on,tstim_off-twin,(tstim_off-twin-tstim_on)/dt);
t = [t;t+twin];
for i=1:length(spks)
    tspk = spks(i).tspk;
    for j=1:length(t)
        nspk(i,j) = length(tspk(tspk>t(1,j) & tspk<t(2,j)));
    end
end
nspk_mu = mean(nspk); t = mean(t);

tbeg_acc = prs.tbeg_acc; exc=0; sup=0; latencies.on_time = NaN;
%on transient
indx_beg = find(t>tbeg_acc,1); indx_end = find(t>tbeg_acc+.3,1);
for indx = indx_beg:indx_end
    r_mu = mean(nspk_mu(indx-3:indx+3));
    r_thresh1 = mean(nspk_mu(indx-10:indx-4)) + 2*std(nspk_mu(indx-10:indx-4));
    r_thresh2 = mean(nspk_mu(indx-10:indx-4)) - 2*std(nspk_mu(indx-10:indx-4));
    if (r_thresh1>1 && r_mu>r_thresh1), exc=exc+1; else exc=0; end
    if (r_thresh2>1 && r_mu<r_thresh2), sup=sup+1; else sup=0; end
    if exc==5, latencies.on_time = t(indx-4); latencies.on_type = 'exc';
        break;
    elseif sup==5, latencies.on_time = t(indx-4); latencies.on_type = 'sup';
        break;
    end
end

tbeg_dec = prs.tbeg_dec; exc=0; sup=0; latencies.off_time = NaN;
%off transient
indx_beg = find(t>tbeg_dec,1); indx_end = find(t>tbeg_dec+.3,1);
for indx = indx_beg:indx_end
    r_mu = mean(nspk_mu(indx-3:indx+3));
    r_thresh1 = mean(nspk_mu(indx-10:indx-4)) + 2*std(nspk_mu(indx-10:indx-4));
    r_thresh2 = mean(nspk_mu(indx-10:indx-4)) - 2*std(nspk_mu(indx-10:indx-4));
    if r_mu>r_thresh1, exc=exc+1; else exc=0; end
    if r_mu<r_thresh2, sup=sup+1; else sup=0; end
    if exc==5, latencies.off_time = t(indx-4); latencies.off_type = 'exc';
        break;
    elseif sup==5, latencies.off_time = t(indx-4); latencies.off_type = 'sup';
        break;
    end
end
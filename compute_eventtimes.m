function events = compute_eventtimes(spks,prs)

tstim_on = prs.tstim_on; tstim_off = prs.tstim_off;
twin = prs.twin; dt = prs.dt;
twin_current = prs.twin_current; twin_previous = prs.twin_previous;
t = linspace(tstim_on,tstim_off-twin,(tstim_off-twin-tstim_on)/dt);

if ~isfield(spks,'t') % spiketrain
    t = [t;t+twin];
    for i=1:length(spks)
        tspk = spks(i).tspk;
        for j=1:length(t)
            nspk(i,j) = length(tspk(tspk>t(1,j) & tspk<t(2,j)));
        end
    end
    nspk_mu = mean(nspk); t = mean(t);
else % lfp
    for i=1:length(spks)
        v = spks(i).v; t2 = spks(i).t;
        v2(i,:) = interp1(t2,v,t);
    end
    nspk_mu = mean(v2);
end

tbeg_acc = prs.tbeg_acc; exc=0; sup=0; 
events.rise.t_on = NaN; events.rise.type_on = 'nan';
%on transient
indx_beg = find(t>tbeg_acc,1); indx_end = find(t>tbeg_acc+.25,1);
winsize_current = twin_current/dt;
winsize_previous = twin_previous/dt;
for indx = indx_beg:indx_end
    r_mu = mean(nspk_mu(round(indx-winsize_current/2):round(indx+winsize_current/2)));
    r_thresh1 = mean(nspk_mu(round(indx-winsize_previous):indx)) + ...
        2*std(nspk_mu(round(indx-winsize_previous):indx));
    r_thresh2 = mean(nspk_mu(round(indx-winsize_previous):indx)) - ...
        2*std(nspk_mu(round(indx-winsize_previous):indx));
    if ~isfield(spks,'t') % spiketrain
        if (r_thresh1>1 && r_mu>r_thresh1), exc=exc+1; else exc=0; end
        if (r_thresh2>1 && r_mu<r_thresh2), sup=sup+1; else sup=0; end
    else % lfp
        if (r_mu>r_thresh1), exc=exc+1; else exc=0; end
        if (r_mu<r_thresh2), sup=sup+1; else sup=0; end
    end
    if exc==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
        break;
    elseif sup==3, events.rise.t_on = t(indx); events.rise.type_on = 'sup';
        break;
    end
end

tbeg_dec = prs.tbeg_dec; exc=0; sup=0; 
events.rise.t_off = NaN; events.rise.type_off = 'nan';
%off transient
indx_beg = find(t>tbeg_dec,1); indx_end = find(t>tbeg_dec+.25,1);
for indx = indx_beg:indx_end
    r_mu = mean(nspk_mu(round(indx-winsize_current/2):round(indx+winsize_current/2)));
    r_thresh1 = mean(nspk_mu(round(indx-winsize_previous):round(indx-winsize_current/2))) + ...
        2*std(nspk_mu(round(indx-winsize_previous):round(indx-winsize_current/2)));
    r_thresh2 = mean(nspk_mu(round(indx-winsize_previous):round(indx-winsize_current/2))) - ...
        2*std(nspk_mu(round(indx-winsize_previous):round(indx-winsize_current/2)));
    if ~isfield(spks,'t') % spiketrain
        if (r_thresh1>1 && r_mu>r_thresh1), exc=exc+1; else exc=0; end
        if (r_thresh2>1 && r_mu<r_thresh2), sup=sup+1; else sup=0; end
    else % lfp
        if (r_mu>r_thresh1), exc=exc+1; else exc=0; end
        if (r_mu<r_thresh2), sup=sup+1; else sup=0; end
    end
    if exc==3, events.rise.t_off = t(indx); events.rise.type_off = 'exc';
        break;
    elseif sup==3, events.rise.t_off = t(indx); events.rise.type_off = 'sup';
        break;
    end
end

events.latency.t_on = NaN; events.latency.type_on{1} = 'nan';
% onset latency
nspk0.mu = nanmean(nspk_mu(t>0 & t<=tbeg_acc));
nspk0.sig = nanstd(nspk_mu(t>0 & t<=tbeg_acc));
exc_indx = (nspk_mu(t>tbeg_acc & t<=prs.tbeg_dec) > nspk0.mu + 2*nspk0.sig);
for i=1:sum(exc_indx)
    indx = find(exc_indx,i); indx = indx(end);
    if indx<=(length(exc_indx)-10) && all(exc_indx(indx:indx+0.1/dt))
        t2 = t(t>tbeg_acc & t<=prs.tbeg_dec);
        events.latency.t_on = t2(indx);
        events.latency.type_on{1} = 'exc';
        break;
    end
end
sup_indx = (nspk_mu(t>tbeg_acc & t<=prs.tbeg_dec) < nspk0.mu - 2*nspk0.sig);  
for i=1:sum(sup_indx)
    indx = find(sup_indx,i); indx = indx(end);
    if indx<=(length(sup_indx)-10) && all(sup_indx(indx:indx+0.1/dt))
        t2 = t(t>tbeg_acc & t<=prs.tbeg_dec);
        if isnan(events.latency.t_on)
            events.latency.t_on = t2(indx);
            events.latency.type_on{1} = 'sup';
        else
            events.latency.t_on(2) = t2(indx);
            events.latency.type_on{2} = 'sup';
        end
        break;
    end
end

events.peak.t_on = NaN; events.peak.type_on = 'nan';
% time-to-peak
if ~isnan(events.latency.t_on)
    nspk2 = nspk_mu(t>events.latency.t_on(1) & t<events.latency.t_on(1)+0.2);
    t2 = t(t>events.latency.t_on(1) & t<events.latency.t_on(1)+0.2);
    if strcmp(events.latency.type_on{1},'exc')
        events.peak.t_on = t2(find(nspk2==max(nspk2),1));
        events.peak.type_on = 'exc';
    else
        events.peak.t_on = t2(find(nspk2==min(nspk2),1));
        events.peak.type_on = 'sup';
    end
end
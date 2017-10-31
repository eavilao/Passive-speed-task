function plottrials(lfps,multiunits,singleunits,prs)

t_beg = prs.tstim_on - 0.1;
t_end = prs.tstim_off;
nch = length(lfps);
nunits = length(singleunits); clrs = brewermap(nunits,'paired');
ntrls = length(lfps(1).wave);

for i=1:ntrls
    figure(1); set(gcf,'Position',[100 100 400 400]);
    figure(2); set(gcf,'Position',[600 100 400 400]);
    for k=1:nch
       figure(1); hold on;
       plot(lfps(k).wave(i).t,lfps(k).wave(i).v - 0.1*k,'k');
       figure(2); hold on;
       tspk = multiunits(k).spks(i).tspk;
       tspk = tspk(tspk>t_beg & tspk<t_end);
       line([tspk' ; tspk'], [zeros(1,length(tspk)) ; ones(1,length(tspk))] - k,'Color',[.5 .5 .5]);
    end
    for k=1:nunits
        figure(2); hold on;
        ch_no = singleunits(k).channel_no;
        tspk = singleunits(k).spks(i).tspk;
        tspk = tspk(tspk>t_beg & tspk<t_end);
        line([tspk' ; tspk'], [zeros(1,length(tspk)) ; ones(1,length(tspk))] - ch_no,...
            'Color',clrs(k,:),'Linewidth',2);
    end
    figure(1); xlim([t_beg t_end]);
    figure(2); xlim([t_beg t_end]);
    waitforbuttonpress; close all;
end
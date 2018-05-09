% Hartigan dip test for unimodality

% linearspeed
units_width = [experiments(1).singleunits(:).spkwf]; 
spk_width = [units_width(:).width]; 
centers = [50:25:600]; 

figure; hold on
histogram(spk_width,centers, 'Facecolor', 'none'); 
set(gca, 'TickDir', 'out', 'Fontsize', 18, 'xTick', [0 200 400 600])
box off; 
vline(prctile(spk_width',33)); % 300 narrow
vline(prctile(spk_width',66)); % 475 broad

%% plot psth
% gather data
cond = {'ves','vis','com'};
pop = experiments(1).populations(1);
prs = default_prs('linearspeed',0);


% indx narrow and broad
narrow_ves = pop.ves.exc.indx & spk_width < 300; broad_ves = pop.ves.exc.indx & spk_width > 475; % pctile 33 and 66
narrow_vis = pop.vis.exc.indx & spk_width < 300; broad_vis = pop.vis.exc.indx & spk_width > 475;
narrow_com = pop.com.exc.indx & spk_width < 300; broad_com = pop.com.exc.indx & spk_width > 475;

units = experiments(1).singleunits(broad_com);
for i=1:length(units)
    [nstim(i),ntime(i)] = size(units(i).com.rate_pst);
end
nunits = 1:length(units);
for j=1:min(nstim)
            for i=nunits
                rate_pst(i,j,:) = smooth_pst(units(i).ves.rate_pst(j,1:min(ntime)),prs.dt,prs.tsmooth_pop);
                rate_pst0(i,j,:) = smooth_pst(units(i).null.rate_pst(1,1:min(ntime)),prs.dt,prs.tsmooth_pop);
            end   
            broad_ves_rate_mu(j,:) = squeeze(mean(rate_pst(nunits,j,:)));
            broad_ves_rate_sig(j,:) = squeeze(std(rate_pst(nunits,j,:)))/sqrt(length(units));
            broad_ves_rate0_mu(j,:) = squeeze(mean(rate_pst0(nunits,j,:)));
            broad_ves_rate0_sig(j,:) = squeeze(std(rate_pst0(nunits,j,:)))/sqrt(length(units));
end

t = pop.com.exc.rate_pst.time;

% plot psth
    figure; hold on;
    colorscale = 1/size(broad_vis_rate_mu,1):1/size(broad_vis_rate_mu,1):1;
    plot(t,broad_ves_rate0_mu(1,:),'Color','k','linewidth',2);
    for j=1:size(broad_ves_rate_mu,1)
        plot(t,broad_ves_rate_mu(j,:),'Color',colorscale(j)*([0 0 1]),'linewidth',2);
    end
    set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
        'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
    ylim([5 30])

 % plot tuning curve 
 % gather
 s = units(1).ves.stim;
 for i=nunits
      r(i,:) = [units(i).com.rate_avg.mu];
      e(i,:) = [units(i).com.rate_avg.sig];
      r0(i) = [units(i).null.rate_avg.mu];
 end

 % plot
        hold on;
        for i=1:size(r,1)
            errorbar(s,mean(r),mean(e)/2,'Color',[0 0 1],'Linewidth',2);
        end
        hline(mean(r0),'--k');
        set(gca,'xlim',[min(s)-min(unique(diff(s))) max(s)+min(unique(diff(s)))],'XTick',s,...
            'XTickLabel',s*1e2,'TickDir','Out','Fontsize',16);

%%
% pick only exc neurons
resp_ves_exc = experiments(1).populations(1).ves.exc.indx;
resp_vis_exc = experiments(1).populations(1).vis.exc.indx;
ves_vis_exc = resp_ves_exc & resp_vis_exc;

ves_width = spk_width(resp_ves_exc); 
vis_width = spk_width(resp_vis_exc);
ves_vis_width = spk_width(ves_vis_exc);


% angspeed
units_ang = [experiments(2).singleunits(:).spkwf]; 
spk_width_ang = [units_ang(:).width]; 

hist(spk_width_ang,20)

% pick only exc neurons

resp_ves_exc = experiments(1).populations(1).ves.exc.indx;
resp_vis_exc = experiments(1).populations(1).vis.exc.indx;
ves_vis_exc = resp_ves_exc & resp_vis_exc;

ves_width = spk_width(resp_ves_exc); 
vis_width = spk_width(resp_vis_exc);
ves_vis_width = spk_width(ves_vis_exc);
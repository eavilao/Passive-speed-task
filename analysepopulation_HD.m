function results = analysepopulation_HD(data,unitname,modality,prs)

thisdata = [data{1}.(modality)];
nunits = length(thisdata);
switch unitname
    case {'singleunits','multiunits'}
        % population average psth
        for i=1:nunits, [nstim(i),ntime(i)] = size(thisdata(i).rate_pst); end
        for j=1:min(nstim)
            for i=1:nunits
                rate_pst(i,j,:) = smooth_pst(thisdata(i).rate_pst(j,1:min(ntime)),prs.dt,prs.tsmooth_pop);
            end
%             for i=1:nunits
%                 rate_pst(i,j,:) = thisdata(i).rate_pst(j,1:min(ntime));
%             end
            results.all.rate_pst.mu(j,:) = squeeze(mean(rate_pst(:,j,:)));
            results.all.rate_pst.sig(j,:) = squeeze(std(rate_pst(:,j,:)))/sqrt(nunits);
        end
        results.all.rate_pst.time = thisdata(1).time(1:min(ntime));
        
        % population average tuning curves and nspk
        for i=1:nunits
            rate_avg(i,:) = [thisdata(i).rate_avg.mu];
            
            % nspk, choice, heading
            results.all.nspk(:,i) = thisdata(i).npsk_choice_heading(1:179,1);
            results.all.choice(:,i) = thisdata(i).npsk_choice_heading(1:179,2);
            results.all.heading(:,i) = thisdata(i).npsk_choice_heading(1:179,3);
            
        end
        results.all.rate_avg.stim = thisdata(1).stim;
        results.all.rate_avg.mu = mean(rate_avg);
        results.all.rate_avg.sig = std(rate_avg)/sqrt(nunits);
        
        % stim-response correlation
        y=rate_avg; x=repmat(thisdata(1).stim,[size(y,1),1]); [r,p] = corr(x(:),y(:));
        results.all.stats.corr_stimrate.r = r; results.all.stats.corr_stimrate.p = p;
        results.all.rate_avg = y; results.all.stim = x; 
        
        % discindx stats
        for i=1:nunits, d(i) = thisdata(i).stats.discindx; end
        results.all.discindx = d;
        

        % preferred/anti-preferred stimulus
        for i=1:nunits, stim_pref(i) = thisdata(i).stats.stim.pref; end
        for i=1:nunits, stim_antipref(i) = thisdata(i).stats.stim.antipref; end
        results.all.stim_pref = stim_pref; results.all.stim_antipref = stim_antipref;

        % rise time
        for i=1:nunits
            t_on(i) = thisdata(i).stats.events.rise.t_on;
            type_on{i} = thisdata(i).stats.events.rise.type_on;
        end
        for i=1:nunits
            t_off(i) = thisdata(i).stats.events.rise.t_off;
            type_off{i} = thisdata(i).stats.events.rise.type_off;
        end
        results.all.rise.t_on = t_on; results.all.rise.t_off = t_off;
        results.all.rise.type_on = type_on; results.all.rise.type_off = type_off;
        results.all.stats.rise.t_on.mu = nanmedian(results.all.rise.t_on);
        results.all.stats.rise.t_on.sig = nanstd(results.all.rise.t_on)/sqrt(nunits);
        results.all.stats.rise.t_off.mu = nanmedian(results.all.rise.t_off);
        results.all.stats.rise.t_off.sig = nanstd(results.all.rise.t_off)/sqrt(nunits);

        
        % latency
        clear t_on;
        for i=1:nunits, t_on(i) = thisdata(i).stats.events.latency.t_on(1); end
        
        results.all.latency.t_on = t_on;
        results.all.stats.latency.t_on.mu = nanmedian(results.all.latency.t_on);
        results.all.stats.latency.t_on.sig = nanstd(results.all.latency.t_on)/sqrt(nunits);
        
        % peak time
        clear t_on;
        for i=1:nunits, t_on(i) = thisdata(i).stats.events.peak.t_on; end
        
        results.all.peak.t_on = t_on;
        results.all.stats.peak.t_on.mu = nanmedian(results.all.peak.t_on);
        results.all.stats.peak.t_on.sig = nanstd(results.all.peak.t_on)/sqrt(nunits);

    case 'lfps'
        % population average psth
        results.time = thisdata(1).time;
        for i=1:nunits, [nstim(i),ntime(i)] = size(thisdata(i).wave_pst); end
        for j=1:min(nstim)
            for i=1:nunits
                wave_pst(i,j,:) = thisdata(i).wave_pst(j,1:min(ntime));
            end
            results.all.wave_pst.mu(j,:) = squeeze(mean(wave_pst(:,j,:)));
            results.all.wave_pst.sig(j,:) = squeeze(std(wave_pst(:,j,:)))/sqrt(nunits);
        end
        
        % peak time
        for i=1:nunits, results.all.SEP.t_peak(i,:) = [thisdata(i).stats.SEP.t_peak]; end
        for i=1:nunits, results.all.SEP.v_peak(i,:) = [thisdata(i).stats.SEP.amplitude]; end
        for i=1:nunits, results.all.SEP.type_peak{i} = [thisdata(i).stats.SEP.type_peak]; end
        
        % preferred/anti-preferred stimulus
        for i=1:nunits, stim_pref(i) = thisdata(i).stats.stim.pref; end
        for i=1:nunits, stim_antipref(i) = thisdata(i).stats.stim.antipref; end
        results.all.stim_pref = stim_pref; results.all.stim_antipref = stim_antipref;
end

% rise time
for i=1:nunits
    t_on(i) = thisdata(i).stats.events.rise.t_on;
    type_on{i} = thisdata(i).stats.events.rise.type_on;
end
for i=1:nunits
    t_off(i) = thisdata(i).stats.events.rise.t_off;
    type_off{i} = thisdata(i).stats.events.rise.type_off;
end
results.all.rise.t_on = t_on; results.all.rise.t_off = t_off;
results.all.rise.type_on = type_on; results.all.rise.type_off = type_off;
results.all.stats.rise.t_on.mu = nanmedian(results.all.rise.t_on);
results.all.stats.rise.t_on.sig = nanstd(results.all.rise.t_on)/sqrt(nunits);
results.all.stats.rise.t_off.mu = nanmedian(results.all.rise.t_off);
results.all.stats.rise.t_off.sig = nanstd(results.all.rise.t_off)/sqrt(nunits);

% latency
clear t_on;
for i=1:nunits, t_on(i) = thisdata(i).stats.events.latency.t_on(1); end
results.all.latency.t_on = t_on;
results.all.stats.latency.t_on.mu = nanmedian(results.all.latency.t_on);
results.all.stats.latency.t_on.sig = nanstd(results.all.latency.t_on)/sqrt(nunits);

% peak time
clear t_on;
for i=1:nunits, t_on(i) = thisdata(i).stats.events.peak.t_on; end
results.all.peak.t_on = t_on;
results.all.stats.peak.t_on.mu = nanmedian(results.all.peak.t_on);
results.all.stats.peak.t_on.sig = nanstd(results.all.peak.t_on)/sqrt(nunits);
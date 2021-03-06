function results = analysepopulation_1DAzi(data,unitname,modality,prs)

thisdata = [data{1}.(modality)]; thisdata0 = [data{1}.('null')];
nunits = length(thisdata);
switch unitname
    case {'singleunits','multiunits'}
        % index neurons based on responsiveness
        for i=1:nunits
            tuned(i)=thisdata(i).stats.flags.tuning;
            exc(i)=thisdata(i).stats.flags.exc;
            sup(i)=thisdata(i).stats.flags.sup;
            resp(i)=exc(i)|sup(i);
            unresp(i)=~resp(i);
            untuned(i)=~tuned(i) & resp(i);
            tuned_exc(i)=tuned(i) & exc(i);
            tuned_sup(i)=tuned(i) & sup(i);
        end
        results.tuned.indx = tuned;
        results.untuned.indx = untuned;
        results.resp.indx = resp;
        results.unresp.indx = unresp;
        results.exc.indx = exc;
        results.sup.indx = sup;
        results.tuned_exc.indx = tuned_exc;
        results.tuned_sup.indx = tuned_sup;
        
        % percentages
        results.tuned.N_neurons = sum(results.tuned.indx);
        results.tuned.f_neurons = sum(results.tuned.indx)/numel(results.tuned.indx);
        results.exc.N_neurons = sum(results.exc.indx);
        results.exc.f_neurons = sum(results.exc.indx)/numel(results.exc.indx);
        results.sup.N_neurons = sum(results.sup.indx);
        results.sup.f_neurons = sum(results.sup.indx)/numel(results.sup.indx);
        results.resp.N_neurons = sum(results.resp.indx);
        results.resp.f_neurons = sum(results.resp.indx)/numel(results.resp.indx);
        
       
        %% extra stuff to align responses to preferred heading for the 1DAzi case
        for i=1:nunits
            s0 = thisdata(i).stim(ceil(end/2)+1); ds = diff(thisdata(i).stim); ds = ds(1);
            % find preferred stimulus
            pref_stim = thisdata(i).stats.stim.pref;
            % if pref_stim~=180, shift
            thisdata(i).rate_avg = circshift(thisdata(i).rate_avg,(s0 - pref_stim)/ds);
            thisdata(i).rate_pst = circshift(thisdata(i).rate_pst,(s0 - pref_stim)/ds);
            % duplicate 0 to 360
            thisdata(i).stim(end+1) = thisdata(i).stim(end)+ds; % stim
            thisdata(i).rate_avg(end+1) = thisdata(i).rate_avg(1); % tuning
            thisdata(i).rate_pst(end+1,:) = thisdata(i).rate_pst(1,:); % pst
        end
        
        %% store matrix of population responses to preferred direction (not smoothed)
%         for i=1:nunits
%             results.all.rate_pst.aligned(i,:) = thisdata(i).rate_pst(5,:);  %5 --> aligned to preferred heading
%         end
        %% population average psth
        for i=1:nunits, [nstim(i),ntime(i)] = size(thisdata(i).rate_pst); end
        for j=1:min(nstim)
            for i=1:nunits
                rate_pst(i,j,:) = smooth_pst(thisdata(i).rate_pst(j,1:min(ntime)),prs.dt,prs.tsmooth_pop);
                rate_pst0(i,j,:) = smooth_pst(thisdata0(i).rate_pst(1,1:min(ntime)),prs.dt,prs.tsmooth_pop);
                if j == 5
                    results.all.rate_pst.aligned(i,:) = rate_pst(i,j,:); % store matrix of population responses to preferred direction
                end
            end
            %             for i=1:nunits
            %                 rate_pst(i,j,:) = thisdata(i).rate_pst(j,1:min(ntime));
            %                 rate_pst0(i,j,:) = thisdata0(i).rate_pst(1,1:min(ntime));
%             end


            results.all.rate_pst.mu(j,:) = squeeze(mean(rate_pst(:,j,:)));
            results.all.rate_pst.sig(j,:) = squeeze(std(rate_pst(:,j,:)))/sqrt(nunits);
            results.all.rate_pst0.mu = squeeze(mean(rate_pst0(:,j,:)));
            results.all.rate_pst0.sig = squeeze(std(rate_pst0(:,j,:)))/sqrt(nunits);
            
            
            results.tuned.rate_pst.mu(j,:) = squeeze(mean(rate_pst(tuned,j,:)));
            results.tuned.rate_pst.sig(j,:) = squeeze(std(rate_pst(tuned,j,:)))/sqrt(sum(tuned));
            results.tuned.rate_pst0.mu = squeeze(mean(rate_pst0(tuned,j,:)));
            results.tuned.rate_pst0.sig = squeeze(std(rate_pst0(tuned,j,:)))/sqrt(sum(tuned));
            
            results.untuned.rate_pst.mu(j,:) = squeeze(mean(rate_pst(untuned,j,:)));
            results.untuned.rate_pst.sig(j,:) = squeeze(std(rate_pst(untuned,j,:)))/sqrt(sum(untuned));
            results.untuned.rate_pst0.mu = squeeze(mean(rate_pst0(untuned,j,:)));
            results.untuned.rate_pst0.sig = squeeze(std(rate_pst0(untuned,j,:)))/sqrt(sum(untuned));
            
            results.resp.rate_pst.mu(j,:) = squeeze(mean(rate_pst(resp,j,:)));
            results.resp.rate_pst.sig(j,:) = squeeze(std(rate_pst(resp,j,:)))/sqrt(sum(resp));
            results.resp.rate_pst0.mu = squeeze(mean(rate_pst0(resp,j,:)));
            results.resp.rate_pst0.sig = squeeze(std(rate_pst0(resp,j,:)))/sqrt(sum(resp));
            
            results.unresp.rate_pst.mu(j,:) = squeeze(mean(rate_pst(unresp,j,:)));
            results.unresp.rate_pst.sig(j,:) = squeeze(std(rate_pst(unresp,j,:)))/sqrt(sum(unresp));
            results.unresp.rate_pst0.mu = squeeze(mean(rate_pst0(unresp,j,:)));
            results.unresp.rate_pst0.sig = squeeze(std(rate_pst0(unresp,j,:)))/sqrt(sum(unresp));
            
            results.exc.rate_pst.mu(j,:) = squeeze(mean(rate_pst(exc,j,:)));
            results.exc.rate_pst.sig(j,:) = squeeze(std(rate_pst(exc,j,:)))/sqrt(sum(exc));
            results.exc.rate_pst0.mu = squeeze(mean(rate_pst0(exc,j,:)));
            results.exc.rate_pst0.sig = squeeze(std(rate_pst0(exc,j,:)))/sqrt(sum(exc));
            
            results.sup.rate_pst.mu(j,:) = squeeze(mean(rate_pst(sup,j,:)));
            results.sup.rate_pst.sig(j,:) = squeeze(std(rate_pst(sup,j,:)))/sqrt(sum(sup));
            results.sup.rate_pst0.mu = squeeze(mean(rate_pst0(sup,j,:)));
            results.sup.rate_pst0.sig = squeeze(std(rate_pst0(sup,j,:)))/sqrt(sum(sup));
            
            results.tuned_exc.rate_pst.mu(j,:) = squeeze(mean(rate_pst(tuned_exc,j,:)));
            results.tuned_exc.rate_pst.sig(j,:) = squeeze(std(rate_pst(tuned_exc,j,:)))/sqrt(sum(tuned_exc));
            results.tuned_exc.rate_pst0.mu = squeeze(mean(rate_pst0(tuned_exc,j,:)));
            results.tuned_exc.rate_pst0.sig = squeeze(std(rate_pst0(tuned_exc,j,:)))/sqrt(sum(tuned_exc));
            
            results.tuned_sup.rate_pst.mu(j,:) = squeeze(mean(rate_pst(tuned_sup,j,:)));
            results.tuned_sup.rate_pst.sig(j,:) = squeeze(std(rate_pst(tuned_sup,j,:)))/sqrt(sum(tuned_sup));
            results.tuned_sup.rate_pst0.mu = squeeze(mean(rate_pst0(tuned_sup,j,:)));
            results.tuned_sup.rate_pst0.sig = squeeze(std(rate_pst0(tuned_sup,j,:)))/sqrt(sum(tuned_sup));
        end
        results.all.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.tuned.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.untuned.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.resp.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.unresp.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.exc.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.sup.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.tuned_exc.rate_pst.time = thisdata(1).time(1:min(ntime));
        results.tuned_sup.rate_pst.time = thisdata(1).time(1:min(ntime));
        
%         % PCA of psth
%         tindx = thisdata(1).time>prs.tstim_on &
%         thisdata(1).time<prs.tstim_off+0.4; 
%         X=squeeze(mean(rate_pst(:,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.all.pca_pst = Y1; results.all.pca_var = Y3;
%         X=squeeze(mean(rate_pst(tuned,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.tuned.pca_pst = Y1; results.tuned.pca_var = Y3;
%         X=squeeze(mean(rate_pst(untuned,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.untuned.pca_pst = Y1; results.untuned.pca_var = Y3;
%         X=squeeze(mean(rate_pst(resp,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.resp.pca_pst = Y1; results.resp.pca_var = Y3;
%         X=squeeze(mean(rate_pst(unresp,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.unresp.pca_pst = Y1; results.unresp.pca_var = Y3;
%         X=squeeze(mean(rate_pst(exc,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.exc.pca_pst = Y1; results.exc.pca_var = Y3;
%         X=squeeze(mean(rate_pst(sup,:,tindx),2));
%         [Y1,~,Y3]=pca(X); results.sup.pca_pst = Y1; results.sup.pca_var = Y3;
        
        % population average tuning curves
        for i=1:nunits
            rate_avg(i,:) = [thisdata(i).rate_avg.mu]; 
            rate_avg0(i,:) = [thisdata0(i).rate_avg.mu];   
        end
        results.all.rate_avg.stim = thisdata(1).stim;
        results.all.rate_avg.mu = mean(rate_avg);
        results.all.rate_avg.sig = std(rate_avg)/sqrt(nunits);
        results.all.rate_avg0.mu = mean(rate_avg0);
        results.all.rate_avg0.sig = std(rate_avg0)/sqrt(nunits);
        
        results.tuned.rate_avg.stim = thisdata(1).stim;
        results.tuned.rate_avg.mu = mean(rate_avg(results.tuned.indx,:));
        results.tuned.rate_avg.sig = std(rate_avg(results.tuned.indx,:))/sqrt(sum(results.tuned.indx));
        results.tuned.rate_avg0.mu = mean(rate_avg0(results.tuned.indx,:));
        results.tuned.rate_avg0.sig = std(rate_avg0(results.tuned.indx,:))/sqrt(sum(results.tuned.indx));
        
        results.untuned.rate_avg.stim = thisdata(1).stim;
        results.untuned.rate_avg.mu = mean(rate_avg(results.untuned.indx,:));
        results.untuned.rate_avg.sig = std(rate_avg(results.untuned.indx,:))/sqrt(sum(results.untuned.indx));
        results.untuned.rate_avg0.mu = mean(rate_avg0(results.untuned.indx,:));
        results.untuned.rate_avg0.sig = std(rate_avg0(results.untuned.indx,:))/sqrt(sum(results.untuned.indx));

        results.resp.rate_avg.stim = thisdata(1).stim;
        results.resp.rate_avg.mu = mean(rate_avg(results.resp.indx,:));
        results.resp.rate_avg.sig = std(rate_avg(results.resp.indx,:))/sqrt(sum(results.resp.indx));;
        results.resp.rate_avg0.mu = mean(rate_avg0(results.resp.indx,:));
        results.resp.rate_avg0.sig = std(rate_avg0(results.resp.indx,:))/sqrt(sum(results.resp.indx));
        
        results.unresp.rate_avg.stim = thisdata(1).stim;
        results.unresp.rate_avg.mu = mean(rate_avg(results.unresp.indx,:));
        results.unresp.rate_avg.sig = std(rate_avg(results.unresp.indx,:))/sqrt(sum(results.unresp.indx));
        results.unresp.rate_avg0.mu = mean(rate_avg0(results.unresp.indx,:));
        results.unresp.rate_avg0.sig = std(rate_avg0(results.unresp.indx,:))/sqrt(sum(results.unresp.indx));
        
        results.exc.rate_avg.stim = thisdata(1).stim;
        results.exc.rate_avg.mu = mean(rate_avg(results.exc.indx,:));
        results.exc.rate_avg.sig = std(rate_avg(results.exc.indx,:))/sqrt(sum(results.exc.indx));
        results.exc.rate_avg0.mu = mean(rate_avg0(results.exc.indx,:));
        results.exc.rate_avg0.sig = std(rate_avg0(results.exc.indx,:))/sqrt(sum(results.exc.indx));
        
        results.sup.rate_avg.stim = thisdata(1).stim;
        results.sup.rate_avg.mu = mean(rate_avg(results.sup.indx,:));
        results.sup.rate_avg.sig = std(rate_avg(results.sup.indx,:))/sqrt(sum(results.sup.indx));
        results.sup.rate_avg0.mu = mean(rate_avg0(results.sup.indx,:));
        results.sup.rate_avg0.sig = std(rate_avg0(results.sup.indx,:))/sqrt(sum(results.sup.indx));
        
        results.tuned_exc.rate_avg.stim = thisdata(1).stim;
        results.tuned_exc.rate_avg.mu = mean(rate_avg(results.tuned_exc.indx,:));
        results.tuned_exc.rate_avg.sig = std(rate_avg(results.tuned_exc.indx,:))/sqrt(sum(results.tuned_exc.indx));
        results.tuned_exc.rate_avg0.mu = mean(rate_avg0(results.tuned_exc.indx,:));
        results.tuned_exc.rate_avg0.sig = std(rate_avg0(results.tuned_exc.indx,:))/sqrt(sum(results.tuned_exc.indx));
        
        results.tuned_sup.rate_avg.stim = thisdata(1).stim;
        results.tuned_sup.rate_avg.mu = mean(rate_avg(results.tuned_sup.indx,:));
        results.tuned_sup.rate_avg.sig = std(rate_avg(results.tuned_sup.indx,:))/sqrt(sum(results.tuned_sup.indx));
        results.tuned_sup.rate_avg0.mu = mean(rate_avg0(results.tuned_sup.indx,:));
        results.tuned_sup.rate_avg0.sig = std(rate_avg0(results.tuned_sup.indx,:))/sqrt(sum(results.tuned_sup.indx));
        
        % stim-response correlation
        y=rate_avg; x=repmat(thisdata(1).stim,[size(y,1),1]); [r,p] = corr(x(:),y(:));
        results.all.stats.corr_stimrate.r = r; results.all.stats.corr_stimrate.p = p;
        y=rate_avg(tuned,:); x=repmat(thisdata(1).stim,[size(y,1),1]); if ~isempty(y), [r,p] = corr(x(:),y(:)); else r = nan; p = nan; end
        results.tuned.stats.corr_stimrate.r = r; results.tuned.stats.corr_stimrate.p = p;
        y=rate_avg(untuned,:); x=repmat(thisdata(1).stim,[size(y,1),1]); if ~isempty(y), [r,p] = corr(x(:),y(:)); else r = nan; p = nan; end
        results.untuned.stats.corr_stimrate.r = r; results.untuned.stats.corr_stimrate.p = p;
        y=rate_avg(resp,:); x=repmat(thisdata(1).stim,[size(y,1),1]); if ~isempty(y), [r,p] = corr(x(:),y(:)); else r = nan; p = nan; end
        results.resp.stats.corr_stimrate.r = r; results.resp.stats.corr_stimrate.p = p;
        y=rate_avg(unresp,:); x=repmat(thisdata(1).stim,[size(y,1),1]); if ~isempty(y), [r,p] = corr(x(:),y(:)); else r = nan; p = nan; end
        results.unresp.stats.corr_stimrate.r = r; results.unresp.stats.corr_stimrate.p = p;
        y=rate_avg(exc,:); x=repmat(thisdata(1).stim,[size(y,1),1]); if ~isempty(y), [r,p] = corr(x(:),y(:)); else r = nan; p = nan; end
        results.exc.stats.corr_stimrate.r = r; results.exc.stats.corr_stimrate.p = p;
        if sum(sup)>0
            y=rate_avg(sup,:); x=repmat(thisdata(1).stim,[size(y,1),1]); [r,p] = corr(x(:),y(:));
        else
            r=nan; p=nan;
        end
        results.sup.stats.corr_stimrate.r = r; results.sup.stats.corr_stimrate.p = p;
        
        % discindx stats
        for i=1:nunits, d(i) = thisdata(i).stats.discindx; end
        results.all.discindx = d;
        results.unresp.discindx = d(unresp);
        results.resp.discindx = d(resp);
        results.tuned.discindx = d(tuned);
        results.untuned.discindx = d(untuned);
        results.tuned_exc.discindx = d(tuned_exc);
        results.tuned_sup.discindx = d(tuned_sup);
        if ~isempty(d(resp)),[results.resp.stats.discindx.p,results.resp.stats.discindx.h]=ranksum(d(resp),d(unresp));
        else
            results.resp.stats.discindx.p = nan;
            results.resp.stats.discindx.h = nan;
        end
        if ~isempty(d(tuned)) & ~isempty(d(resp & untuned)),[results.tuned.stats.discindx.p,results.tuned.stats.discindx.h]=ranksum(d(tuned),d(resp & untuned)); 
        else
            results.tuned.stats.discindx.p= nan; 
            results.tuned.stats.discindx.h = nan; 
        end
        results.exc.discindx = d(exc);
        if ~isempty(d(exc)) [results.exc.stats.discindx.p,results.exc.stats.discindx.h]=ranksum(d(exc),d(unresp)); 
        else
            results.exc.stats.discindx.p = nan; 
            results.exc.stats.discindx.h = nan; 
        end
        results.sup.discindx = d(sup);
        if sum(sup)>0
            [results.sup.stats.discindx.p,results.sup.stats.discindx.h]=ranksum(d(sup),d(unresp));
        else
            results.sup.stats.discindx.p = nan;
            results.sup.stats.discindx.h = nan;
        end
                
        % preferred/anti-preferred stimulus
        for i=1:nunits, stim_pref(i) = thisdata(i).stats.stim.pref; end
        for i=1:nunits, stim_antipref(i) = thisdata(i).stats.stim.antipref; end
        results.all.stim_pref = stim_pref; results.all.stim_antipref = stim_antipref;
        results.tuned.stim_pref = stim_pref(tuned); results.tuned.stim_antipref = stim_antipref(tuned);
        results.untuned.stim_pref = stim_pref(untuned); results.untuned.stim_antipref = stim_antipref(untuned);
        results.tuned_exc.stim_pref = stim_pref(tuned_exc); results.tuned_exc.stim_antipref = stim_antipref(tuned_exc);
        results.tuned_sup.stim_pref = stim_pref(tuned_sup); results.tuned_sup.stim_antipref = stim_antipref(tuned_sup);        
        results.resp.stim_pref = stim_pref(resp); results.resp.stim_antipref = stim_antipref(resp);
        results.unresp.stim_pref = stim_pref(unresp); results.unresp.stim_antipref = stim_antipref(unresp);
        results.exc.stim_pref = stim_pref(exc); results.exc.stim_antipref = stim_antipref(exc);
        results.sup.stim_pref = stim_pref(sup); results.sup.stim_antipref = stim_antipref(sup);
        
        % tuning
        for i=1:nunits, tuning_h(i) = thisdata(i).stats.flags.tuning; end
        for i=1:nunits, tuning_p(i) = thisdata(i).stats.pvals.tuning; end
        results.all.stats.tuning.h = tuning_h; results.all.stats.tuning.p = tuning_p;
        results.resp.stats.tuning.h = tuning_h(resp); results.resp.stats.tuning.p = tuning_p(resp);
        results.unresp.stats.tuning.h = tuning_h(unresp); results.unresp.stats.tuning.p = tuning_p(unresp);
        results.exc.stats.tuning.h = tuning_h(exc); results.exc.stats.tuning.p = tuning_p(exc);
        results.sup.stats.tuning.h = tuning_h(sup); results.sup.stats.tuning.p = tuning_p(sup);
        
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

        results.tuned.rise.t_on = t_on(tuned); results.tuned.rise.t_off = t_off(tuned);
        results.tuned.rise.type_on = type_on(tuned); results.tuned.rise.type_off = type_off(tuned);
        results.tuned.stats.rise.t_on.mu = nanmedian(results.tuned.rise.t_on);
        results.tuned.stats.rise.t_on.sig = nanstd(results.tuned.rise.t_on)/sqrt(sum(tuned));
        results.tuned.stats.rise.t_off.mu = nanmedian(results.tuned.rise.t_off);
        results.tuned.stats.rise.t_off.sig = nanstd(results.tuned.rise.t_off)/sqrt(sum(tuned));
        
        results.untuned.rise.t_on = t_on(untuned); results.untuned.rise.t_off = t_off(untuned);
        results.untuned.rise.type_on = type_on(untuned); results.tuned.rise.type_off = type_off(untuned);
        results.untuned.stats.rise.t_on.mu = nanmedian(results.untuned.rise.t_on);
        results.untuned.stats.rise.t_on.sig = nanstd(results.untuned.rise.t_on)/sqrt(sum(untuned));
        results.untuned.stats.rise.t_off.mu = nanmedian(results.untuned.rise.t_off);
        results.untuned.stats.rise.t_off.sig = nanstd(results.untuned.rise.t_off)/sqrt(sum(untuned));
        
        results.resp.rise.t_on = t_on(resp); results.resp.rise.t_off = t_off(resp);
        results.resp.rise.type_on = type_on(resp); results.resp.rise.type_off = type_off(resp);
        results.resp.stats.rise.t_on.mu = nanmedian(results.resp.rise.t_on);
        results.resp.stats.rise.t_on.sig = nanstd(results.resp.rise.t_on)/sqrt(sum(resp));
        results.resp.stats.rise.t_off.mu = nanmedian(results.resp.rise.t_off);
        results.resp.stats.rise.t_off.sig = nanstd(results.resp.rise.t_off)/sqrt(sum(resp));
        
        results.unresp.rise.t_on = t_on(unresp); results.unresp.rise.t_off = t_off(unresp);
        results.unresp.rise.type_on = type_on(unresp); results.unresp.rise.type_off = type_off(unresp);
        results.unresp.stats.rise.t_on.mu = nanmedian(results.unresp.rise.t_on);
        results.unresp.stats.rise.t_on.sig = nanstd(results.unresp.rise.t_on)/sqrt(sum(unresp));
        results.unresp.stats.rise.t_off.mu = nanmedian(results.unresp.rise.t_off);
        results.unresp.stats.rise.t_off.sig = nanstd(results.unresp.rise.t_off)/sqrt(sum(unresp));
        
        results.exc.rise.t_on = t_on(exc); results.exc.rise.t_off = t_off(exc);
        results.exc.rise.type_on = type_on(exc); results.exc.rise.type_off = type_off(exc);
        results.exc.stats.rise.t_on.mu = nanmedian(results.exc.rise.t_on);
        results.exc.stats.rise.t_on.sig = nanstd(results.exc.rise.t_on)/sqrt(sum(exc));
        results.exc.stats.rise.t_off.mu = nanmedian(results.exc.rise.t_off);
        results.exc.stats.rise.t_off.sig = nanstd(results.exc.rise.t_off)/sqrt(sum(exc));
        
        results.sup.rise.t_on = t_on(sup); results.sup.rise.t_off = t_off(sup);
        results.sup.rise.type_on = type_on(sup); results.sup.rise.type_off = type_off(sup);
        results.sup.stats.rise.t_on.mu = nanmedian(results.sup.rise.t_on);
        results.sup.stats.rise.t_on.sig = nanstd(results.sup.rise.t_on)/sqrt(sum(sup));
        results.sup.stats.rise.t_off.mu = nanmedian(results.sup.rise.t_off);
        results.sup.stats.rise.t_off.sig = nanstd(results.sup.rise.t_off)/sqrt(sum(sup));
        
        % latency
        clear t_on;
        for i=1:nunits, t_on(i) = thisdata(i).stats.events.latency.t_on(1); end
        
        results.all.latency.t_on = t_on;
        results.all.stats.latency.t_on.mu = nanmedian(results.all.latency.t_on);
        results.all.stats.latency.t_on.sig = nanstd(results.all.latency.t_on)/sqrt(nunits);
        
        results.tuned.latency.t_on = t_on(tuned);
        results.tuned.stats.latency.t_on.mu = nanmedian(results.tuned.latency.t_on);
        results.tuned.stats.latency.t_on.sig = nanstd(results.tuned.latency.t_on)/sqrt(sum(tuned));
        
        results.untuned.latency.t_on = t_on(untuned);
        results.untuned.stats.latency.t_on.mu = nanmedian(results.untuned.latency.t_on);
        results.untuned.stats.latency.t_on.sig = nanstd(results.untuned.latency.t_on)/sqrt(sum(untuned));
        
        results.resp.latency.t_on = t_on(resp);
        results.resp.stats.latency.t_on.mu = nanmedian(results.resp.latency.t_on);
        results.resp.stats.latency.t_on.sig = nanstd(results.resp.latency.t_on)/sqrt(sum(resp));
        
        results.unresp.latency.t_on = t_on(unresp);
        results.unresp.stats.latency.t_on.mu = nanmedian(results.unresp.latency.t_on);
        results.unresp.stats.latency.t_on.sig = nanstd(results.unresp.latency.t_on)/sqrt(sum(unresp));
        
        results.exc.latency.t_on = t_on(exc);
        results.exc.stats.latency.t_on.mu = nanmedian(results.exc.latency.t_on);
        results.exc.stats.latency.t_on.sig = nanstd(results.exc.latency.t_on)/sqrt(sum(exc));
        
        results.sup.latency.t_on = t_on(sup);
        results.sup.stats.latency.t_on.mu = nanmedian(results.sup.latency.t_on);
        results.sup.stats.latency.t_on.sig = nanstd(results.sup.latency.t_on)/sqrt(sum(sup));
        
        % peak time
        clear t_on;
        for i=1:nunits, t_on(i) = thisdata(i).stats.events.peak.t_on; end
        
        results.all.peak.t_on = t_on;
        results.all.stats.peak.t_on.mu = nanmedian(results.all.peak.t_on);
        results.all.stats.peak.t_on.sig = nanstd(results.all.peak.t_on)/sqrt(nunits);

        results.tuned.peak.t_on = t_on(tuned);
        results.tuned.stats.peak.t_on.mu = nanmedian(results.tuned.peak.t_on);
        results.tuned.stats.peak.t_on.sig = nanstd(results.tuned.peak.t_on)/sqrt(sum(tuned));
        
        results.untuned.peak.t_on = t_on(untuned);
        results.untuned.stats.peak.t_on.mu = nanmedian(results.untuned.peak.t_on);
        results.untuned.stats.peak.t_on.sig = nanstd(results.untuned.peak.t_on)/sqrt(sum(untuned));
        
        results.resp.peak.t_on = t_on(resp);
        results.resp.stats.peak.t_on.mu = nanmedian(results.resp.peak.t_on);
        results.resp.stats.peak.t_on.sig = nanstd(results.resp.peak.t_on)/sqrt(sum(resp));
        
        results.unresp.peak.t_on = t_on(unresp);
        results.unresp.stats.peak.t_on.mu = nanmedian(results.unresp.peak.t_on);
        results.unresp.stats.peak.t_on.sig = nanstd(results.unresp.peak.t_on)/sqrt(sum(unresp));
        
        results.exc.peak.t_on = t_on(exc);
        results.exc.stats.peak.t_on.mu = nanmedian(results.exc.peak.t_on);
        results.exc.stats.peak.t_on.sig = nanstd(results.exc.peak.t_on)/sqrt(sum(exc));
        
        results.sup.peak.t_on = t_on(sup);
        results.sup.stats.peak.t_on.mu = nanmedian(results.sup.peak.t_on);
        results.sup.stats.peak.t_on.sig = nanstd(results.sup.peak.t_on)/sqrt(sum(sup));
        
    case 'lfps'
        % population average psth
        results.time = thisdata(1).time;
        for i=1:nunits, [nstim(i),ntime(i)] = size(thisdata(i).wave_pst); end
        for j=1:min(nstim)
            for i=1:nunits
                wave_pst(i,j,:) = thisdata(i).wave_pst(j,1:min(ntime));
                wave_pst0(i,j,:) = thisdata0(i).wave_pst(1,1:min(ntime));
            end
            results.all.wave_pst.mu(j,:) = squeeze(mean(wave_pst(:,j,:)));
            results.all.wave_pst.sig(j,:) = squeeze(std(wave_pst(:,j,:)))/sqrt(nunits);
            results.all.wave_pst0.mu = squeeze(mean(wave_pst0(:,j,:)));
            results.all.wave_pst0.sig = squeeze(std(wave_pst0(:,j,:)))/sqrt(nunits);
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
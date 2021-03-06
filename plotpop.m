function plotpop(exp_name,unit_type,pop,pop_type,plottype)

switch exp_name
    case 'linearspeed'
        cond = {'ves','vis','com'};
    case 'angularspeed'
        cond = {'ves','vis'};
    case '1DAzi'
        cond = {'ves','vis','com'};
    case 'HD'
        cond = {'ves','vis','com'};
end

switch plottype
    case 'psth'
        switch exp_name
            case 'linearspeed'
                % gather data
                t = pop.(cond{1}).(pop_type).rate_pst.time;
                for i=1:length(cond)
                    r(i,:,:) = pop.(cond{i}).(pop_type).rate_pst.mu;
                    r0(i,:) = pop.(cond{i}).(pop_type).rate_pst0.mu;
                end
                % plot
                for i=1:size(r,1)
                    figure; hold on;
                    colorscale = 1/size(r,2):1/size(r,2):1;
                    plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                    for j=1:size(r,2)
                        plot(t,squeeze(r(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                    end
                    set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                        'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                    title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                        ' ' num2str(pop.unitname)]);
                end
            case 'angularspeed'
                % gather
                stim = pop.(cond{1}).(pop_type).rate_avg.stim;
                t = pop.(cond{1}).(pop_type).rate_pst.time;
                for i=1:length(cond)
                    r_ccw(i,:,:) = pop.(cond{i}).(pop_type).rate_pst.mu(stim<0,:);
                    r_cw(i,:,:) = pop.(cond{i}).(pop_type).rate_pst.mu(stim>0,:);
                    r0(i,:) = pop.(cond{i}).(pop_type).rate_pst0.mu;
                end
                % plot
                for i=1:size(r_ccw,1)
                    figure; hold on;
                    colorscale = 1/size(r_ccw,2):1/size(r_ccw,2):1;
                    plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                    for j=1:size(r_ccw,2)
                        plot(t,squeeze(r_ccw(i,j,:)),'Color',colorscale(6-j)*([1 2 3]==i),'linewidth',2);
                    end
                    set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                        'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                    title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                        ' ' num2str(pop.unitname) '_{ccw}']);
                    figure; hold on;
                    colorscale = 1/size(r_cw,2):1/size(r_cw,2):1;
                    plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                    for j=1:size(r_cw,2)
                        plot(t,squeeze(r_cw(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                    end
                    set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                        'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                    title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                        ' ' num2str(pop.unitname) '_{cw}']);
                end
                
            case '1DAzi'
                % gather % To plot for specific population type
                t = pop.(cond{1}).(pop_type).rate_pst.time;
                for i=1:length(cond)
                    r(i,:,:) = pop.(cond{i}).(pop_type).rate_pst.mu;
                    r0(i,:) = pop.(cond{i}).(pop_type).rate_pst0.mu;
                end
                pos_psth = ...
                    [.7 .4 .25 .25; ...
                    .7 .725 .25 .25; ...
                    .375 .725 .25 .25; ...
                    .05 .725 .25 .25; ...
                    .05 .4 .25 .25; ...
                    .05 .075 .25 .25; ...
                    .375 .075 .25 .25; ...
                    .7 .075 .25 .25];
                s = pop.(cond{1}).(pop_type).rate_avg.stim;
                
                %%% plot all directions (not useful)
                %                 for i=1:size(r,1)
                %                     %figure; hold on;
                %                     colorscale = 1/size(r,2):1/size(r,2):1;
                %                     %plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                %                     for j=1:8 %j=1:size(r,2)
                %                         axes('Position',pos_psth(j,:)); hold on
                %                         plot(t,squeeze(r(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                %                         plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                %                         set(gca,'xlim',[-0.1 3],'XTick',[0 0.1 0.5 2.5 2.9]-0.1,...
                %                         'XTickLabel',[0 0.1 0.5 2.5 2.9],'TickDir','Out','Fontsize',16);
                %                     end
                %                     %set(gca,'xlim',[-0.1 3],'XTick',[0 0.1 0.5 2.5 2.9]-0.1,...
                %                     %   'XTickLabel',[0 0.1 0.5 2.5 2.9],'TickDir','Out','Fontsize',16);
                %                     title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                %                         ' ' num2str(pop.unitname)]);
                %                 end
                
                
                % Plot only tuned direction
                for i=1:size(r,1)
                    figure; hold on;
                    colorscale = 1/size(r,2):1/size(r,2):1;
                    %plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                    for j=5 %tuned direction shifted %j=1:size(r,2)
                        plot(t,squeeze(r(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                        plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                        set(gca,'xlim',[-0.1 3],'XTick',[0 0.1 0.5 2.5 2.9]-0.1,...
                            'XTickLabel',[0 0.1 0.5 2.5 2.9],'TickDir','Out','Fontsize',16);
                    end
                    title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                        ' ' num2str(pop.unitname)]);
                    xlabel('Time (s)'); ylabel('Firing rate (spks/s)');
                end
                
                
                
                
                
                
                
            case 'HD'
                % gather
                t = pop.(cond{1}).(pop_type).rate_pst.time;
                for i=1:length(cond)
                    r(i,:,:) = pop.(cond{i}).(pop_type).rate_pst.mu;
                    r0(i,:) = pop.(cond{i}).(pop_type).rate_pst0.mu;
                end
                % plot
                for i=1:size(r,1)
                    figure; hold on;
                    colorscale = 1/size(r,2):1/size(r,2):1;
                    plot(t,squeeze(r0(i,:)),'Color','k','linewidth',2);
                    for j=1:size(r,2)
                        plot(t,squeeze(r(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                    end
                    set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                        'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                    title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                        ' ' num2str(pop.unitname)]);
                end
        end
        
    case 'tuning'
        % gather
        s = pop.(cond{1}).(pop_type).rate_avg.stim;
        for i=1:length(cond)
            r(i,:) = pop.(cond{i}).(pop_type).rate_avg.mu;
            e(i,:) = pop.(cond{i}).(pop_type).rate_avg.sig;
            r0(i) = pop.(cond{i}).(pop_type).rate_avg0.mu;
        end
        % plot
        figure; hold on;
        for i=1:size(r,1)
            errorbar(s,r(i,:),e(i,:)/2,'Color',[1 2 3]==i,'Linewidth',2);
        end
        hline(mean(r0),'--k');
        set(gca,'xlim',[min(s)-min(unique(diff(s))) max(s)+min(unique(diff(s)))],'XTick',s,...
            'XTickLabel',s*1e2,'TickDir','Out','Fontsize',16);
        title([num2str(pop.groupname) ' ' num2str(pop_type) ...
            ' ' num2str(pop.unitname)]);
        
    case 'tuning_1DAzi'
        % gather for specific pop_type
        s = pop.(cond{1}).(pop_type).rate_avg.stim;
        for i=1:length(cond)
            r(i,:) = pop.(cond{i}).(pop_type).rate_avg.mu;
            e(i,:) = pop.(cond{i}).(pop_type).rate_avg.sig;
            r0(i) = pop.(cond{i}).(pop_type).rate_avg0.mu;
        end
        % plot
        figure; hold on;
        for i=1:size(r,1)
            errorbar(s,r(i,:),e(i,:)/2,'Color',[1 2 3]==i,'Linewidth',2);
        end
        hline(mean(r0),'--k');
        set(gca,'xlim',[min(s)-unique(diff(s)) max(s)+unique(diff(s))],'XTick',s,...
            'XTickLabel',[-180 -135 -90 -45 0 45 90 135 180],'TickDir','Out','Fontsize',16);
        title([num2str(pop.groupname) ' ' num2str(pop_type) ...
            ' ' num2str(pop.unitname)]);
        xlabel('Azimuth (0=preferred)'); ylabel('Firing rate (spks/s)')
        
    case 'polar_1DAzi'
        s = pop.(cond{1}).(pop_type).rate_avg.stim;
        for i=1:length(cond)
            r(i,:) = pop.(cond{i}).(pop_type).rate_avg.mu;
            e(i,:) = pop.(cond{i}).(pop_type).rate_avg.sig;
            r0(i) = pop.(cond{i}).(pop_type).rate_avg0.mu;
        end
        
        for i=1:size(r,1)
            polarplot(s*pi/180,r(i,:),'Color',[1 2 3]==i,'Linewidth',2);
            hold on
        end
        
        
    case 'prefstim'
        % gather
        s = pop.(cond{1}).(pop_type).rate_avg.stim;
        for i=1:length(cond)
            ps = pop.(cond{i}).(pop_type).stim_pref;
            y(i,:) = hist(ps,s);
        end
        % plot
        figure; hold on;
        for i=1:size(y,1)
            plot(s,y(i,:),'Color',[1 2 3]==i,'Linewidth',2);
        end
        set(gca,'xlim',[min(s)-unique(diff(s)) max(s)+unique(diff(s))],'XTick',s,...
            'XTickLabel',s*1e2,'TickDir','Out','Fontsize',16);
        title([num2str(pop.groupname) ' ' num2str(pop_type) ...
            ' ' num2str(pop.unitname)]);
        
    case 'prefstim_1DAzi'
        % gather
        s = pop.(cond{1}).(pop_type).rate_avg.stim(1:8);
        for i=1:length(cond)
            ps = pop.(cond{i}).(pop_type).stim_pref;
            y(i,:) = hist(ps,s((1:8)));
        end
        % plot
        figure; hold on;
        for i=1:size(y,1)
            plot(s(1:8),y(i,:),'Color',[1 2 3]==i,'Linewidth',2);
        end
        set(gca,'xlim',[min(s)-unique(diff(s)) max(s)+unique(diff(s))],'XTick',s,...
            'XTickLabel',s,'TickDir','Out','Fontsize',16);
        title([num2str(pop.groupname) ' ' num2str(pop_type) ...
            ' ' num2str(pop.unitname)]);
        xlabel('Preferred direction (deg)'); ylabel('No. of neurons');
        
        
    case 'discindx'
        % gather
        dcntr = 0.05:0.05:1;
        %         for i=1:length(cond)
        %             d0 = pop.(cond{i}).untuned.discindx; med_d0(i) = median(d0);
        %             d = pop.(cond{i}).(pop_type).discindx; med_d(i) = median(d);
        %             y0(i,:) = hist(d0,dcntr); y(i,:) = hist(d,dcntr);
        %             %             z0(i,:) = ecdfhist(ecdf(d0),dcntr); z(i,:) = ecdfhist(ecdf(d),dcntr);
        %             %             hold on;
        %             %             [F,X] = ecdf(d0); plot(X,F,'Color',0.65*([1 2 3]==i),'Linewidth',2);
        %             %             [F,X] = ecdf(d); plot(X,F,'Color',([1 2 3]==i),'Linewidth',2);
        %             %             axis([0 1 0 1]);
        %         end
        
        for i=1:length(cond)
            %             d0 = pop.(cond{i}).unresp.discindx; med_d0(i) = nanmedian(d0); % median
            %             d = pop.(cond{i}).(pop_type).discindx; med_d(i) = nanmedian(d); % median
            d0 = pop.(cond{i}).unresp.discindx; med_d0(i) = nanmean(d0); %mean
            d = pop.(cond{i}).(pop_type).discindx; med_d(i) = nanmean(d); % mean
            
            y0(i,:) = hist(d0,dcntr); y(i,:) = hist(d,dcntr);
            %             z0(i,:) = ecdfhist(ecdf(d0),dcntr); z(i,:) = ecdfhist(ecdf(d),dcntr);
            %             hold on;
            %             [F,X] = ecdf(d0); plot(X,F,'Color',0.65*([1 2 3]==i),'Linewidth',2);
            %             [F,X] = ecdf(d); plot(X,F,'Color',([1 2 3]==i),'Linewidth',2);
            %             axis([0 1 0 1]);
            
        end
        
        
        % plot pdf
        for i=1:length(cond)
            figure; hold on;
            plot(dcntr,y0(i,:),'k','Linewidth',2);
            ha = area(dcntr,y(i,:),'FaceColor',[1 2 3]==i,'Linewidth',2);
            set(get(ha,'Children'),'FaceAlpha',0.5,'EdgeColor','none');
            set(gca,'xlim',[0 1],'XTick',[0 0.5 1],'XTickLabel',[0 0.5 1], 'yTick',[0 30 60], 'YTickLabel',[0 30 60],'TickDir','Out','Fontsize',16);
            vline(med_d0(i),'--k'); vline(med_d(i),'--r');
            title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                ' ' num2str(pop.unitname)]);
        end
    case 'discindx_1DAzi'
        % gather
        dcntr = 0.05:0.05:1;
        %         for i=1:length(cond)
        %             d0 = pop.(cond{i}).untuned.discindx; med_d0(i) = median(d0);
        %             d = pop.(cond{i}).(pop_type).discindx; med_d(i) = median(d);
        %             y0(i,:) = hist(d0,dcntr); y(i,:) = hist(d,dcntr);
        %             %             z0(i,:) = ecdfhist(ecdf(d0),dcntr); z(i,:) = ecdfhist(ecdf(d),dcntr);
        %             %             hold on;
        %             %             [F,X] = ecdf(d0); plot(X,F,'Color',0.65*([1 2 3]==i),'Linewidth',2);
        %             %             [F,X] = ecdf(d); plot(X,F,'Color',([1 2 3]==i),'Linewidth',2);
        %             %             axis([0 1 0 1]);
        %         end
        
        
        for i=1:length(cond)
            % all neurons
            d = pop.(cond{i}).(pop_type).discindx; med_d(i) = nanmedian(d);
            d_modality(i,:)=pop.(cond{i}).(pop_type).discindx;
            y(i,:) = hist(d,dcntr);
            %tuned
            d_tuned = pop.(cond{i}).tuned.discindx; med_d_tuned(i) = nanmedian(d_tuned);
            d_modality_tuned(i,:)=pop.(cond{i}).(pop_type).discindx;
            y_tuned(i,:) = hist(d_tuned,dcntr);
            
            %             z0(i,:) = ecdfhist(ecdf(d0),dcntr); z(i,:) = ecdfhist(ecdf(d),dcntr);
            %             hold on;
            %             [F,X] = ecdf(d0); plot(X,F,'Color',0.65*([1 2 3]==i),'Linewidth',2);
            %             [F,X] = ecdf(d); plot(X,F,'Color',([1 2 3]==i),'Linewidth',2);
            %             axis([0 1 0 1]);
            
        end
        
        
        
        % plot
        for i=1:length(cond)
            figure; hold on;
            plot(dcntr,y(i,:),'k','Linewidth',2); %all
            plot(dcntr,y_tuned(i,:),'k','Linewidth',2); %tuned
            ha = area(dcntr,y_tuned(i,:),'FaceColor',[1 2 3]==i,'Linewidth',2);
            set(get(ha,'Children'),'FaceAlpha',0.5,'EdgeColor','none');
            set(gca,'xlim',[0 1],'XTick',[0 0.5 1],'XTickLabel',[0 0.5 1],'TickDir','Out','Fontsize',16);
            vline(med_d(i),'--k'); vline(med_d_tuned(i),'--c');
            title([num2str(cond{i}) ' : ' num2str(pop.groupname) ' ' num2str(pop_type) ...
                ' ' num2str(pop.unitname)]);
            ylabel('No. of neurons'); xlabel('DDI');
        end
        
        
        
        
    case 'compare_rmax' % For angular speed
        units = pop{1};
        for i=1:length(units)
            rves(i,:) = [units(i).ves.rate_avg.mu] - units(i).null.rate_avg.mu;
            rvis(i,:) = [units(i).vis.rate_avg.mu] - units(i).null.rate_avg.mu;
        end
        vesresp = pop{2}.ves.resp.indx; vestuned = pop{2}.ves.tuned.indx;
        visresp = pop{2}.vis.resp.indx; vistuned = pop{2}.vis.tuned.indx;
        hold on; plot(0.01:100,0.01:100,'--k');
        scatter(rves(vesresp & ~vestuned,1),rves(vesresp & ~vestuned,end),'or');
        scatter(rves(vestuned,1),rves(vestuned,end),'or','filled');
        scatter(rvis(visresp & ~vistuned,1),rvis(visresp & ~vistuned,end),'og');
        scatter(rvis(vistuned,1),rvis(vistuned,end),'og','filled');
        set(gca,'XScale','log','YScale','log','TickDir','Out');
        axis([1 100 1 100]);
        xlabel('ccw'); ylabel('cw');
        
    case 'compare_modal'
        switch exp_name
            case 'linearspeed'
                units = pop{1};
                for i=1:length(units)
                    r_ves(i) = max(abs([units(i).ves.rate_avg.mu] - units(i).null.rate_avg.mu));
                    r_vis(i) = max(abs([units(i).vis.rate_avg.mu] - units(i).null.rate_avg.mu));
                    r_com(i) = max(abs([units(i).com.rate_avg.mu] - units(i).null.rate_avg.mu));
                end
                %indx
                resp_ves_exc = pop{2}.ves.exc.indx;
                resp_ves_sup = pop{2}.ves.sup.indx;
                resp_vis_exc = pop{2}.vis.exc.indx;
                resp_vis_sup = pop{2}.vis.sup.indx;
                
                ves_vis_exc = resp_ves_exc & resp_vis_exc;
                ves_vis_sup = resp_ves_sup & resp_vis_sup;
                
                
                
                figure; hold on;
                scatter(r_com(ves_vis_exc | 0)./r_vis(ves_vis_exc | 0), r_com(ves_vis_exc | 0)./r_ves(ves_vis_exc | 0),'ok','filled');
                scatter(r_com(ves_vis_sup | 0)./r_vis(ves_vis_sup | 0), r_com(ves_vis_sup | 0)./r_ves(ves_vis_sup | 0),'ok');
                plot(0.1:10,0.1:10,'--k');
                set(gca,'TickDir', 'out', 'xtick', [0 1 2 3], 'ytick', [0 1 2 3], 'FontSize',18);
                xlabel('(R_{com} - R_{0})/(R_{vis} - R_{0})'); ylabel('(R_{com} - R_{0})/(R_{ves} - R_{0})'); axis ([0 3 0 3]);
                %axis equal
                
                %[p,h] = signrank(r_com(ves_vis_exc | ves_vis_sup)./r_ves(ves_vis_exc | ves_vis_sup),r_com(ves_vis_exc | ves_vis_sup)./r_vis(ves_vis_exc | ves_vis_sup));
                [p,h] = ttest(r_com(ves_vis_exc | ves_vis_sup)./r_ves(ves_vis_exc | ves_vis_sup),r_com(ves_vis_exc | ves_vis_sup)./r_vis(ves_vis_exc | ves_vis_sup));
                %text(2,1,['p = ' num2str(p)]);
                
                % y axis histogram
                figure; hold on;
                hist(r_com(ves_vis_exc | ves_vis_sup)./r_ves(ves_vis_exc | ves_vis_sup),30);
                xlim ([0 3])
                set(gca,'ylim', [0 8], 'TickDir', 'out','ytick', [0 8],'xtick', [0 1 2 3],'XDir', 'reverse','FontSize', 18); xlabel('r_{ves}/r_{com}','FontSize', 18);
                vline(mean(r_com(ves_vis_exc | ves_vis_sup)./r_ves(ves_vis_exc | ves_vis_sup),'omitnan'));
                % x axis histogram
                figure; hold on;
                hist(r_com(ves_vis_exc | ves_vis_sup)./r_vis(ves_vis_exc | ves_vis_sup),30);
                xlim ([0 3])
                set(gca, 'ylim', [0 14],'TickDir', 'out','ytick', [0 14],'xtick', [0 1 2 3],'FontSize', 18); xlabel('r_{vis}/r_{com}','FontSize', 18);
                vline(mean(r_com(ves_vis_exc | ves_vis_sup)./r_vis(ves_vis_exc | ves_vis_sup),'omitnan'));
                % diagonal histogram
                figure; hold on;
                hist(r_com(ves_vis_exc | ves_vis_sup)./r_ves(ves_vis_exc | ves_vis_sup) - r_com(ves_vis_exc | ves_vis_sup)./r_vis(ves_vis_exc | ves_vis_sup),linspace(-3,3,25))
                xlabel('diagonal'); set(gca, 'TickDir', 'out', 'ytick', [0 10],'FontSize', 18); vline(0); xlim ([-2.423 2.423]);
                vline(mean(r_com(ves_vis_exc | ves_vis_sup)./r_ves(ves_vis_exc | ves_vis_sup) - r_com(ves_vis_exc | ves_vis_sup)./r_vis(ves_vis_exc | ves_vis_sup),'omitnan'));
                
                
                
            case 'angularspeed'
                units = pop{1};
                plot(0.1:100,0.1:100,'--k');
                set(gca,'XScale','Log'); set(gca,'YScale','Log');
                xlabel('r_{ves} (spk/s)'); ylabel('r_{vis} (spk/s)'); axis ([0 100 0 100]);
                for i=1:length(units)
                    r_ves(i) = max([units(i).ves.rate_avg.mu]) - units(i).null.rate_avg.mu;
                    r_vis(i) = max([units(i).vis.rate_avg.mu]) - units(i).null.rate_avg.mu;
                end
                resp_ves = pop{2}.ves.exc.indx; %pop{2}.ves.resp.indx;
                resp_vis = pop{2}.vis.exc.indx; %pop{2}.vis.resp.indx;
                figure; hold on;
                scatter(r_ves(resp_ves),r_vis(resp_ves),'or','filled');
                scatter(r_ves(resp_vis),r_vis(resp_vis),'og','filled');
                scatter(r_ves(resp_ves & resp_vis),r_vis(resp_ves & resp_vis),'oc','filled');
            case '1DAzi'
                units = pop{1};
                pop = pop{2};
                
                for  i=1:length(units)
                    r_ves(i) = max([units(i).ves.rate_avg.mu]); %- units(i).null.rate_avg.mu;
                    r_vis(i) = max([units(i).vis.rate_avg.mu]); %- units(i).null.rate_avg.mu;
                    r_com(i) = max([units(i).com.rate_avg.mu]); %- units(i).null.rate_avg.mu;
                    
                    r_ves_vis(i) = (r_ves(i)+r_vis(i))/2;
                    
                    aligned_r_ves(i)= max(pop.ves.all.rate_pst.aligned(i,:));
                    aligned_r_vis(i)= max(pop.vis.all.rate_pst.aligned(i,:));
                    aligned_r_com(i)= max(pop.com.all.rate_pst.aligned(i,:));
                end
                %                 resp_ves = pop{2}.ves.all.indx;
                %                 resp_vis = pop{2}.vis.all.indx;
                %                 resp_com = pop{2}.com.all.indx;
                
                
                
                % plot
                figure; hold on;
                scatter(r_ves./r_com,...
                    r_vis./r_com,'ok');
                plot(0.1:10,0.1:10,'--k');
                set(gca,'XScale','Log'); set(gca,'YScale','Log', 'TickDir', 'out');
                xlabel('r_{ves}/r_{com}'); ylabel('r_{vis}/r_{com}'); axis ([0.1 10 0.1 10]);
                [h,p] = ttest(r_ves./r_com,r_vis./r_com);
                text(2,1,['p = ' num2str(p)]);
                vline(1); hline(1);
                scatter(nanmean(r_ves(r_com~=0)./r_com(r_com~=0)), nanmean(r_vis(r_com~=0)./r_com(r_com~=0)), 'or', 'filled');
                
                % plot mean ves+vis vs com
                figure; hold on;
                scatter(r_ves_vis,r_com,'ok','filled');
                plot(0.1:1000,0.1:1000,'--k');
                set(gca,'XScale','Log'); set(gca,'YScale','Log', 'TickDir', 'out');
                xlabel('r_{vesVis}'); ylabel('r_{com}'); axis ([0.4 200 0.4 200]);
                [h,p] = ttest(r_com,r_ves_vis);
                text(2,100,['p = ' num2str(p)]);
                vline(1); hline(1);
                
                
                
            case 'HD'
                units = pop{1};
                for i=1:length(units)
                    r_ves(i) = max([units(i).ves.rate_avg.mu]); %- units(i).null.rate_avg.mu;
                    r_vis(i) = max([units(i).vis.rate_avg.mu]); %- units(i).null.rate_avg.mu;
                    r_com(i) = max([units(i).com.rate_avg.mu]); %- units(i).null.rate_avg.mu;
                end
                resp_ves = pop{2}.ves.resp.indx;
                resp_vis = pop{2}.vis.resp.indx;
                resp_com = pop{2}.com.resp.indx;
                figure; hold on;
                scatter(r_ves(resp_ves & resp_vis)./r_com(resp_ves & resp_vis),...
                    r_vis(resp_ves & resp_vis)./r_com(resp_ves & resp_vis),'ok','filled');
                plot(0.1:10,0.1:10,'--k');
                set(gca,'XScale','Log'); set(gca,'YScale','Log');
                xlabel('r_{ves}/r_{com}'); ylabel('r_{vis}/r_{com}'); axis ([0.4 3 0.4 3]);
                [p,h] = ranksum(r_ves(resp_ves & resp_vis)./r_com(resp_ves & resp_vis),r_vis(resp_ves & resp_vis)./r_com(resp_ves & resp_vis));
                text(2,1,['p = ' num2str(p)]);
        end
        
    case 'compare_var'  %beta
        units = pop{1};
        for  i=1:length(units)
            r_ves(i) = max([units(i).ves.rate_avg.mu]); %- units(i).null.rate_avg.mu;
            r_vis(i) = max([units(i).vis.rate_avg.mu]); %- units(i).null.rate_avg.mu;
            r_com(i) = max([units(i).com.rate_avg.mu]); %- units(i).null.rate_avg.mu;
        end
        r_ves_var_med= nanmedian(r_ves./(r_vis+r_ves));
        r_vis_var_med= nanmedian(r_vis./(r_vis+r_ves));
        
        
        figure;
        hist(r_ves./(r_vis+r_ves))
        title('ves/vis+ves')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 14, 'XTick', [0 0.5 1]);
        figure;
        hist(r_vis./(r_vis+r_ves));
        title('vis/vis+ves')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 14, 'XTick', [0 0.5 1]);
        
        
        
        
        
    case 'compare_ssi'
        pop = pop{2};
        ves_tuned_exc = pop.ves.tuned_exc.indx;
        vis_tuned_exc = pop.vis.tuned_exc.indx;
        both_tuned = ves_tuned_exc & vis_tuned_exc;
        figure; hold on;
        scatter(pop.ves.all.discindx(ves_tuned_exc),pop.vis.all.discindx(ves_tuned_exc),'or','filled');
        scatter(pop.ves.all.discindx(vis_tuned_exc),pop.vis.all.discindx(vis_tuned_exc),'og','filled');
        scatter(pop.ves.all.discindx(both_tuned),pop.vis.all.discindx(both_tuned),'oc','filled');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
        xlabel('vestibular SSI'); ylabel('visual SSI');
        set(gca, 'TickDir', 'out', 'XTick', [0 0.5 1], 'YTick', [0 0.5 1]);
        %             figure; hold on;
        %             scatter(pop.ves.all.discindx(ves_tuned),pop.com.all.discindx(ves_tuned),'or','filled');
        %             scatter(pop.ves.all.discindx(vis_tuned),pop.com.all.discindx(vis_tuned),'og','filled');
        %             scatter(pop.ves.all.discindx(both_tuned),pop.com.all.discindx(both_tuned),'oc','filled');
        %             axis([0 1 0 1]);
        %             figure; hold on;
        %             scatter(pop.vis.all.discindx(ves_tuned),pop.com.all.discindx(ves_tuned),'or','filled');
        %             scatter(pop.vis.all.discindx(vis_tuned),pop.com.all.discindx(vis_tuned),'og','filled');
        %             scatter(pop.vis.all.discindx(both_tuned),pop.com.all.discindx(both_tuned),'oc','filled');
        %             axis([0 1 0 1]);
        
    case 'compare_ssi_all'
        pop = pop{2};
        figure; hold on;
        scatter(pop.vis.all.discindx,pop.ves.all.discindx,'ok','filled');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
        xlabel('visual SSI'); ylabel('vestibular SSI');
        set(gca, 'TickDir', 'out', 'XTick', [0 0.5 1], 'YTick', [0 0.5 1], 'FontSize', 16);
        
    case 'compare_ssi_1DAzi'
        pop = pop{2};
        figure; hold on;
        scatter(pop.ves.all.discindx,pop.vis.all.discindx,'ok','filled');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
        xlabel('vestibular DDI'); ylabel('visual DDI');
        set(gca, 'TickDir', 'out', 'XTick', [0 0.5 1], 'YTick', [0 0.5 1],'FontSize', 16);
        [p,h] = ttest(pop.ves.all.discindx,pop.vis.all.discindx);
        text(0.8,0.2,['p = ' num2str(h)]);
        hold off
        figure; hold on;
        scatter(pop.ves.all.discindx,pop.com.all.discindx,'ok','filled');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
        xlabel('vestibular DDI'); ylabel('combined DDI');
        set(gca, 'TickDir', 'out', 'XTick', [0 0.5 1], 'YTick', [0 0.5 1],'FontSize', 16);
        [p,h] = ttest(pop.ves.all.discindx,pop.com.all.discindx);
        text(0.8,0.2,['p = ' num2str(h)]);
        hold off
        figure; hold on;
        scatter(pop.vis.all.discindx,pop.com.all.discindx,'ok','filled');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
        xlabel('visual DDI'); ylabel('combined DDI');
        set(gca, 'TickDir', 'out', 'XTick', [0 0.5 1], 'YTick', [0 0.5 1], 'FontSize', 16);
        [p,h] = ttest(pop.vis.all.discindx,pop.com.all.discindx);
        text(0.8,0.2,['p = ' num2str(h)]);
        hold off
        
        
    case 'corr_spatial'
        switch unit_type
            case {'multiunits','lfps'}
                % gather
                for i=1:length(pop)
                    r(i,:,:) = pop(i).r;
                    x(i,:) = pop(i).x;
                end
                r = squeeze(mean(r));
                x = mean(x);
                % plot
                figure;
                imagesc(x(:),x(:),r,[0 1]); colormap hot; colorbar;
                set(gca,'XTick',2:2:16,'TickDir','Out','Fontsize',16);
                figure; hold on; colors = 'rgmbk';
                for i=1:length(x)
                    plot(x,r(i,:),'Color',char((mod(i,4)~=0)*colors(5) + (mod(i,4)==0)*colors(ceil(i/4))),...
                        'LineWidth',0.5 + 2*(mod(i,4)==0));
                end
                set(gca,'xlim',[1 16],'XTick',2:2:16,'TickDir','Out','Fontsize',16);
                xlabel('channel #'); ylabel('correlation coefficient');
                corr2.x = x; corr2.r = r; corr2 = corr_ch2dist(corr2,0.1); dx = corr2.dx; r_mu = corr2.r_mu;
                figure; plot(dx,r_mu,'Color','k','Linewidth',3);
                set(gca,'xlim',[min(dx) max(dx)],'TickDir','Out','Fontsize',16);
                xlabel('inter-channel distance (mm)'); ylabel('correlation coefficient');
            case 'singleunits'
                % gather
                r_mu = pop.all.r;
                dx = pop.all.dx;
                figure; plot(dx,r_mu,'Color','k','Linewidth',3);
                set(gca,'xlim',[min(dx) max(dx)],'TickDir','Out','Fontsize',16);
                xlabel('inter-channel distance (mm)'); ylabel('correlation coefficient');
        end
    case 'corr_spatiotemporal'
        % gather
        for i=1:length(pop)
            r(i,:,:,:) = pop(i).r;
            x(i,:) = pop(i).x;
            t(i,:) = pop(i).t;
        end
        r = squeeze(mean(r)); x = mean(x); t = mean(t);
        plot3(repmat(x,[length(t) 1]),repmat(t,[length(x) 1])',squeeze(r(1,:,:))');
        set(gca,'xlim',[min(x) max(x)],'TickDir','Out','Fontsize',16);
        xlabel('channel #'); ylabel('time-lag (s)');
    case 'distance'
        switch unit_type
            case {'singleunits','multiunits'}
                dx = pop.all.dx;
                r = pop.all.r_nspk;
                hold on;
                errorbar(dx,r.mu,r.std,'Linewidth',2,'Color','k','Marker','o','LineStyle','none');
            case 'lfps'
                dx = pop.all.dx;
                r = pop.all.r;
                %                 rvis = pop.all.vis.r.mu;
                hold on;
                errorbar(dx,r.mu,r.std,'Linewidth',2,'Color','k','Marker','o','LineStyle','none');
                %                 plot(dx,rvis,'Linewidth',3,'Color','g');
        end
    case 'clustering'
        multiunits = pop{1};
        singleunits = pop{2};
        sessions = pop{3};
        for i=1:length(singleunits)
            ves_discindx1(i) = singleunits(i).ves.stats.discindx;
            ves_discindx2(i) = multiunits(i).ves.stats.discindx;
            vis_discindx1(i) = singleunits(i).vis.stats.discindx;
            vis_discindx2(i) = multiunits(i).vis.stats.discindx;
        end
        ves_discindx1 = ves_discindx1(~isnan(ves_discindx1));
        ves_discindx2 = ves_discindx2(~isnan(ves_discindx1));
        % discindx
        % vestibular
        npairs = length(ves_discindx1);
        for i=1:100
            indx = randperm(npairs); indx = indx(1:100);
            d(i)=corr(ves_discindx1(indx)',ves_discindx2(indx)');
            indx1 = randperm(npairs); indx1 = indx1(1:100);
            indx2 = randperm(npairs); indx2 = indx2(1:100);
            d0(i)=corr(ves_discindx1(indx1)',ves_discindx2(indx2)');
        end
        rves_d.mu = median(d); rves_d.std = iqr(d);
        rves_d0.mu = median(d0); rves_d0.std = iqr(d0);
        % visual
        npairs = length(vis_discindx1);
        for i=1:100
            indx = randperm(npairs); indx = indx(1:100);
            d(i)=corr(vis_discindx1(indx)',vis_discindx2(indx)');
            indx1 = randperm(npairs); indx1 = indx1(1:100);
            indx2 = randperm(npairs); indx2 = indx2(1:100);
            d0(i)=corr(vis_discindx1(indx1)',vis_discindx2(indx2)');
        end
        rvis_d.mu = median(d); rvis_d.std = iqr(d);
        rvis_d0.mu = median(d0); rvis_d0.std = iqr(d0);
        
        rves_s=[]; rves_n=[]; rves_s0=[]; rves_n0=[];
        rvis_s=[]; rvis_n=[]; rvis_s0=[]; rvis_n0=[];
        for i=1:length(sessions)
            if ~isempty(sessions(i).singleunits)
                rves_s = [rves_s sessions(i).singleunits.clustering.ves.r_sig];
                rves_n = [rves_n sessions(i).singleunits.clustering.ves.r_noise];
                rvis_s = [rves_s sessions(i).singleunits.clustering.vis.r_sig];
                rvis_n = [rves_n sessions(i).singleunits.clustering.vis.r_noise];
                rves_s0 = [rves_s0 sessions(i).singleunits.clustering_shuffled.ves.r_sig];
                rves_n0 = [rves_n0 sessions(i).singleunits.clustering_shuffled.ves.r_noise];
                rvis_s0 = [rves_s0 sessions(i).singleunits.clustering_shuffled.vis.r_sig];
                rvis_n0 = [rves_n0 sessions(i).singleunits.clustering_shuffled.vis.r_noise];
            end
        end
        rves_s = nanmean(rves_s); rves_n = nanmean(rves_n);
        rves_s0 = nanmean(rves_s0); rves_n0 = nanmean(rves_n0);
        rvis_s = nanmean(rvis_s); rvis_n = nanmean(rvis_n);
        rvis_s0 = nanmean(rvis_s0); rvis_n0 = nanmean(rvis_n0);
        
        hold on;
        errorbarxy(rves_d.mu,rves_d0.mu,rves_d.std/2,...
            rves_d0.std/2,rves_d.std/2,rves_d0.std/2,'r'); hold on;
        plot(rves_d.mu,rves_d0.mu,'sr','MarkerFaceColor','r');
        hold on; errorbarxy(rvis_d.mu,rvis_d0.mu,rvis_d.std/2,...
            rvis_d0.std/2,rvis_d.std/2,rvis_d0.std/2,'g','g'); hold on;
        plot(rvis_d.mu,rvis_d0.mu,'sg','MarkerFaceColor','g');
        axis([-.1 .4 -.1 .4]); hline(0,'k'); vline(0,'k');
        plot(0:.1:.4,0:.1:.4,'--k');
    case 'clustering_discindx'
        multiunits = pop{1};
        singleunits = pop{2};
        for i=1:length(singleunits)
            ves_discindx1(i) = singleunits(i).ves.stats.discindx;
            ves_tuned1(i) = singleunits(i).ves.stats.flags.tuning;
            ves_rmax1(i) = max([singleunits(i).ves.rate_avg.mu]) - singleunits(i).null.rate_avg.mu;
            ves_discindx2(i) = multiunits(i).ves.stats.discindx;
            ves_tuned2(i) = multiunits(i).ves.stats.flags.tuning;
            ves_rmax2(i) = max([multiunits(i).ves.rate_avg.mu]) - multiunits(i).null.rate_avg.mu;
            vis_discindx1(i) = singleunits(i).vis.stats.discindx;
            vis_tuned1(i) = singleunits(i).vis.stats.flags.tuning;
            vis_rmax1(i) = max([singleunits(i).vis.rate_avg.mu]) - singleunits(i).null.rate_avg.mu;
            vis_discindx2(i) = multiunits(i).vis.stats.discindx;
            vis_tuned2(i) = multiunits(i).vis.stats.flags.tuning;
            vis_rmax2(i) = max([multiunits(i).vis.rate_avg.mu]) - multiunits(i).null.rate_avg.mu;
        end
        ves_tuned = ves_tuned1 | ves_tuned2;
        vis_tuned = vis_tuned1 | vis_tuned2;
        
        figure; hold on;
        indx = ~isnan(ves_discindx1) & ~isnan(ves_discindx2) & ~ves_tuned;
        scatter(ves_discindx1(indx),ves_discindx2(indx),'.r');
        scatter(ves_discindx1(ves_tuned),ves_discindx2(ves_tuned),'.r');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
        figure; hold on;
        indx = ~isnan(vis_discindx1) & ~isnan(vis_discindx2) & ~vis_tuned;
        scatter(vis_discindx1(indx),vis_discindx2(indx),'.g');
        scatter(vis_discindx1(vis_tuned),vis_discindx2(vis_tuned),'.g');
        axis([0 1 0 1]); plot(0:1,0:1,'--k');
    case 'allpsth'
        switch unit_type
            case {'singleunits','multiunits'}
                units = pop{1};
                pop = pop{2};
                for i=1:length(cond)
                    indx = pop.(cond{i}).(pop_type).indx;
                    theseunits = units(indx); nunits = length(theseunits);
                    for k=1:nunits, nt(k) = length(theseunits(k).(cond{i}).time); end; nt = min(nt);
                    t = theseunits(i).(cond{i}).time(1:nt); clear r;
                    for j=1:length(theseunits)
                        r(j,:) = mean(theseunits(j).(cond{i}).rate_pst(:,1:nt));
                    end
                    t_on = pop.(cond{i}).(pop_type).rise.t_on;
                    type_on = pop.(cond{i}).(pop_type).rise.type_on;
                    t_off = pop.(cond{i}).(pop_type).rise.t_off;
                    type_off = pop.(cond{i}).(pop_type).rise.type_off;
                    t_on(strcmp(type_on,'sup')) = NaN;
                    t_off(strcmp(type_off,'sup')) = NaN;
                    
                    [r,t,t_on,t_off] = sort_risetimes(r,t,t_on,t_off,cond{i});
                    off_units = ~isnan(t_off); on_units = isnan(t_off);
                    % plot colormap
                    B = goodcolormap('wr');
                    figure; set(gcf,'Position',[100 200 300 300]);
                    hold on; colormap(B');
                    imagesc(t,1:nunits,r,[0,1]);
                    scatter(t_on,1:nunits,5,'k','filled');
                    scatter(t_off,1:nunits,5,'k','filled');
                    set(gca,'xlim',[-0.1 1.5],'ylim',[1 nunits],'XTick',[0 0.4 1.1],'YTick',[],...
                        'XTickLabel',[],'YTickLabel',[],'TickDir','Out','Fontsize',16);
                    % plot average
                    figure; set(gcf,'Position',[100 200 300 150]);
                    hold on;
                    plot(t,mean(r(off_units,:)),'k','Linewidth',3);
                    plot(t,mean(r(on_units,:)),'Color',[.7 .7 .7],'Linewidth',3);
                    set(gca,'xlim',[-0.1 1.5],'ylim',[-.2 1],'XTick',[0 0.4 1.1],'YTick',[],...
                        'XTickLabel',[],'YTickLabel',[],'TickDir','Out','Fontsize',16);
                end
        end
        
    case 'allpsth_1DAzi'
        switch unit_type
            case {'singleunits','multiunits'}
                units = pop{1};
                pop = pop{2};
                for i=1:length(cond)
                    theseunits = pop; nunits = 1:length(pop.(cond{i}).all.rate_pst.aligned(:,1));
                    for k=1:length(nunits), nt(k) = length(theseunits.(cond{i}).all.rate_pst.aligned(k,:)); end; nt = min(nt);
                    t = units(i).(cond{i}).time(1:nt); clear r;
                    for j=1:length(nunits)
                        r(j,:) = pop.(cond{i}).all.rate_pst.aligned(j,1:nt);
                    end
                    % normalize and sort
                    [r,t] = smooth_colormap_1DAzi(r,t,cond{i});
                    % plot colormap
                    B = goodcolormap('bwr');
                    figure; set(gcf,'Position',[100 200 300 300]);
                    %hold on; colormap(B');
                    hold on; colormap(pink);
                    imagesc(t,1:nunits,r,[0,1]);
                    %scatter(t_on,1:nunits,5,'k','filled');
                    %scatter(t_off,1:nunits,5,'k','filled');
                    set(gca,'xlim',[-0.1 3],'XTick',[0 0.1 0.5 2.5 2.9]-0.1,'ylim',[1 nunits(end)],...
                        'XTickLabel',[0 0.1 0.4 2.4 2.8],'YTickLabel',[],'TickDir','Out','Fontsize',16);
                    xlabel('Time (s)'); ylabel('Neuron');
                    title(cond(i))
                    
                    % plot average
                    figure; set(gcf,'Position',[100 200 300 150]);
                    hold on;
                    plot(t,nanmean(r(:,:)),'k','Linewidth',3);
                    set(gca,'xlim',[-0.1 3],'ylim',[-.2 1],'XTick',[0 0.1 0.5 2.5 2.9]-0.1,...
                        'XTickLabel',[],'YTickLabel',[],'TickDir','Out','Fontsize',16);
                end
        end
    case 'lfps'
        sessions = pop{1};
        t = sessions(1).lfps.time; nt = numel(t); nsessions = length(sessions);
        for i=1:length(cond)
            clear r;
            for j=1:nsessions
                v = sessions(j).lfps.(cond{i});
                r(j,:) = v(1:nt);
            end
            [r,t] = normalise_psth(r,t);
            t_on = pop{2}.(cond{i}).all.rise.t_on;
            t_on = reshape(t_on,[length(t_on)/nsessions nsessions]); t_on = nanmean(t_on);
            t_off = pop{2}.(cond{i}).all.rise.t_off;
            t_off = reshape(t_off,[length(t_off)/nsessions nsessions]); t_off = nanmean(t_off);
            % plot colormap
            B = goodcolormap('bwr');
            figure; set(gcf,'Position',[100 200 300 300]);
            hold on; colormap(B');
            imagesc(t,1:nsessions,r,[-1,1]);
            %                     scatter(t_on,1:nsessions,5,'k','filled');
            %                     scatter(t_off,1:nsessions,5,'k','filled');
            set(gca,'xlim',[-0.1 1.5],'ylim',[1 nsessions],'XTick',[0 0.4 1.1],'YTick',[],...
                'XTickLabel',[],'YTickLabel',[],'TickDir','Out','Fontsize',16);
            % plot average
            figure; set(gcf,'Position',[100 200 300 150]);
            hold on;
            plot(t,mean(r),'k','Linewidth',3);
            set(gca,'xlim',[-0.1 1.5],'ylim',[-1 0.2],'XTick',[0 0.4 1.1],'YTick',[],...
                'XTickLabel',[],'YTickLabel',[],'TickDir','Out','Fontsize',16);
        end
        
    case 'hist_times'
        %% peak response
        %gather tuned neurons
        for i=1:length(cond)
            peak_tuned{i,:} = pop.(cond{i}).tuned.peak.t_on;
        end
        
        % all neurons
        for i=1:length(cond)
            peak_all{i,:} = pop.(cond{i}).all.peak.t_on;
        end
        
        % peak tuned neurons
        figure; hold on;
        for i=1:length(cond)
            h1 = histfit(peak_tuned{i,:},10,'kernel');
            alpha(0.25);
            set(h1(1),'FaceColor',[1 2 3]==i, 'EdgeColor',[1 2 3]==i);
            set(h1(2),'Color',[1 2 3]==i);
            title('tuned neurons')
            endt op
            
            % peak all neurons
            figure; hold on;
            for i=1:length(cond)
                h1 = histfit(peak_all{i,:},10,'kernel');
                alpha(0);
                set(h1(1),'FaceColor',[1 2 3]==i, 'EdgeColor',[1 2 3]==i);
                set(h1(2),'Color',[1 2 3]==i);
            end
            set(gca,'xlim',[0 2.9], 'TickDir', 'out', 'FontSize', 18,'XTick',[0 0.1 0.5 2.5 2.9]);
            xlabel('Peak response time (s)'); ylabel('Nr. of neurons')
            %% latency
            %
            %         % all neurons
            %         for i=1:length(cond)
            %             latency_all{i,:} = pop.(cond{i}).all.latency.t_on;
            %         end
            %
            %          % latency all neurons
            %         figure; hold on;
            %         for i=1:length(cond)
            %             h1 = histfit(latency_all{i,:},10,'kernel');
            %             alpha(0);
            %             set(h1(1),'FaceColor',[1 2 3]==i, 'EdgeColor',[1 2 3]==i);
            %             set(h1(2),'Color',[1 2 3]==i);
            %         end
            %         set(gca,'xlim',[0 2.9], 'TickDir', 'out', 'FontSize', 18,'XTick',[0 0.1 0.5 2.5 2.9]);
            %         xlabel('Latency (s)'); ylabel('Nr. of neurons')
        end
        
    case 'neuronal_threshold'
        
        %gather neuronal threshold
        load('psycho_neuro.mat');
        ps = psycho_neuro;
        for u=1:length(ps)
            for i=1:length(cond)
                if numel(cell2mat(ps(u).Thresh_neu))==3
                    thresh(u,i) = cell2mat(ps(u).Thresh_neu(i));
                else
                    continue;
                end
            end
            
        end
        
        %plot threshold distributions for different conditions
        
        figure; hold on;
        histogram(thresh(:,1), 50, 'EdgeColor', 'r', 'DisplayStyle', 'stairs');
        histogram(thresh(:,2), 50, 'EdgeColor', 'g', 'DisplayStyle', 'stairs');
        histogram(thresh(:,3), 50, 'EdgeColor', 'b', 'DisplayStyle', 'stairs');
        set(gca, 'TickDir', 'out','xlim',([1.5 180]),'ylim',([0 20]),'FontSize', 18);
        xlabel('neuronal threshold'); ylabel('nr of neurons')
        
        % plot scatter ves vs com and vis vs com TODO remove >300
        ves = thresh(:,1); 
        vis = thresh(:,2); 
        com = thresh(:,3); 
        [h_ves, p_ves] = ttest(ves,com);
        [h_vis, p_vis] = ttest(vis,com);
        
        % plot 
        figure; hold on
        plot(ves, com, '.r', 'MarkerSize', 18);
        plot(vis, com, '.g','MarkerSize', 18);
        set(gca,'TickDir', 'out','xlim',([0 200]), 'ylim', [0 200], 'FontSize', 18)
        plot(0:200, 0:200, 'k')
        xlabel ('Vis | Ves'); ylabel('Com')
        title('neuronal threshold')
        
    case 'psycho_threshold'
        
        %gather neuronal threshold
        load('psycho_neuro.mat');
        ps = psycho_neuro;
        for u=1:length(ps)
            for i=1:length(cond)
                if numel(cell2mat(ps(u).Thresh_neu))==3
                    thresh(u,i) = cell2mat(ps(u).Thresh_psy(i));
                else
                    continue;
                end
            end
            
        end
        
        %plot threshold distributions for different conditions
        
        figure; hold on;
        histogram(thresh(:,1), 50, 'EdgeColor', 'r', 'DisplayStyle', 'stairs');
        histogram(thresh(:,2), 50, 'EdgeColor', 'g', 'DisplayStyle', 'stairs');
        histogram(thresh(:,3), 50, 'EdgeColor', 'b', 'DisplayStyle', 'stairs');
        set(gca, 'TickDir', 'out','xlim',([1.5 180]),'ylim',([0 20]),'FontSize', 18);
        xlabel('neuronal threshold'); ylabel('nr of neurons')
        
    case 'psychometric_all'
        
        load('psycho_neuro.mat');
        ps = psycho_neuro;
      
        h{1} = 'ro';  f{1} = 'r-';
        h{2} = 'gd';  f{2} = 'g-';
        h{3} = 'bs';  f{3} = 'b-';
        
        %get mean function
        for i=1:length(cond)
            for u=1:length(ps)
                unique_headings = ps(u).fit_data_psycho_cum{1}(:,1);
                xi = min(unique_headings) : 0.1 : max(unique_headings);
                fit_vals{i,u,:} = cum_gaussfit(ps(u).psy_perf{i},xi); 
            end
            
        end
        
        for i=1:length(cond)
            figure; hold on;
            for u=1:length(ps)
                unique_headings = ps(u).fit_data_psycho_cum{1}(:,1);
                xi = min(unique_headings) : 0.1 : max(unique_headings);
                beta = [0, 1.0];
                %   plot data in logarithmical space instead of linspace
                %plot(unique_headings, pCorrect(i,:), h{i}, xi, cum_gaussfit(psy_perf{i}, xi), f{i}, 'MarkerSize', 6, 'Linewidth', 1.5);
                %plot(unique_headings, ps(u).fit_data_psycho_cum{i}(:,2), h{i}, xi, cum_gaussfit(ps(u).psy_perf{i}, xi), f{i}, 'MarkerSize', 6, 'Linewidth', 1.5);
                plot(xi, cum_gaussfit(ps(u).psy_perf{i}, xi),'Linewidth', 1.5, 'Color', [1 2 3]==i)
                set(gca, 'TickDir', 'out', 'ylim',([0,1]), 'YTick', [0 0.5 1], 'FontSize', 16);
                xlabel('Heading Angles');
                ylabel('Rightward Choices');
                box off
            end
            title('Psychometric function all')
        end
        
    case 'neurometric_all'
        
        load('psycho_neuro.mat');
        ps = psycho_neuro;
        
        h{1} = 'ro';  f{1} = 'r-'; 
        h{2} = 'gd';  f{2} = 'g-';  
        h{3} = 'bs';  f{3} = 'b-'; 
        for i=1:length(cond)
            figure; hold on;
            counter =0; 
            for u=1:length(ps)
                
                unique_headings = ps(u).fit_data_neuro_cum{1}(:,1);
                xi = min(unique_headings) : 0.1 : max(unique_headings);
                %plot(neu_heading, ps(u).Neuro_correct{i}, h{i},xi,cum_gaussfit(ps(u).neu_perf{i},xi),f{i},'MarkerSize', 6, 'Linewidth', 1.5);
                if ~isempty(ps(u).neu_perf) & ~(numel(ps(u).neu_perf)<3)
                    plot(xi,cum_gaussfit(ps(u).neu_perf{i},xi),'Linewidth', 1.5, 'Color', [1 2 3]==i)
                    set(gca, 'TickDir', 'out', 'ylim',([0,1]), 'xlim',([-50,50]), 'YTick', [0 0.5 1], 'FontSize', 16);
                    xlabel('Heading Angles');
                    box off
                    xlabel('Heading Angles');
                    ylabel('Rightward Choices');
                    counter = counter+1; 
                end
            end
            title(['Neurometric function all n= ' num2str(counter)])
        end
        
    case 'choice_corr'
        
        % corr
        for i=1:length(cond)
            heading = pop.(cond{i}).all.heading;
            choice = pop.(cond{i}).all.choice;
            spk = pop.(cond{i}).all.nspk;
            for j=1:length(pop.ves.all.stim)
                [corr_spk_headRHO(i,j),corr_spk_headPVAL(i,j)] = corr(spk(:,j),heading(:,j));
                [corr_spk_choiceRHO(i,j),corr_spk_choicePVAL(i,j)] = corr(spk(:,j),choice(:,j));
            end
        end
        
        % plot spk_heading
        figure; hold on;
        for i=1:length(cond)
            plot(corr_spk_headRHO(i,:),'.','Color', [1 2 3]==i, 'MarkerSize', 12)
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        ylabel('corr'); xlabel('neuron')
        title('corr spk headRHO')
        
        
        % plot spk_choice
        figure; hold on;
        for i=1:length(cond)
            plot(corr_spk_choiceRHO(i,:),'.','Color', [1 2 3]==i, 'MarkerSize', 12)
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        ylabel('corr'); xlabel('neuron')
        title('corr spk choiceRHO')
        
        
        % hist spk_heading corr
        figure; hold on;
        histogram(corr_spk_headRHO(1,:),'EdgeColor', 'r', 'FaceColor', 'r');
        histogram(corr_spk_headRHO(2,:),'EdgeColor', 'g', 'FaceColor', 'g');
        histogram(corr_spk_headRHO(3,:),'EdgeColor', 'b', 'FaceColor', 'b');
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        xlabel('Corr'); ylabel('nr of neurons')
        title('spk head corr');
        
        % hist spk_choice corr
        figure; hold on;
        histogram(corr_spk_choiceRHO(1,:),'EdgeColor', 'r', 'FaceColor', 'r');
        histogram(corr_spk_choiceRHO(2,:),'EdgeColor', 'g', 'FaceColor', 'g');
        histogram(corr_spk_choiceRHO(3,:),'EdgeColor', 'b', 'FaceColor', 'b');
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        xlabel('Corr'); ylabel('nr of neurons')
        title('spk choice corr');
        
        % plot spks vs heading (108 vis, 115 com, 93 com)
        figure; hold on;
        for i=1:length(cond)
            spk = pop.(cond{i}).all.nspk(:,93);
            heading = pop.(cond{i}).all.heading(:,93);
            if strcmp(cond{i}, 'ves')
                spk_mu = pop.(cond{i}).all.rate_avg(93,:);
                unique_headings = unique(pop.ves.all.stim)';
                plot(heading, spk+0.5,'.', 'Color', [1 2 3]==i, 'MarkerSize', 12);
                plot(unique_headings, spk_mu, 'Color',[1 2 3]==i, 'LineWidth', 2);
            else
                spk_mu = pop.(cond{i}).all.rate_avg(93,:);
                unique_headings = unique(pop.ves.all.stim)';
                plot(heading, spk,'.', 'Color', [1 2 3]==i, 'MarkerSize', 12);
                plot(unique_headings, spk_mu, 'Color',[1 2 3]==i, 'LineWidth', 2);
            end
            set(gca, 'TickDir', 'out', 'FontSize', 18);
            xlabel('Heading (deg)'); ylabel('spks & firing rate (spks/s)')
        end
        
        % plot just mean spks and heading
        figure; hold on; 
        for i=1:length(cond)
            spk = pop.(cond{i}).all.nspk(:,93);
            heading = pop.(cond{i}).all.heading(:,93); 
            spk_mu = pop.(cond{i}).all.rate_avg(93,:);
            unique_headings = unique(pop.ves.all.stim)';
            plot(unique_headings, spk_mu, 'Color',[1 2 3]==i, 'LineWidth', 2);
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        xlabel('Heading (deg)'); ylabel('spks & firing rate (spks/s)')
        
        
        
    case 'partial_corr'
        
        % gather spk_choice
        for j=1:length(pop.ves.all.stim)
            for i=1:length(cond)
                spk_choice(i,j) = abs(pop.(cond{i}).all.parCorr(1,2,j));
            end
        end
        
        % gather spk_heading
        for j=1:length(pop.ves.all.stim)
            for i=1:length(cond)
                spk_heading(i,j) = abs(pop.(cond{i}).all.parCorr(1,3,j));
            end
        end
        
        % plot distributions
        % spk_choice
        figure; hold on;
        for i=1:length(cond)
            h1 = histfit(spk_choice(i,:),20,'kernel'); 
            set(h1(1),'FaceColor', [1 2 3] ==i, 'EdgeColor', [1 2 3] ==i);
            set(h1(2),'Color',[1 2 3] ==i);
            alpha(0.25)
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        xlabel('abs parCorr spk choice');ylabel('nr of neurons');
        
        % spk_heading
         figure; hold on;
        for i=1:length(cond)
            h1 = histfit(spk_heading(i,:),20,'kernel'); 
            set(h1(1),'FaceColor', [1 2 3] ==i, 'EdgeColor', [1 2 3] ==i);
            set(h1(2),'Color',[1 2 3] ==i);
            alpha(0.25)
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        xlabel('abs parCorr spk heading');ylabel('nr of neurons');
        
         %% plot distributions no fit
        % spk_choice
        figure; hold on;
        for i=1:length(cond)
            hist(spk_choice(i,:),20);
            [h1,x1] = hist(spk_choice(i,:),20);
            plot(x1,h1, 'Color', [1 2 3]==i,'LineWidth', 3);
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        xlabel('abs parCorr spk choice');ylabel('nr of neurons');
        
        % spk_heading
        figure; hold on;
        for i=1:length(cond)
            hist(spk_heading(i,:),20);
            [h1,x1] = hist(spk_heading(i,:),20);
            plot(x1,h1, 'Color', [1 2 3]==i, 'LineWidth', 3);
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        xlabel('abs parCorr spk heading');ylabel('nr of neurons');
        
    case 'partial_corr_thresh'
        % gather parCorr choice_spk
        for j=1:length(pop.ves.all.stim)
            for i=1:length(cond)
                spk(i,j) = pop.(cond{i}).all.parCorr(1,2,j);
            end
        end
        
        %gather neuronal threshold  
        load('psycho_neuro.mat');
        ps = psycho_neuro;
        for u=1:length(ps)
            for i=1:length(cond)
                if numel(cell2mat(ps(u).Thresh_neu))==3
                    thresh(u,i) = cell2mat(ps(u).Thresh_neu(i));
                else
                    continue;
                end
            end
            
        end
        
        % plot neuronal thresholds vs partial choice corr
        figure; hold on;
        for i=1:length(cond)
            plot(spk(i,:), thresh(:,i),'.', 'Color', [1 2 3] == i, 'MarkerSize', 12); 
        end
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        xlabel('corr'); ylabel('neuronal threshold')
end
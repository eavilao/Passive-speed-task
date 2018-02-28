function plotunit(unit,unit_num,exp_name,plottype)

monk_id = unit.monk_id;
session_id = unit.session_id;
channel_no = unit.channel_no;

switch exp_name
    case 'linearspeed'
        cond = {'ves','vis','com'};
    case 'angularspeed'
        cond = {'ves','vis'};
    case '1DAzi'
        cond = {'ves','vis','com'};
        %cond = {'ves','vis'};
    case 'HD'
        cond = {'ves','vis','com'};
end

switch exp_name
    case 'linearspeed'
        switch plottype
            case 'psth'
                % gather
                t = unit.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:) = mean(unit.(cond{i}).rate_pst);
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    plot(t,r(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num)]);
            case 'tuning'
                % gather
                s = unit.(cond{1}).stim;
                for i=1:length(cond)
                    r(i,:) = [unit.(cond{i}).rate_avg.mu];
                    e(i,:) = [unit.(cond{i}).rate_avg.sig];
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    errorbar(s,r(i,:),e(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[min(s)-unique(diff(s)) max(s)+unique(diff(s))],'XTick',s,...
                    'XTickLabel',s*1e2,'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num)]);
            case 'psth_sig'
                %gather
                t = unit.(cond{1}).time;
                pos_psth = ...
                    [.15 .4 .15 .2; ...
                    .35 .4 .15 .2; ...
                    .55 .4 .15 .2; ...
                    .75 .4 .15 .2;];
                axes('Position',pos_psth(2,:)); hold on;
                for i=1;
                    for l=1:4
                        shadedErrorBar(t,unit.(cond{i}).rate_pst(l,:),unit.(cond{i}).rate_pst_sig(l,:),...
                            {'Color',[0 0 0] + 0.25*(l-1)*[1 0 0]});
                    end
                end
                yl = ylim; ylim([0 yl(2)]); xlim([-0.1 1.5]); yl2 = ylim;
                
                axes('Position',pos_psth(3,:)); hold on;
                for i=2;
                    for l=1:4
                        shadedErrorBar(t,unit.(cond{i}).rate_pst(l,:),unit.(cond{i}).rate_pst_sig(l,:),...
                            {'Color',[0 0 0] + 0.25*(l-1)*[0 1 0]});
                    end
                end
                yl = ylim; ylim([0 yl(2)]); xlim([-0.1 1.5]); yl3 = ylim;
                
                axes('Position',pos_psth(4,:)); hold on;
                for i=3;
                    for l=1:4
                        shadedErrorBar(t,unit.(cond{i}).rate_pst(l,:),unit.(cond{i}).rate_pst_sig(l,:),...
                            {'Color',[0 0 0] + 0.25*(l-1)*[0 0 1]});
                    end
                end
                yl = ylim; ylim([0 yl(2)]); xlim([-0.1 1.5]); yl4 = ylim;
                ymax = ceil(max(max([yl2;yl3;yl4]))/10)*10;
                ymaxYl2 = ceil(max(max(yl2))/10)*10;
                ax = get(gcf,'Children');
                set(ax(1),'Ylim',[0 ymax],'YTick',[0 ymax/2 ymax],'XTick',[-0.1 0 0.4 1.1 1.5],'TickDir','Out','XTickLabel',{'0',[],[],[], '1.6'}); %vline(0.4,'--k'); vline(1.1,'--k');
                set(ax(2),'Ylim',[0 ymax],'YTick',[0 ymax/2 ymax],'XTick',[-0.1 0 0.4 1.1 1.5],'TickDir','Out','XTickLabel',{'0',[],[],[], '1.6'}); %vline(0.4,'--k'); vline(1.1,'--k');
                set(ax(3),'Ylim',[0 ymax],'YTick',[0 ymax/2 ymax],'XTick',[-0.1 0 0.4 1.1 1.5],'TickDir','Out','XTickLabel',{'0',[],[],[], '1.6'}); %vline(0.4,'--k'); vline(1.1,'--k');
                
        end
    case 'angularspeed'
        switch plottype
            case 'psth'
                % gather
                t = unit.(cond{1}).time; s = unit.(cond{1}).stim;
                for i=1:length(cond)
                    r_ccw(i,:) = mean(unit.(cond{i}).rate_pst(s<0,:));    % ccw speeds
                    r_cw(i,:) = mean(unit.(cond{i}).rate_pst(s>0,:));     % cw speeds
                end
                % plot cw
                figure;hold on;
                for i=1:size(r_cw,1)
                    plot(t,r_cw(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                yticks = get(gca,'YTick'); ylim([0 max(abs(yticks))]);
                set(gca,'YTickLabel',abs(get(gca,'YTick')));
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num) 'cw']);
                % plot ccw
                figure; hold on;
                for i=1:size(r_ccw,1)
                    plot(t,r_ccw(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                yticks = get(gca,'YTick'); ylim([0 max(abs(yticks))]);
                set(gca,'YTickLabel',abs(get(gca,'YTick')));
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num), 'ccw']);
                
            case 'tuning'
                % gather
                s = unit.(cond{1}).stim;
                for i=1:length(cond)
                    r(i,:) = [unit.(cond{i}).rate_avg.mu];
                    e(i,:) = [unit.(cond{i}).rate_avg.sig];
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    errorbar(s,r(i,:),e(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[min(s)-unique(diff(s(s<0))) max(s)+unique(diff(s(s>0)))],'XTick',s,...
                    'XTickLabel',s,'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num)]);
            case 'psth_sig'
                %gather
                t = unit.(cond{1}).time;
                pos_psth(1,:,:) = ...
                    [.15 .6 .15 .2; ...
                    .35 .6 .15 .2; ...
                    .55 .6 .15 .2];
                pos_psth(2,:,:) = ...
                    [.15 .2 .15 .2; ...
                    .35 .2 .15 .2; ...
                    .55 .2 .15 .2];
                axes('Position',pos_psth(1,1,:)); hold on;
                for i=1;
                    for l=1:5 % anticlockwise
                        shadedErrorBar(t,-(unit.(cond{i}).rate_pst(l,:)),-(unit.(cond{i}).rate_pst_sig(l,:)),...
                            {'Color',[0 0 0] + 0.25*(5-l)*[1 0 0]});
                    end
                    for l=6:10 % clockwise
                        shadedErrorBar(t,unit.(cond{i}).rate_pst(l,:),unit.(cond{i}).rate_pst_sig(l,:),...
                            {'Color',[0 0 0] + 0.25*(l-6)*[1 0 0]});
                    end
                end
                ylim([-20 20]);
                xlim([-0.1 1.5]);
                
                
                axes('Position',pos_psth(2,2,:)); hold on;
                for i=2;
                    for l=1:5 % anticlockwise
                        shadedErrorBar(t,-(unit.(cond{i}).rate_pst(l,:)),-(unit.(cond{i}).rate_pst_sig(l,:)),...
                            {'Color',[0 0 0] + 0.25*(5-l)*[0 1 0]});
                    end
                    for l=6:10 % clockwise
                        shadedErrorBar(t,unit.(cond{i}).rate_pst(l,:),unit.(cond{i}).rate_pst_sig(l,:),...
                            {'Color',[0 0 0] + 0.25*(l-6)*[0 1 0]});
                    end
                    
                end
                ylim([-20 20]);
                xlim([-0.1 1.5]);
                
        end
    case '1DAzi'
        switch plottype
            case 'psth'
                % gather
                t = unit.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:) = mean(unit.(cond{i}).rate_pst);
                    r_tun(i,:) = [unit.(cond{i}).rate_avg.mu];
                    e_tun(i,:) = [unit.(cond{i}).rate_avg.sig];
                end
                
                
                r_tun(1,end+1)= r_tun(1,1);  %Quick dirty trick...to remove extra-repeated direction
                r_tun(2,9)= r_tun(2,1);
                r_tun(3,9)= r_tun(3,1);
                r_tun0 = repmat(unit.null.rate_avg.mu,1,9);
                
                
                pos_tuning = [.35 .375 .3 .3];
                pos_psth = ...
                    [.7 .4 .25 .25; ...
                    .7 .725 .25 .25; ...
                    .375 .725 .25 .25; ...
                    .05 .725 .25 .25; ...
                    .05 .4 .25 .25; ...
                    .05 .075 .25 .25; ...
                    .375 .075 .25 .25; ...
                    .7 .075 .25 .25];
                s = unit.(cond{1}).stim;
                s_tun = s; s_tun(end+1)=360;
                
                %get ylimits
                for i=1:size(r,1)
                    ybounds(i,:) = max(max(unit.(cond{i}).rate_pst(:,:)));
                end
                ybounds = ceil(max(abs(ybounds)));
                
                % Polar plot
                for i=1:size(r,1)
                    polarplot(s_tun*pi/180,r_tun(i,:),'Color',[1 2 3]==i,'Linewidth',2);
                    set(gca, 'Position', pos_tuning, 'ThetaTick', [0    45    90   135   180   225   270   315 ], 'FontSize', 12, 'Box', 'off')
                    set(gcf, 'Position', [275 49 1219 948])
                    hold on
                end
                %polarplot(s_tun*pi/180 ,r_tun0(1,:),'Color','k','Linewidth',2);
                polarplot(s_tun(1):0.01:s_tun(end),repelem(r_tun0(1),size(s_tun(1):0.01:s_tun(end),2)),'Color','k','Linewidth',2);
                
                % Plot - every direction separately
                %hFig = figure; set(gcf,'Position',[50 50 750 750]);
                %hFig = figure; set(gcf,'Position',[331 49 1181 948]);
                
                for l=1:length(unit.(cond{1}).stim)
                    axes('Position',pos_psth(l,:)); hold on
                    for i=1:size(r,1)
                        plot(t,unit.(cond{i}).rate_pst(l,:),'Color',i==[1 2 3], 'linewidth',2)
                        set(gca,'xlim',[-0.1 3],'XTick',[0 0.5 2.5 2.9]-0.1,...
                            'XTickLabel',[0 0.5 2.5 2.9],'ylim', [0 ybounds],'YTick',[0 ybounds/2 ybounds], 'TickDir','Out','Fontsize',18);
                        %                         set(gca,'xlim',[-0.1 3],'XTick',[0 0.1 0.5 2.5 2.9]-0.1,...
                        %                             'XTickLabel',[0 0.1 0.5 2.5 2.9],'ylim', [0 52],'YTick',[0 26 52], 'TickDir','Out','Fontsize',18);
                        if l == 3
                            title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num)]);
                        end
                        
                        vline(0,'g')
                    end
                end
                
            case 'tuning'
                % gather
                s = unit.(cond{1}).stim;
                for i=1:length(cond)
                    r(i,:) = [unit.(cond{i}).rate_avg.mu];
                    e(i,:) = [unit.(cond{i}).rate_avg.sig];
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    errorbar(s,r(i,:),e(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[min(s)-unique(diff(s)) max(s)+unique(diff(s))],'XTick',s,...
                    'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num)]);
                
        end
    case 'HD'
        switch plottype
            case 'psth'
                % gather
                t = unit.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:) = mean(unit.(cond{i}).rate_pst);
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    plot(t,r(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_u' num2str(unit_num)]);
            case 'psychometric'
                %gather
                unique_headings = unique(unit.ves.headings);
                for i=1:length(cond)
                    all_headings{i,:} = unit.(cond{i}).headings;
                    choice{i,:} = unit.(cond{i}).choice;
                end
                
                for i=1:length(cond)
                    for j=1:length(unique_headings)
                        id{j,:} = all_headings{i} == unique_headings(j);
                        nTrials(i,j) = sum(id{j,:});
                        nCorrect(i,j)= sum(all_headings{i} == unique_headings(j) & choice{i} == 0);
                        pCorrect(i,j) = nCorrect(i,j)/nTrials(i,j);
                        if unique_headings(j)<0
                            pCorrect(i,j)= 1-pCorrect(i,j);
                        end
                         fit_data_psycho_cum{i}(j,1)=unique_headings(j); 
                         fit_data_psycho_cum{i}(j,2)=pCorrect(i,j);
                         fit_data_psycho_cum{i}(j,3)=nTrials(i,j);
                    end
                end
                
%                 Plot data points only
%                 figure; hold on;
%                 for i=1:length(cond)
%                     plot(headings,pCorrect(i,:),'.','Color',i==[1 2 3],'MarkerSize',15)
%                     set(gca, 'TickDir','Out','Fontsize',16, 'ylim',([0 1]), 'yTick', ([0 0.5 1])); 
%                     box off
%                     ylabel('Proportion rightward');
%                     xlabel('Heading direction (deg)')
%                 end
                 
                 % fit psychometric function using Wichman's MLE method to estimate threshold and bias(same as TEMPO GUI)
                 for i=1:length(cond)
                     wichman_psy = pfit(fit_data_psycho_cum{i},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
                     Thresh_psy{i} = wichman_psy.params.est(2);
                     Bias_psy{i} = wichman_psy.params.est(1);
                     psy_perf{i} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
                 end
                 
                 %plot psychometric with fit 
                 h{1} = 'ro';  f{1} = 'r-';  g{1} = 'ro-';
                 h{2} = 'gd';  f{2} = 'g-';  g{2} = 'gd-';
                 h{3} = 'bs';  f{3} = 'b-';  g{3} = 'bs-';
                 figure; hold on;
                 for i=1:length(cond)
                     xi = min(unique_headings) : 0.1 : max(unique_headings);
                     beta = [0, 1.0];
                     %   plot data in logarithmical space instead of linspace
                     plot(unique_headings, pCorrect(i,:), h{i}, xi, cum_gaussfit(psy_perf{i}, xi), f{i}, 'MarkerSize', 6, 'Linewidth', 1.5);
                     set(gca, 'TickDir', 'out', 'ylim',([0,1]), 'YTick', [0 0.5 1], 'FontSize', 16);
                     xlabel('Heading Angles');
                     ylabel('Rightward Choices');
                     box off
                 end
                 title('Psychometric function')
                 annotation('textbox', [0.55 0.33 0.4 0.1], 'string', ['Threshold (1=Ves, 2=Vis, 3=Com)= ' num2str([Thresh_psy{:}])])
                 annotation('textbox',  [0.55 0.2 0.4 0.1], 'string', ['Bias (1=Ves, 2=Vis, 3=Com)= ' num2str([Bias_psy{:}])])
                 
                 % neurometric function
                 %gather
                 for i=1:length(cond)
                     for j=1:length(unique_headings)
                     r{i,j} = [unit.(cond{i}).nspk(j)'];
                     r_avg(i,:) = [unit.(cond{i}).rate_avg.mu]; 
                     end
                 end
                 
                
    
                 % decide whether ves and vis is congruent tuning. Fit line by linear
                 % regression first and compare the sign of each condition to decide whether
                 % congruent or opposite, this is used to check whether congruent cells lead
                 % to better neuronal performance in combined condition, and vice versa
                 for i=1:length(cond)
                     [rr,pp] = corrcoef(unique_headings,r_avg(i,:));
                     line_re{i} = rr(1,2);
                     line_p{i} = pp(1,2);
                 end
                 
                 % calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
                 for i=1:length(cond)
                     for j=1:length(unique_headings)-1
                         trials_n = all_headings{i} == unique_headings(j);
                         fit_data_neuro_cum{i}(j,3) = sum(trials_n);
                         if j < (1+length(unique_headings))/2
                             Neuro_correct{i}(j) =  rocN(r{i,length(unique_headings)-j+1}{:},r{i,j}{:},100 );
                         else
                             Neuro_correct{i}(j) =  rocN(r{i,length(unique_headings)-j}{:}, r{i,(j+1)}{:},100);
                         end
                         if line_re{i} > 0
                             Neuro_correct{i}(j) = 1 - Neuro_correct{i}{j};
                         end
                     end
                 end
                 
                 for i=1:length(cond)
                     fit_data_neuro_cum{i}(:,1) = unique_headings(unique_headings~=0);
                     fit_data_neuro_cum{i}(:,2) = Neuro_correct{i}(:);
                     wichman_neu = pfit(fit_data_neuro_cum{i}(2:7,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
                     Thresh_neu{i} = wichman_neu.params.est(2);
                     % negative and positive infinite value means flat tuning
                     if Thresh_neu{i}<0 | Thresh_neu{i}> 300
                         Thresh_neu{i} = 300;
                         wichman_neu.params.est(2) = 300;
                     end
                     Bias_neu{i} = wichman_neu.params.est(1);
                     neu_perf{i} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
                 end
                 
                 %plot neurometric function
                 h{1} = 'ro';  f{1} = 'r-';  g{1} = 'ro-';
                 h{2} = 'gd';  f{2} = 'g-';  g{2} = 'gd-';
                 h{3} = 'bs';  f{3} = 'b-';  g{3} = 'bs-';
                 
                 figure; hold on;
                 for i=1:length(cond)
                     neu_heading = unique_headings(unique_headings~=0);
                     xi = min(unique_headings) : 0.1 : max(unique_headings);
                     plot(neu_heading, Neuro_correct{i}, h{i},xi,cum_gaussfit(neu_perf{i},xi),f{i},'MarkerSize', 6, 'Linewidth', 1.5);
                     set(gca, 'TickDir', 'out', 'ylim',([0,1]), 'YTick', [0 0.5 1], 'FontSize', 16);
                     xlabel('Heading Angles');
                     box off
                     xlabel('Heading Angles');
                     ylabel('Rightward Choices');
                     title ('Neurometric function')
                 end
        end
        
end
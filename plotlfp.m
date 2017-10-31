function plotlfp(lfp,lfp_num,exp_name,plottype)

monk_id = lfp.monk_id;
session_id = lfp.session_id;
channel_no = lfp.channel_no;

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

switch exp_name
    case 'linearspeed'
        switch plottype
            case 'psth'
                % gather
                t = lfp.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:) = mean(lfp.(cond{i}).wave_pst);
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    plot(t,r(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num)]);
            case 'allpsth'
                % gather
                t = lfp.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:,:) = lfp.(cond{i}).wave_pst;
                end
                % plot
                figure; set(gcf,'Position',[200 100 300 700]);
                colorscale = 1/size(r,2):1/size(r,2):1;
                for i=1:size(r,1)
                    subplot(size(r,1),1,i); hold on;
                    for j=1:size(r,2)
                        plot(t',squeeze(r(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                    end
                    set(gca,'xlim',[-0.1 1.5],'ylim',[floor(min(r(:))*10)/10 ceil(max(r(:))*10)/10],...
                        'XTick',[0.5 1.2]-0.1,'XTickLabel',[],'TickDir','Out','Fontsize',16);
                end
                suptitle(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num)]);
        end
    case 'angularspeed'
        switch plottype
            case 'psth'
                % gather
                t = lfp.(cond{1}).time;
                for i=1:length(cond)
                    stim = lfp.(cond{i}).stim;
                    r_ccw(i,:) = mean(lfp.(cond{i}).wave_pst(stim<0,:));
                    r_cw(i,:) = mean(lfp.(cond{i}).wave_pst(stim>0,:));
                end
                % plot
                figure; hold on;
                for i=1:size(r_ccw,1)
                    plot(t,r_ccw(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num) '_{ccw}']);                figure; hold on;
                for i=1:size(r_cw,1)
                    plot(t,r_cw(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num) '_{cw}']);
            case 'psth_average'
                % gather
                t = lfp.(cond{1}).time;
                for i=1:length(cond)
                    stim = lfp.(cond{i}).stim;
                    clear r_ccw1 r_cw1;
                    for j=1:length(lfp_num)
                        r_ccw1(j,:) = mean(lfp_num(j).(cond{i}).wave_pst(stim<0,:));
                        r_cw1(j,:) = mean(lfp_num(j).(cond{i}).wave_pst(stim>0,:));
                    end
                    r_ccw(i,:) = squeeze(mean(r_ccw1));
                    r_cw(i,:) = squeeze(mean(r_cw1));
                end
                % plot
                figure; hold on;
                for i=1:length(cond)
                    plot(t,r_ccw(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) '_{ccw}']); figure; hold on;
                for i=1:length(cond)
                    plot(t,r_cw(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) '_{cw}']);
            case 'allpsth'
                % gather
                t = lfp.(cond{1}).time;  s = lfp.(cond{1}).stim;
                for i=1:length(cond)
                    r_ccw(i,:,:) = lfp.(cond{i}).wave_pst(s<0,:);    % ccw speeds
                    r_cw(i,:,:) = lfp.(cond{i}).wave_pst(s>0,:);     % cw speeds
                end
                % plot
                figure; set(gcf,'Position',[200 100 700 700]);
                colorscale = 1/size(r_cw,2):1/size(r_cw,2):1;
                for i=1:size(r_cw,1)
                    for j=1:size(r_cw,2)
                        subplot(size(r_cw,1),2,2*i-1); hold on;
                        plot(t',squeeze(r_cw(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                        set(gca,'xlim',[-0.1 1.5],'ylim',[floor(min(r_cw(:))*10)/10 ceil(max(r_cw(:))*10)/10],...
                            'XTick',[0.5 1.2]-0.1,'XTickLabel',[],'TickDir','Out','Fontsize',16);
                        subplot(size(r_ccw,1),2,2*i); hold on;
                        plot(t',squeeze(r_ccw(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                        set(gca,'xlim',[-0.1 1.5],'ylim',[floor(min(r_cw(:))*10)/10 ceil(max(r_cw(:))*10)/10],...
                            'XTick',[0.5 1.2]-0.1,'XTickLabel',[],'TickDir','Out','Fontsize',16);
                    end
                end
                suptitle(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num)]);
        end
         case '1DAzi'
        switch plottype
            case 'psth'
                % gather
                t = lfp.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:) = mean(lfp.(cond{i}).wave_pst);
                end
                % plot
                hold on;
                for i=1:size(r,1)
                    plot(t,r(i,:),'Color',[1 2 3]==i,'linewidth',2);
                end
                set(gca,'xlim',[-0.1 1.5],'XTick',[0 0.1 0.5 1.2 1.6]-0.1,...
                    'XTickLabel',[0 0.1 0.5 1.2 1.6],'TickDir','Out','Fontsize',16);
                title(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num)]);
            case 'allpsth'
                % gather
                t = lfp.(cond{1}).time;
                for i=1:length(cond)
                    r(i,:,:) = lfp.(cond{i}).wave_pst;
                end
                % plot
                figure; set(gcf,'Position',[200 100 300 700]);
                colorscale = 1/size(r,2):1/size(r,2):1;
                for i=1:size(r,1)
                    subplot(size(r,1),1,i); hold on;
                    for j=1:size(r,2)
                        plot(t',squeeze(r(i,j,:)),'Color',colorscale(j)*([1 2 3]==i),'linewidth',2);
                    end
                    set(gca,'xlim',[-0.1 1.5],'ylim',[floor(min(r(:))*10)/10 ceil(max(r(:))*10)/10],...
                        'XTick',[0.5 1.2]-0.1,'XTickLabel',[],'TickDir','Out','Fontsize',16);
                end
                suptitle(['m' num2str(monk_id) 's' num2str(session_id) 'ch' num2str(channel_no) '\_l' num2str(lfp_num)]);
        end
end
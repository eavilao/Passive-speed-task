%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for conflict heading discrimination task
%--	06/01/14 BL
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%save(['C:\Documents and Settings\Heide\My Documents\Matlab\adapt_hd_backdoor_files\' FILE '.mat']);

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);
temp_mask_status = data.moog_params(MASK_STATUS,:,MOOG);
temp_mask_radius = data.moog_params(MASK_RADIUS,:,MOOG);
temp_microstim = data.moog_params(MICROSTIM,:,MOOG);
temp_conflict_angle = data.moog_params(CONFLICT_ANGLE,:,MOOG);

trials = 1:length(temp_heading);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);
mask_status= temp_mask_status( select_trials );
mask_radius= temp_mask_radius( select_trials );
microstim = temp_microstim(select_trials );
conflict_angle = temp_conflict_angle(select_trials);

vest_heading = heading+conflict_angle/2; %%% obtain the real heading angle for the MOOG
heading = vest_heading;

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');
unique_heading_nonzero = unique_heading(unique_heading~=0);
unique_conflict_angle = munique(conflict_angle');


one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence)*length(unique_conflict_angle);

repetition = floor( length(heading)/one_repetition ); % take minimum repetition

% whether to plot performance over time
overtimeplot = 1;  % compute and plot
%overtimeplot = 0;  % not compute and plot

%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
LEFT = 1;
RIGHT = 2;
for i= 1 : length(total_trials) 
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
     %   choice(i) = RIGHT;
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end

correct_rate = [];
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_conflict_angle)    
         for i = 1:length(unique_heading)
             trials_select =logical( (heading == unique_heading(i)) & (conflict_angle == unique_conflict_angle(k)) & (motion_coherence==unique_motion_coherence(c)) ) ;
             rightward_trials = (trials_select & (choice == RIGHT) );
             rightward_rate = 1*sum(rightward_trials) / sum(trials_select);  
             fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);  
             fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
             fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select); 
         end 

         % the correct rate does not take coherence into account,temporarily 05/29/09
         trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_conflict_angle(k))  ) ;
         trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_conflict_angle(k))  ) ;
         trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_conflict_angle(k)) ); %exclude 0 headings
         correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);

         aa = find(fit_data_psycho_cum{c,k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
         fit_valid{c,k}(:,1) = fit_data_psycho_cum{c,k}(aa,1); 
         fit_valid{c,k}(:,2) = fit_data_psycho_cum{c,k}(aa,2);
         fit_valid{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
         
    end
end

% fit_valid = CrashorDistorted( fit_valid, FILE ); %%%% some data points may crash the pfit function. 6/1/14 BL

%%%%%% use Wichman's MLE method to estimate threshold and bias
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_conflict_angle)    
        wichman_psy = pfit(fit_valid{c,k}(1:end,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
        Thresh_psy{c,k} = wichman_psy.params.est(2);
        Bias_psy{c,k} = wichman_psy.params.est(1);
        psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];

%     %   similar way to fit data
%         [bb,tt] = cum_gaussfit_max1(fit_valid{c,k});
%         Thresh_psy{c,k} = tt;
%         Bias_psy{c,k} = bb;
%         psy_perf{c,k} =[bb,tt];
    end
end
% added by GY 12-04-07
% now this is the prediction when there are three stimuli conditions 
if length(unique_stim_type) ==3
    Thresh_pred = sqrt( Thresh_psy{1}^2*Thresh_psy{2}^2/(Thresh_psy{1}^2+Thresh_psy{2}^2) );
end
% this is the output, you can use it for plot of example cells
xi = min(unique_heading) : 0.1 : max(unique_heading);
for k = 1:length(unique_stim_type)  
    yi{k} = cum_gaussfit(psy_perf{k}, xi);
end
if length(unique_stim_type) ==3
    yi_pred = cum_gaussfit([Bias_psy{3},Thresh_pred], xi); % smoothed line for prediction with PSE at actual combined condition
end





% plot psychometric function here
symbo{1} = 'bo';    symbo{2} = 'ro';    symbo{3} = 'go';   symbo{4} = 'co';    symbo{5} = 'ko';    symbo{6} = 'c*';   symbo{7} = 'g*';    symbo{8} = 'r*';    symbo{9} = 'b*';
fitline{1} = 'b-'; fitline{2} = 'r-'; fitline{3} = 'g-'; fitline{4} = 'c-'; fitline{5} = 'k-'; fitline{6} = 'c--'; fitline{7} = 'g--'; fitline{8} = 'r--'; fitline{9} = 'b--';  


figure(2);
set(2,'Position', [200,100 700,800], 'Name', ['Heading Discrimination-Conflict (' FILE(1:9) ')']);
axes('position',[0.2,0.47, 0.6,0.4] );
% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];
xi = min(unique_heading) : 0.1 : max(unique_heading);
for c = 1:length(unique_motion_coherence) % different coherence level
        subplot(length(unique_motion_coherence), 2, 2*(c-1)+1);  
        p1 = plot(fit_valid{c,1}(:,1), fit_valid{c,1}(:,2), symbo{1},  xi, cum_gaussfit(psy_perf{c,1}, xi),  fitline{1}); hold on;
        p2 = plot(fit_valid{c,end}(:,1), fit_valid{c,end}(:,2), symbo{end},  xi, cum_gaussfit(psy_perf{c,end}, xi),  fitline{end});
        p3 = plot(fit_valid{c,2}(:,1), fit_valid{c,2}(:,2), symbo{2},  xi, cum_gaussfit(psy_perf{c,2}, xi),  fitline{2} );
        p4 = plot(fit_valid{c,end-1}(:,1), fit_valid{c,end-1}(:,2), symbo{end-1},  xi, cum_gaussfit(psy_perf{c,end-1}, xi),  fitline{end-1});
        p5 = plot(fit_valid{c,3}(:,1), fit_valid{c,3}(:,2), symbo{3},  xi, cum_gaussfit(psy_perf{c,3}, xi),  fitline{3} );
        p6 = plot(fit_valid{c,end-2}(:,1), fit_valid{c,end-2}(:,2), symbo{end-2},  xi, cum_gaussfit(psy_perf{c,end-2}, xi),  fitline{end-2} );
        p7 = plot(fit_valid{c,4}(:,1), fit_valid{c,4}(:,2), symbo{4},  xi, cum_gaussfit(psy_perf{c,4}, xi),  fitline{4} );
        p8 = plot(fit_valid{c,end-3}(:,1), fit_valid{c,end-3}(:,2), symbo{end-3},  xi, cum_gaussfit(psy_perf{c,end-3}, xi),  fitline{end-3} );
        p9 = plot(fit_valid{c,5}(:,1), fit_valid{c,5}(:,2), symbo{5},  xi, cum_gaussfit(psy_perf{c,5}, xi),  fitline{5} );
           
        xlabel('Heading Angles');   
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        xlim([-25 25]);
        title(['coherence = ', num2str(unique_motion_coherence(c))]); 
        if c==1
            legend([p1(2) p3(2) p5(2) p7(2) p9(2) p8(2) p6(2) p4(2) p2(2)], 'angle= -40deg', 'angle= -20deg', 'angle= -10deg', 'angle= -5deg', 'angle= 0deg', 'angle= 5deg', 'angle= 10deg', 'angle= 20deg', 'angle= 40deg');
        end
        hold off;
        
        subplot(length(unique_motion_coherence), 2, 2*(c-1)+2); hold on; 
        bias_conflict = [];
        for kk = 1:length(unique_conflict_angle)
            bias_conflict(kk) = psy_perf{c, kk}(1);
        end
        plot(unique_conflict_angle(1:length(unique_conflict_angle)), bias_conflict, 'ks-');
        xlabel('Disparity (deg)');
        ylabel('Bias (deg)');
        title(['coherence = ', num2str(unique_motion_coherence(c))]);
        set(gca,'XTick',[-40 -20 0  20 40]);
        xlim([-50 50]);
        hold off;
        
end

% filename = ['c:\data_analysis\conflict_heading\', FILE(1:9), '.mat'];      %%   for Jasper, Quigley, Gunnison 
% filename = ['c:\data_analysis\conflict_heading\', FILE(1:10), '.mat'];     %%   for Yosemite  
% save(filename, 'fit_valid', 'unique_heading', 'unique_conflict_angle', 'unique_motion_coherence');


%---------------------------------------------------------------------------------------
return;
%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for slant discrimination task
%--	11/8/11 Ari
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
disp('Slant Discrimination')
TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%save(['C:\Documents and Settings\Heide\My Documents\Matlab\adapt_hd_backdoor_files\' FILE '.mat']);
%get the column of values for azimuth and elevation and stim_type
% % temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
% % temp_elevation = data.moog_params(ELEVATION,:,MOOG);
% % temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% % temp_heading   = data.moog_params(HEADING, :, MOOG); 
% % temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
% % temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
% % temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
% % temp_total_trials = data.misc_params(OUTCOME, :);
% % temp_mask_status = data.moog_params(MASK_STATUS,:,MOOG);
% % temp_mask_radius = data.moog_params(MASK_RADIUS,:,MOOG);
% % temp_microstim = data.moog_params(MICROSTIM,:,MOOG);

temp_total_trials = data.misc_params(OUTCOME, :);
temp_heading   = data.dots_params(DOTS_HGRAD_MAG, :, 1 );

trials = 1:length(temp_heading);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% % stim_type = temp_stim_type( select_trials );
% % heading = temp_heading( select_trials );
% % amplitude= temp_amplitude( select_trials );
% % num_sigmas= temp_num_sigmas( select_trials );
% % motion_coherence = temp_motion_coherence(select_trials);
% % total_trials = temp_total_trials( select_trials);
% % mask_status= temp_mask_status( select_trials );
% % mask_radius= temp_mask_radius( select_trials );
% % microstim = temp_microstim(select_trials );

total_trials = temp_total_trials( select_trials);
heading = temp_heading( select_trials );

% unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
% unique_amplitude = munique(amplitude');
% unique_num_sigmas = munique(num_sigmas');
% unique_motion_coherence = munique(motion_coherence');
% unique_mask_status = munique(mask_status');
% unique_mask_radius = munique(mask_radius');
% unique_microstim = munique(microstim');
% unique_heading_nonzero = unique_heading(unique_heading~=0);

if 1%%length(unique_motion_coherence)==1
   one_repetition = length(unique_heading);%*length(unique_stim_type);
else
   one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence)-length(unique_heading); 
end
repetition = floor( length(heading)/one_repetition ); % take minimum repetition

% whether to plot performance over time
% overtimeplot = 1;  % compute and plot
overtimeplot = 0;  % not compute and plot

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
for c = 1:1%length(unique_motion_coherence) % different coherence level
    for k = 1:1%length(unique_stim_type)    
         for i = 1:length(unique_heading)
%              if unique_stim_type(k) == 1 % for vestibular condition, take all the data regardless of visual coherence
%                  trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) ) ;
%              else 
                 trials_select =logical( (heading == unique_heading(i)));% & (stim_type==unique_stim_type(k)) & (motion_coherence==unique_motion_coherence(c)) ) ;
%              end
             rightward_trials = (trials_select & (choice == RIGHT) );
             rightward_rate = 1*sum(rightward_trials) / sum(trials_select);  
             fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);  
             fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
             fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select); 
         end 
%          halfheading = length(unique_heading_nonzero/2);
%          for j = 1: halfheading
%              trials_left = find( (heading==unique_heading_nonzero(halfheading+1-j)) & (choice==LEFT) & (stim_type==unique_stim_type(k))  ) ;
%              trials_right  = find( (heading==unique_heading_nonzero(halfheading+j)) & (choice==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
%              trials_all = find( ((heading==unique_heading_nonzero(halfheading+1-j)|(heading==unique_heading_nonzero(halfheading+j)) & (stim_type==unique_stim_type(k)) ); 
%              correct_rate(k,j) = (length(trials_right)+length(trials_left))/length(trials_all);
%              % for later weibull fit
%              fit_valid_weibull{c,k}(j,1) = unique_heading_nonzero(halfheading+j);
%              fit_valid_weibull{c,k}(j,2) = correct_rate(k,j);
%              fit_valid_weibull{c,k}(j,3) = fit_data_psycho_cum{c,k}(aa,3);
%          end
         % the correct rate does not take coherence into account,temporarily 05/29/09
         trials_rightward = find( (heading > 0) & (choice==RIGHT));% & (stim_type==unique_stim_type(k))  ) ;
         trials_leftward  = find( (heading < 0) & (choice==LEFT));% & (stim_type==unique_stim_type(k))  ) ;
         trials_all = find( ((heading < 0)|(heading > 0)));% & (stim_type==unique_stim_type(k)) ); %exclude 0 headings
         correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);

         aa = find(fit_data_psycho_cum{c,k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
         fit_valid{c,k}(:,1) = fit_data_psycho_cum{c,k}(aa,1); 
         fit_valid{c,k}(:,2) = fit_data_psycho_cum{c,k}(aa,2);
         fit_valid{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
         
%          % for later weibull fit use
%          fit_valid_weibull{c,k}(:,1) = unique_heading( unique_heading>0) ); 
%          fit_valid_weibull{c,k}(:,2) = correct_rate(k,:);
%          fit_valid_weibull{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
    end
end
%%%%%% use Wichman's MLE method to estimate threshold and bias
 
figure(2); title(FILE)
wichman_psy = pfit(fit_valid{1,1},'shape','cumulative gaussian','runs',2000,'n_intervals',1,'cuts',[0.159 0.5 0.841],'plot_opt','plot without stats','sens',0, 'compute_stats',1,'verbose',0);
Thresh_psy = wichman_psy.params.est(2)
Bias_psy = wichman_psy.params.est(1)
psy_perf = [wichman_psy.params.est(1),wichman_psy.params.est(2)]

text(-0.03,1,strcat('Threshold = ','',num2str(Thresh_psy)))
text(-0.03,0.9,strcat('Bias = ','',num2str(Bias_psy)))
%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for heading discrimination Gabor task
%--	Jing 10/19/2016
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_Gabor(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; 

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
temp_phase = data.moog_params(PHASE,:,MOOG);

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
phase = temp_phase(select_trials );

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');
unique_heading_nonzero = unique_heading(unique_heading~=0);
unique_phase = munique(phase');

one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_phase);
repetition = floor( length(heading)/one_repetition ); % take minimum repetition

% whether to plot performance over time
overtimeplot = 1;  % compute and plot

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
for c = 1:length(unique_phase) % different FP location
    for k = 1:length(unique_stim_type)    
         for i = 1:length(unique_heading)
            trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (phase==unique_phase(c)) ) ;
            rightward_trials = (trials_select & (choice == RIGHT) );
            rightward_rate = 1*sum(rightward_trials) / sum(trials_select);  
            fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);  
            fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
            fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select); 
         end 

         trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_stim_type(k)) & (phase==unique_phase(c)) ) ;
         trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_stim_type(k)) & (phase==unique_phase(c)) ) ;
         trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k)) & (phase==unique_phase(c))); %exclude 0 headings
         correct_proportion{c,k} = (length(trials_rightward)+length(trials_leftward))/length(trials_all);

         aa = find(fit_data_psycho_cum{c,k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
         fit_valid{c,k}(:,1) = fit_data_psycho_cum{c,k}(aa,1); 
         fit_valid{c,k}(:,2) = fit_data_psycho_cum{c,k}(aa,2);
         fit_valid{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
    end
end
%%%%%% use Wichman's MLE method to estimate threshold and bias
for c = 1:length(unique_phase) % different FP location
    for k = 1:length(unique_stim_type)    
        wichman_psy = pfit(fit_valid{c,k}(1:end,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');  
        Thresh_psy{c,k} = wichman_psy.params.est(2);
        Bias_psy{c,k} = wichman_psy.params.est(1);
        psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
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

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot psychometric, neurometric, CP over time
% % run the slide threshold over time, see whether performance fluctuate across time
% not work for coherence, temporarily 05/29/09
if overtimeplot == 1
    span = 5;  % calculate threshod every ? repeats;
    slide = 1;  % slide threshod with increment of ? repeats;
    BegTrial_shift = BegTrial;
    EndTrial_shift = BegTrial_shift + span*one_repetition-1;
    n=0;
    while EndTrial_shift <= EndTrial    
        n = n + 1;
        select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
        stim_type_shift = temp_stim_type( select_trials_shift );
        heading_shift = temp_heading( select_trials_shift );
        phase_shift = temp_phase(select_trials_shift );
        unique_stim_type_shift = munique(stim_type_shift');
        unique_heading_shift = munique(heading_shift');
        unique_phase_shift = munique(phase_shift');
        total_trials_shift = temp_total_trials( select_trials_shift);       
        condition_shift = stim_type_shift;
        unique_condition_shift = unique_stim_type_shift;
       
        for c = 1:length(unique_phase_shift)
            for k = 1:length(unique_condition_shift)
                for i = 1:length(unique_heading)
                    trials_shift =logical( (heading_shift == unique_heading(i)) & (condition_shift == unique_condition_shift(k)) & (phase_shift == unique_phase_shift(c))) ;
                    correct_trials_shift = (trials_shift & (total_trials_shift == CORRECT) );
                    % make 'S' curve by using the rightward choice for y-axis
                    if sum(trials_shift)>0
                        if ( unique_heading(i) < 0 )
                            correct_rate_shift(i) = 1 - 1*sum(correct_trials_shift) / sum(trials_shift);
                        else
                            correct_rate_shift(i) = 1*sum(correct_trials_shift) / sum(trials_shift);
                        end
                    end
                    Trials_num(i) = sum(trials_shift);
                end
                aa = find(correct_rate_shift >-1 );
                for j = 1:length(aa)
                    fit_data_psycho_cum_shift{c,k}(j, 1) = fit_data_psycho_cum{c,k}(aa(j), 1);
                    fit_data_psycho_cum_shift{c,k}(j, 2) = correct_rate_shift(aa(j));
                    fit_data_psycho_cum_shift{c,k}(j, 3) = Trials_num(aa(j));
                end
                % this fixes a strange error: cum_gaussfit/pfit sometimes fail when pct choices are all 0's or 1's -CRF 8-13-08
                if fit_data_psycho_cum_shift{c,k}(:,2)==0 | fit_data_psycho_cum_shift{c,k}(:,2)==1
                    fit_data_psycho_cum_shift{c,k}(fit_data_psycho_cum_shift{c,k}==0) = 0.001;
                    fit_data_psycho_cum_shift{c,k}(fit_data_psycho_cum_shift{c,k}==1) = 0.999;
                end
                [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{c,k}); % to save time, use a different fit method
                psy_thresh_shift{c}(k,n) = tt;
                psy_bias_shift{c}(k,n) = bb;  % added Bias, CRF 11-5-09
                if n == 11
                    xxx = 0;
                end
            end
        end        
        BegTrial_shift = BegTrial_shift + slide*one_repetition;
        EndTrial_shift = BegTrial_shift + span*one_repetition-1;
   end
end

% plot psychometric function here
symbo{1,1} = 'bo';    symbo{2,1} = 'ro';    symbo{3,1} = 'go'; 
symbo{4,1} = 'ko';    symbo{5,1} = 'mo';    symbo{6,1} = 'co'; 

fitline{1,1} = 'b-';    fitline{2,1} = 'r-';    fitline{3,1} = 'g-'; 
fitline{4,1} = 'k-';    fitline{5,1} = 'm-';    fitline{6,1} = 'c-'; 

figure(2);
set(2,'Position',[200,50 1400,950], 'Name', 'Heading Discrimination-Gabor');
axes('position',[0.06,0.1, 0.4,0.65] );
% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];
xi = min(unique_heading) : 0.1 : max(unique_heading);
for c = 1:length(unique_phase) 
    for k = 1:length(unique_stim_type)      
        plot(unique_heading, fit_valid{c,k}(:,2), symbo{c,k},  xi, cum_gaussfit(psy_perf{c,k}, xi),  fitline{c,k} );
        xlabel('Heading Angles');   
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        hold on;  
    end
    legend_txt{c*2-1} = [num2str(unique_phase(c))];
    legend_txt{c*2} = [''];
end
y=0:0.1:1.0;
x = zeros(size(y));
plot(x,y,'--k');
legend(legend_txt, 'Location','northwest');

% output some text of basic parameters in the figure
axes('position',[0.1,0.65, 0.4,0.3] );
xlim( [0,50] );
ylim( [2,10] );
text(0, 10, FILE);
text(15,10,'phase =');
text(30,10,'repeats =');
text(25,10,num2str(unique_phase) );
text(40,10,num2str(repetition) );
text(10,8.5, '   u                           sigma                   correct rate');

for c = 1:length(unique_phase)   %different FP location
    for k = 1:length(unique_stim_type)
        text(0,9-k-(c-1)*0.5, num2str(unique_stim_type(k))); 
        text(10,9-k-(c-1)*0.5,num2str(Bias_psy{c,k}) );
        text(20,9-k-(c-1)*0.5,num2str(Thresh_psy{c,k}) );
        text(30,9-k-(c-1)*0.5,num2str(correct_proportion{c,k}) );
    end
end
axis off;

% % plot psycho over time
if overtimeplot ==1
    axes('position',[0.55,0.5, 0.4,0.25] );
    for c = 1:length(unique_phase)   
        for k = 1:length(unique_stim_type)
            m_threshold = psy_thresh_shift{c}(k,:);
            index = find(m_threshold > 150);
            m_threshold(index) = 150;
            plot(3:n+2,m_threshold, fitline{c,k});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold on;
            xlabel('Repetition');
            ylabel('Threshold');
            xlim([0+3, n+2]);
        end
    end;
        
    axes('position',[0.55,0.1, 0.4,0.25] );
    for c = 1:length(unique_phase)
        for k = 1:length(unique_stim_type)
            m_bias = psy_bias_shift{c}(k,:);
            index = find(m_bias > 100);
            m_bias(index) = 100;
            index = find(m_bias < -100);
            m_bias(index) = -100;
            plot(3:n+2,m_bias, fitline{c,k});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold on;
            xlabel('Repetition');
            ylabel('Bias');
            xlim([0+3, n+2]);
        end
    end
end
orient tall;

return;
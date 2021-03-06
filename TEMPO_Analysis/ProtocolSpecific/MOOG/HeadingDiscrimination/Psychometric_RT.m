%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for  Reaction Time Task heading discrimination task
%--	01/24/09 Tunde
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_RT(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

using_batch = 0;
% this variable used to remove the option of multiple bins and setting the
% bin number to one when running a batch

if using_batch == 1
    cd Z:\Users\Adhira
end


%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);


% this part parses out 
% temp_event_data = data.event_data(1, :, :);
% temp_event_data = reshape(temp_event_data, size(temp_event_data, 2), size(temp_event_data, 3));


    %for each selected trial, find all of the non-zero elements of event_data; these are the code values
    %then sprintf them to a string and display them in the console window.
    num_trials = size(data.event_data, 3);
    stim_durations = zeros(1, num_trials);
    pre_stim_delay = zeros(1, num_trials);
    trials = 1:length(temp_heading);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
%     stim_durations = stim_durations(select_trials); %Edit AS - for partial set of trials

    
    for i = 1:num_trials
%         stim_start = find(data.event_data(1,:,i) == 4 );
%         stim_end = find(data.event_data(1,:,i) == 5 );
        stim_start = find(data.event_data(1,:,i) == VSTIM_ON_CD );
        stim_end = find(data.event_data(1,:,i) == VSTIM_OFF_CD );
        monk_in_start_window = find(data.event_data(1,:,i) == IN_FIX_WIN_CD );
%         stim_end
%         stim_start
%         monk_in_start_window
        if (isempty(stim_start|stim_end|monk_in_start_window) == 1)%Some trials have this error - AS
            stim_durations(i) = 0; 
            pre_stim_delay(i) = 0;
        else
            stim_durations(i) = (stim_end - stim_start);
            pre_stim_delay(i) = (stim_start - monk_in_start_window);
        end
    end
%     stim_durations = stim_durations(select_trials);
   


stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);


unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');


one_repetition = length(unique_heading)*length(unique_stim_type);
repetition = floor( length(heading)/one_repetition ); % take minimum repetition

% whether to plot performance over time
%overtimeplot = 1;  % compute and plot
overtimeplot = 1;  % not compute and plot

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


% %specifying number of bins to use


% hist 

max_stim_duration = max(stim_durations);
min_stim_duration = min(stim_durations);
% bin_width = (max_stim_duration - min_stim_duration)/num_bins;
% b1s = min_stim_duration;
% b2s = (min_stim_duration + bin_width * 1);
% b3s = (min_stim_duration + bin_width * 2);
% b3e = max_stim_duration;
% bin_sections = [b1s, b2s, b3s, b3e];
% bin_sections_mid = [(b1s + b2s)/2, (b2s + b3s)/2, (b3s + b3e)/2];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % max_stim_duration = max(stim_durations);
% % % % % min_stim_duration = min(stim_durations);
% % % % % relative_bin_range = max_stim_duration - min_stim_duration;
% % % % % 
% % % % % relative_bin_range = relative_bin_range;
% % % % % 
% % % % % t4_relative_bin_range = 1.5^(relative_bin_range); %transforming the relative bin width to a linear scale for dividing into equally sized bins
% % % % % 
% % % % % bin_width_t4_relative_bin_range = t4_relative_bin_range/num_bins; %equally sized linear bin
% % % % % 
% % % % % t4_bin_range = 0:bin_width_t4_relative_bin_range:t4_relative_bin_range; %equally sized linear bin
% % % % % 
% % % % % 
% % % % % bin_sections = logb(t4_bin_range, 1.5); %log sized bins 
% % % % % bin_sections = bin_sections; %removing the infinity point due to log of zero



if using_batch == 0
    %     num_bins = input('Please input your number of bins: ', 's');
    %     num_bins = str2num(num_bins);
    num_bins = 5;
elseif using_batch == 1
    num_bins = 1; %number of bins to use
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % old way of calculating log bins that's incorrect due to non linearity of the log function 
% stimulus_dur_range = max_stim_duration - min_stim_duration;
% 
% log_base = exp(1); %logarithm base
% bins_vec = 1:(num_bins + 1); %vector numbering for bins
% log_spaced_bins = diff(logb(bins_vec, log_base)); %log spaced bins 
% sum_log_spaced_bins = sum(log_spaced_bins); % summing all log spaced bins to find the bin base number
% bin_base = stimulus_dur_range/sum_log_spaced_bins;
% log_spaced_bins = bin_base * log_spaced_bins;
% sum(log_spaced_bins);
% 
% bin_sections = zeros(size(bins_vec));
% bin_sections(1) = min_stim_duration;
% bin_sections(end) = max_stim_duration;

% % filling in the remaining bin sections in the middle
% for i = 1:(length(bin_sections) - 1)
%     bin_sections(i + 1) = bin_sections(i) + log_spaced_bins(i);
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

bin_type = 1; 
% bin_type = 0 -- Log based bins
% bin_type = 1 -- equal trials per bin
if bin_type == 0
    % % % % % % % % % % using log based bins % % % % % % % % % %
    log_max_stim_dur = log(max_stim_duration);
    log_min_stim_dur = log(min_stim_duration);


    log_stim_dur_range = log_max_stim_dur - log_min_stim_dur;
    log_stim_binsize = log_stim_dur_range/num_bins; %size of one bin in the log space

    log_bin_sections = zeros(1, num_bins + 1);
    log_bin_sections(1) = log_min_stim_dur;
    log_bin_sections(end) = log_max_stim_dur;

    % filling in the remaining bin sections in the middle
    for i = 1:(length(log_bin_sections) - 1)
        log_bin_sections(i + 1) =  log_bin_sections(i) + log_stim_binsize;
    end

    % testa = exp(log_bin_sections)
    %
    % diff(testa)
    bin_sections = exp(log_bin_sections);

    %calculationg bin mid points for the logarighmic bin plot
    bin_sections_mid = []; %middle of the bins
    for i = 2:length(bin_sections)
        %     i
        bin_sections_mid(i - 1) = ( bin_sections(i) + bin_sections(i - 1) )/2;
    end

    % bin_width = (max_stim_duration - min_stim_duration)/num_bins;
    % b1s = min_stim_duration;
    % b2s = (min_stim_duration + bin_width * 1);
    % b3s = (min_stim_duration + bin_width * 2);
    % b3e = max_stim_duration;
    % bin_sections = [b1s, b2s, b3s, b3e];
    % bin_sections_mid = [(b1s + b2s)/2, (b2s + b3s)/2, (b3s + b3e)/2];

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
elseif bin_type == 1
    bin_sections_prctile = zeros(1, num_bins + 1);
    for i = 1:(num_bins + 1)
        bin_sections_prctile(i) = (i - 1) * 100/num_bins;       
    end
    
    
    bin_sections = prctile(stim_durations, bin_sections_prctile);
    bin_sections_mid = []; %middle of the bins
    for i = 2:length(bin_sections)
        %     i
        bin_sections_mid(i - 1) = ( bin_sections(i) + bin_sections(i - 1) )/2;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % working on the single psychometric plots% % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %





correct_rate_single = [];
for k = 1:length(unique_stim_type)    
    for i = 1:length(unique_heading)
            for m = 1:length(unique_motion_coherence)
                trials_select_single = logical( (heading == unique_heading(i)) &  ...
                    (stim_type == unique_stim_type(k)) & ...
                    (motion_coherence == unique_motion_coherence(m)) );
                correct_trials_single = (trials_select_single & (choice == RIGHT) );
                correct_rate_single(k, i, m) = 1*sum(correct_trials_single) / sum(trials_select_single);
                if m == 2 && k == 1 % this takes care of the vestibular condition with the other coherence.
                    %i.e. you can only have one vestibular condition  when
                    %there are 2 coherences, no repeats.
                    fit_data_psycho_cum_single{k, m}(i, 1) = fit_data_psycho_cum_single{k, 1}(i, 1);
                    fit_data_psycho_cum_single{k, m}(i, 2) = fit_data_psycho_cum_single{k, 1}(i, 2);
                    fit_data_psycho_cum_single{k, m}(i, 3) = fit_data_psycho_cum_single{k, 1}(i, 3);
                else
                    fit_data_psycho_cum_single{k, m}(i, 1) = unique_heading( i );
                    fit_data_psycho_cum_single{k, m}(i, 2) = correct_rate_single(k, i, m);
                    fit_data_psycho_cum_single{k, m}(i, 3) = sum(trials_select_single);
                end
            end
    end
end








% 
% correct_rate_single = [];
% for k = 1:length(unique_stim_type)    
%      for i = 1:length(unique_heading)
%              trials_select_single =logical( (heading == unique_heading(i)) ) ;
%              correct_trials_single = (trials_select_single & (choice == RIGHT) );
%              correct_rate_single(k, i) = 1*sum(correct_trials_single) / sum(trials_select_single);
%              fit_data_psycho_cum_single{k}(i, 1) = unique_heading( i );
%              fit_data_psycho_cum_single{k}(i, 2) = correct_rate_single(k, i);
%              fit_data_psycho_cum_single{k}(i, 3) = sum(trials_select_single);
%      end         
% end


rt_vest_single = cell(size(unique_motion_coherence));
rt_vis_single = cell(size(unique_motion_coherence));
rt_comb_single = cell(size(unique_motion_coherence));
rt_all_single = cell(size(unique_motion_coherence));



if length(unique_stim_type) == 3 %i.e. all 3 conditions/modalities were run on the monkey
    for k = 1:length(unique_stim_type)
        for m = 1:length(unique_motion_coherence)
            trials_select = logical( (stim_type == unique_stim_type(k)) &  ...
                motion_coherence == unique_motion_coherence(m) );
            if unique_stim_type(k) == 1
                if m == 1
                    rt_vest_single{m} = nonzeros(trials_select .* stim_durations);
                elseif m == 2 
                    rt_vest_single{m} = rt_vest_single{1}; %repeating this to keep the code from crashing as vestibular runs will not have 2 different coherences
                end
            elseif unique_stim_type(k) == 2
                rt_vis_single{m} = nonzeros(trials_select .* stim_durations);
            elseif unique_stim_type(k) == 3
                rt_comb_single{m} = nonzeros(trials_select .* stim_durations);
            end
        end
    end
elseif length(unique_stim_type) == 1 % i.e. only one condition/modality was run
    trials_select = logical( (stim_type==unique_stim_type) );
    if unique_stim_type(k) == 1
        rt_vest_single = nonzeros(trials_select .* stim_durations);
    elseif unique_stim_type(k) == 2
        rt_vis_single = nonzeros(trials_select .* stim_durations);
    elseif unique_stim_type(k) == 3
        rt_comb_single = nonzeros(trials_select .* stim_durations);
    end
elseif length(unique_stim_type) == 2 % i.e. only 2 condtions were run on the monkey
    for k = 1:length(unique_stim_type)
        trials_select = logical( (stim_type==unique_stim_type(k)) );
        if unique_stim_type(k) == 1
            rt_vest_single = nonzeros(trials_select .* stim_durations);
        elseif unique_stim_type(k) == 2
            rt_vis_single = nonzeros(trials_select .* stim_durations);
        elseif unique_stim_type(k) == 3
            rt_comb_single = nonzeros(trials_select .* stim_durations);
        end
    end
end






%%%%%% use Wichman's MLE method to estimate threshold and bias
for k = 1:length(unique_stim_type)
    for m = 1:length(unique_motion_coherence)
        
        wichman_psy_single = pfit(fit_data_psycho_cum_single{k, m},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');
        Thresh_psy_single{k, m} = wichman_psy_single.params.est(2);
        Bias_psy_single{k, m} = wichman_psy_single.params.est(1);
        psy_perf_single{k, m} = [wichman_psy_single.params.est(1),wichman_psy_single.params.est(2)];

        % %   similar way to fit data
        %     [bb,tt] = cum_gaussfit_max1(fit_valid{k});
        %     Thresh_psy{k} = tt;
        %     Bias_psy{k} = bb;
        %     psy_perf{k} =[bb,tt];
    end
end


symbo_NaN{1,1} = 'bx';
symbo_NaN{2,1} = 'rx';
symbo_NaN{3,1} = 'gx';

symbo{1,1} = 'bo';
symbo{2,1} = 'ro';
symbo{3,1} = 'go';

fitline{1,1} = 'b-';
fitline{2,1} = 'r-';
fitline{3,1} = 'g-';

for i = 2:num_bins
    symbo{1,i} = 'bo';
    symbo{2,i} = 'ro';
    symbo{3,i} = 'go';

    fitline{1,i} = 'b-';
    fitline{2,i} = 'r-';
    fitline{3,i} = 'g-';
end




fig_list = cell(size(unique_motion_coherence));
fig_handle_list = cell(size(unique_motion_coherence));
fig6_handle = figure(6);
% if length(unique_motion_coherence) == 1 
%     fig2_handle = figure(2);
%     fig_handle_list{6} = 'fig6_handle';
%     fig_list{1} = 'figure(2)';
% elseif length(unique_motion_coherence) == 2
%     fig2_handle = figure(2);
%     fig_list{1} = 'figure(2)';
%     fig_handle_list{1} = 'fig2_handle';
%            
%     fig4_handle = figure(4);
%     fig_list{2} = 'figure(4)';
%     fig_handle_list{2} = 'fig4_handle';
% end
txt_nums_a = [0.2, 0.6];
set(fig6_handle,'Position', [1 54 1280 1024-150], 'Name', 'Psychometric plots');
xi = min(unique_heading) : 0.1 : max(unique_heading);
axes('position',[0.01,0.95, 0.99,0.05] );
for m = 1:length(unique_motion_coherence)
    text(txt_nums_a(m), 0.2, ['Psychometric Plots ' FILE  '; Coherence: ' num2str(unique_motion_coherence(m)) '%'])
end
axis off


for k = 1:length(unique_stim_type)
        for n = 1:length(unique_motion_coherence)
%             eval([fig_list{n}])
            figure(6)
            subplot(1, 2, n)
            plot(unique_heading, reshape(correct_rate_single(k, :, n), 1, length(unique_heading)), symbo{unique_stim_type(k)}, xi, cum_gaussfit(psy_perf_single{k, n}, xi),  fitline{unique_stim_type(k)} );
            hold on
%             if ~isempty(row_NaNs_all{k}) % if NaNs were replaced then plot X's
%                 for m = 1:size(row_NaNs_all{k}, 1)
%                     plot( fit_data_psycho_cum{k, n}(row_NaNs_all{k}(m), 1), fit_data_psycho_cum{k, n}(row_NaNs_all{k}(m), column_NaNs_all{k}(m)), symbo_NaN{unique_stim_type(k)} )
%                 end
%             end
            xlabel('Heading Angles');
            ylim([0,1]);
            ylabel('Rightward Choices');
            %         title(['Stimulus duration: ' num2str(bin_sections_mid(j)) ' ms']) % middle of bins
%             title([ 'Stimulus duration bin range: ' num2str(bin_sections(j)) 'ms to ' num2str(bin_sections(j + 1)) 'ms' ])
            set(gca, 'YTickMode','auto');
            set(gca, 'xTickMode','auto');
            hold on;
            legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
            legend_txt{k*2} = [''];
        end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    unique_heading_max = max(unique_heading);
    unique_heading_min = min(unique_heading);
    diff_max_min_heading = unique_heading_max - unique_heading_min;
    mid_max_min_heading = diff_max_min_heading/2 + unique_heading_min;
    mid_plus_max_min_heading = (diff_max_min_heading*3/4) + unique_heading_min;


    for k = 1:length(unique_stim_type)
            for m = 1:length(unique_motion_coherence)
%                 eval([fig_list{m}])
                figure(6)
                if unique_stim_type(k) == 1 %vestibular
                    subplot(1, 2, m)
                    %             text(0,8-(k*0 + j), num2str(unique_stim_type(k)));
                    text(mid_max_min_heading, 0.5, [{'\color{black}mu' }]  );
                    text(mid_plus_max_min_heading, 0.5, [{'\color{black}Sigma'  }]  );
                    text(mid_max_min_heading, 0.45, [{'\color{blue}' num2str(Bias_psy_single{k, m}) }]  );
                    text(mid_plus_max_min_heading, 0.45, [{'\color{blue}' num2str(Thresh_psy_single{k, m}) }]  );
                elseif unique_stim_type(k) == 2 %visual
                    %             text(0,8-(k*1.5 + j), num2str(unique_stim_type(k)));
                    subplot(1, 2, m)
                    text(mid_max_min_heading, 0.5, [{'\color{black}mu' }]  );
                    text(mid_plus_max_min_heading, 0.5, [{'\color{black}Sigma'  }]  );
                    text(mid_max_min_heading, 0.4, [{'\color{red}' num2str(Bias_psy_single{k, m}) }]  );
                    text(mid_plus_max_min_heading, 0.4, [{'\color{red}' num2str(Thresh_psy_single{k, m}) }]  );
                elseif unique_stim_type(k) == 3 %combined
                    %             text(0,8-(k*2 + j), num2str(unique_stim_type(k)));
                    subplot(1, 2, m)
                    text(mid_max_min_heading, 0.5, [{'\color{black}mu' }]  );
                    text(mid_plus_max_min_heading, 0.5, [{'\color{black}Sigma'  }]  );
                    text(mid_max_min_heading, 0.35, [{'\color{green}' num2str(Bias_psy_single{k, m}) }]  );
                    text(mid_plus_max_min_heading, 0.35, [{'\color{green}' num2str(Thresh_psy_single{k, m}) }]  );
                    %             text(mid_max_min_heading, 0.35, [{'\color{black}' num2str(Bias_psy{k, j}) }]  );
                    %             text(mid_plus_max_min_heading, 0.35, [{'\color{black}' num2str(Thresh_psy{k, j}) }]  );
                end
            end
    end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
correct_rate = [];
for k = 1:length(unique_stim_type)    
    for i = 1:length(unique_heading)
        for j = 1:num_bins
            for m = 1:length(unique_motion_coherence)
                trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & ...
                    stim_durations >= bin_sections(j) & stim_durations <= bin_sections(j + 1) & ...
                    motion_coherence == unique_motion_coherence(m) ) ;
                correct_trials = (trials_select & (choice == RIGHT) );
                correct_rate(k, j, i, m) = 1*sum(correct_trials) / sum(trials_select);
                if m == 2 && k == 1 % this takes care of the vestibular condition with the other coherence.
                    %i.e. you can only have one vestibular condition  when
                    %there are 2 coherences, no repeats.
                    fit_data_psycho_cum{k, j, m}(i, 1) = fit_data_psycho_cum{k, j, 1}(i, 1);
                    fit_data_psycho_cum{k, j, m}(i, 2) = fit_data_psycho_cum{k, j, 1}(i, 2);
                    fit_data_psycho_cum{k, j, m}(i, 3) = fit_data_psycho_cum{k, j, 1}(i, 3);
                else
                    fit_data_psycho_cum{k, j, m}(i, 1) = unique_heading( i );
                    fit_data_psycho_cum{k, j, m}(i, 2) = correct_rate(k, j, i, m);
                    fit_data_psycho_cum{k, j, m}(i, 3) = sum(trials_select);
                end
            end
        end
    end
% % % commented out at the behest of Yong --- Tunde
%      trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
%      trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_stim_type(k))  ) ;
%      trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k)) ); %exclude 0 headings
%      correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);
     
%      aa = find(fit_data_psycho_cum{k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
%      fit_valid{k}(:,1) = fit_data_psycho_cum{k}(aa,1); 
%      fit_valid{k}(:,2) = fit_data_psycho_cum{k}(aa,2);
%      fit_valid{k}(:,3) = fit_data_psycho_cum{k}(aa,3);
end


ves_fit_data_psycho_cum = zeros(length(unique_heading), 3, num_bins, length(unique_motion_coherence));
vis_fit_data_psycho_cum = zeros(length(unique_heading), 3, num_bins, length(unique_motion_coherence));
comb_fit_data_psycho_cum = zeros(length(unique_heading), 3, num_bins, length(unique_motion_coherence));

for i = 1:num_bins
    for j = 1:length(unique_motion_coherence)
        for k = 1:length(unique_stim_type)
            if length(unique_stim_type) == 3 %i.e. all 3 conditions were run on the monkey
                ves_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                vis_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{2, i, j};
                comb_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{3, i, j};
            elseif length(unique_stim_type) == 1 % only one condition was run
                if unique_stim_type == 1 %vestibular
                    ves_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                elseif unique_stim_type == 2 %visual
                    vis_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                elseif unique_stim_type == 3 %combined
                    comb_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                end
            elseif length(unique_stim_type) == 2 % any 2 conditions were run ---- Added by Adhira 01/05/11
                if unique_stim_type(k) == 1 %vestibular
                    ves_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                elseif unique_stim_type(k) == 2 %visual
                    vis_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                elseif unique_stim_type(k) == 3 %combined
                    comb_fit_data_psycho_cum(:, :, i, j) = fit_data_psycho_cum{1, i, j};
                end
            end
        end
    end
end

computername = getenv('COMPUTERNAME'); 
if strcmp(computername, 'SAC') == 1
disp('output for debugging (Tunde) ....')
     total_trials = zeros(3, 1);
     total_trials(1) = sum(sum(ves_fit_data_psycho_cum(:, 3, :)));
     total_trials(2) = sum(sum(vis_fit_data_psycho_cum(:, 3, :)));
     total_trials(3) = sum(sum(comb_fit_data_psycho_cum(:, 3, :)));
     
     total_trials
     
%       total_trials_per_sec = zeros(num_bins, 1)
      ves_bins_total = sum(ves_fit_data_psycho_cum(:, 3, :)); 
      vis_bins_total = sum(vis_fit_data_psycho_cum(:, 3, :));
      comb_bins_total = sum(comb_fit_data_psycho_cum(:, 3, :));
%num of trial in each bin should be about equal
      bins_total = ves_bins_total + vis_bins_total + comb_bins_total
end

% rt_vest = [];
% rt_vis = [];
% rt_comb = [];
% rt_all = [];


rt_vest = cell(size(unique_motion_coherence));
rt_vis = cell(size(unique_motion_coherence));
rt_comb = cell(size(unique_motion_coherence));
rt_all = cell(size(unique_motion_coherence));


for m = 1:length(unique_motion_coherence)
if length(unique_stim_type) == 3 %i.e. all 3 conditions/modalities were run on the monkey
    for k = 1:length(unique_stim_type)

            trials_select = logical( (stim_type == unique_stim_type(k)) &  ...
                motion_coherence == unique_motion_coherence(m) );
            if unique_stim_type(k) == 1
                if m == 1
                    rt_vest{m} = nonzeros(trials_select .* stim_durations);
                elseif m == 2 
                    rt_vest{m} = rt_vest{1}; %repeating this to keep the code from crashing as vestibular runs will not have 2 different coherences
                end
            elseif unique_stim_type(k) == 2
                rt_vis{m} = nonzeros(trials_select .* stim_durations);
            elseif unique_stim_type(k) == 3
                rt_comb{m} = nonzeros(trials_select .* stim_durations);
            end

    end
elseif length(unique_stim_type) == 1 % i.e. only one condition/modality was run
    trials_select = logical( (stim_type==unique_stim_type) );
    if unique_stim_type(k) == 1
        rt_vest{m} = nonzeros(trials_select .* stim_durations);
    elseif unique_stim_type(k) == 2
        rt_vis{m} = nonzeros(trials_select .* stim_durations);
    elseif unique_stim_type(k) == 3
        rt_comb{m} = nonzeros(trials_select .* stim_durations);
    end
elseif length(unique_stim_type) == 2 % i.e. only 2 condtions were run on the monkey
    for k = 1:length(unique_stim_type)
        trials_select = logical( (stim_type==unique_stim_type(k)) );
        if unique_stim_type(k) == 1
            rt_vest{m} = nonzeros(trials_select .* stim_durations);
        elseif unique_stim_type(k) == 2
            rt_vis{m} = nonzeros(trials_select .* stim_durations);
        elseif unique_stim_type(k) == 3
            rt_comb{m} = nonzeros(trials_select .* stim_durations);
        end
    end
end
end

% rt_all = non_zero

% figure; 
% hist(rt_vest


% correct_rate = [];
% for k = 1:length(unique_stim_type)    
%      for i = 1:length(unique_heading)
%          trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k))  ) ;
%          correct_trials = (trials_select & (choice == RIGHT) );
%          correct_rate(k,i) = 1*sum(correct_trials) / sum(trials_select);  
%          fit_data_psycho_cum{k}(i, 1) = unique_heading( i );  
%          fit_data_psycho_cum{k}(i, 2) = correct_rate(k,i);
%          fit_data_psycho_cum{k}(i, 3) = sum(trials_select);    
%      end   
%      trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
%      trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_stim_type(k))  ) ;
%      trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k)) ); %exclude 0 headings
%      correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);
%      
%      aa = find(fit_data_psycho_cum{k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
%      fit_valid{k}(:,1) = fit_data_psycho_cum{k}(aa,1); 
%      fit_valid{k}(:,2) = fit_data_psycho_cum{k}(aa,2);
%      fit_valid{k}(:,3) = fit_data_psycho_cum{k}(aa,3);
% end
% 

% % Removing NaNs so as to avoid fitting isues with with pfit
row_NaNs_all = cell(size(fit_data_psycho_cum)); %initialization; holding place for all heading directions with NaNs
column_NaNs_all = cell(size(fit_data_psycho_cum)); %initialization; holding place for all heading directions with NaNs

for k = 1:length(unique_stim_type)
    for j = 1:num_bins
        for m = 1:length(unique_motion_coherence)
            temp_fit_data_psycho =  fit_data_psycho_cum{k, j, m};
            temp_fit_data_psycho_NaNs = isfinite(temp_fit_data_psycho);
            [row_NaNs, column_NaNs, value_NaNs] = find(temp_fit_data_psycho_NaNs == 0);
            if sum(value_NaNs) == 0
                %do noting and coninue on with analysis
            else
                disp(['removing NaNs from condition = ',  num2str(k), ', bin = ', num2str(j), ', condition = ', num2str(j), ' data...'])
                %             debugging
                %             if k == 1 && j == 2
                %                 disp('break point')
                %             end

                row_NaNs_all{k, j} = row_NaNs;
                column_NaNs_all{k, j} = column_NaNs;

                diff_row_NaNs = diff(row_NaNs); % shows which NaNs are right next to each other.
                diff_row_NaNs = cat(1, diff_row_NaNs, 0); %add a place holder number at the end of vector

                for i = 1:length(value_NaNs)
                    %original code set NaN values to 1, but now I'm just taking
                    %the averge of the direction before and after that heading
                    %and just plot it as a circle with a cross in it.
                    %               temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) = 0;
                    %               temp_fit_data_psycho(row_NaNs(i), 3) = 1; %setting number of observations equal to one for the NaN situation

                    %          for debugging
                    %           k
                    %           j
                    %           i

                    %calculating averages for the position with NaN
                    if row_NaNs(i) == 1 %i.e. the easiest heading
                        temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) =  0;
                        temp_fit_data_psycho(row_NaNs(i), 3) = 1;
                    elseif row_NaNs(i) == length(unique_heading) % the esasiest heading
                        temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) =  1;
                        temp_fit_data_psycho(row_NaNs(i), 3) = 1;
                    elseif diff_row_NaNs(i) == 1 % if there are NaNs side by side, set the first NaN to a zero or one
                        unique_heading_length = length(unique_heading);
                        heading_middle = unique_heading_length/2 + 0.5;
                        if row_NaNs(i) < heading_middle
                            temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) =  0; %if its a leftward direction, set the value to 0
                        elseif row_NaNs(i) == heading_middle && mod(unique_heading_length, 2) == 1
                            temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) = 0.5; %i.e. if a heading is in the middle of the range and the range is odd, set it to 0.5
                        elseif row_NaNs(i) > heading_middle
                            temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) =  1; % if its a rightward direction, set the value to 1
                        end
                        temp_fit_data_psycho(row_NaNs(i), 3) = 1;
                    else
                        temp_fit_data_psycho(row_NaNs(i), column_NaNs(i)) =  mean([ temp_fit_data_psycho(row_NaNs(i) - 1, column_NaNs(i)), temp_fit_data_psycho(row_NaNs(i) + 1, column_NaNs(i)) ]);
                        temp_fit_data_psycho(row_NaNs(i), 3) = round( mean([ temp_fit_data_psycho(row_NaNs(i) - 1, 3), temp_fit_data_psycho(row_NaNs(i) + 1, 3) ]) );
                    end
                end
                fit_data_psycho_cum{k, j, m} = temp_fit_data_psycho;
                %           disp(['replaced ' num2str(sum(value_NaNs)) ' with zero'])
                disp(['replaced ' num2str(sum(value_NaNs)) ' NaNs with average values'])
            end
        end

    end
end




% replacing all NaNvalues in the correct_rate matrix
for k = 1:length(unique_stim_type)
    for i = 1:length(unique_heading)
        for j = 1:num_bins
            for m = 1:length(unique_motion_coherence)
                correct_rate(k, j, i, m) = fit_data_psycho_cum{k, j, m}(i, 2);
            end
        end
    end
end




%%%%%% use Wichman's MLE method to estimate threshold and bias
for k = 1:length(unique_stim_type)
    for j = 1:num_bins
        for m = 1:length(unique_motion_coherence)
%             k
%             j
%             m
%             if k == 2 && j == 2 && m == 2 
%                 disp('break bitch!!')
%             end

            wichman_psy = pfit(fit_data_psycho_cum{k, j, m},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
            Thresh_psy{k, j, m} = wichman_psy.params.est(2);
            Bias_psy{k, j, m} = wichman_psy.params.est(1);
            psy_perf{k, j, m} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];

            % %   similar way to fit data
            %     [bb,tt] = cum_gaussfit_max1(fit_valid{k});
            %     Thresh_psy{k} = tt;
            %     Bias_psy{k} = bb;
            %     psy_perf{k} =[bb,tt];
        end
    end
end

% added by GY 12-04-07
% now this is the prediction when there are three stimuli conditions 
% if length(unique_stim_type) ==3
%     Thresh_pred = sqrt( Thresh_psy{1}^2*Thresh_psy{2}^2/(Thresh_psy{1}^2+Thresh_psy{2}^2) );
% end
% % this is the output, you can use it for plot of example cells
% xi = min(unique_heading) : 0.1 : max(unique_heading);
% for k = 1:length(unique_stim_type)  
%     yi{k} = cum_gaussfit(psy_perf{k}, xi);
% end
% if length(unique_stim_type) ==3
%     yi_pred = cum_gaussfit([Bias_psy{3},Thresh_pred], xi); % smoothed line for prediction with PSE at actual combined condition
% end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot psychometric, neurometric, CP over time
% % run the slide threshold over time, see whether performance fluctuate across time
% % % commented out at the behest of Yong --- Tunde
% if overtimeplot == 1
%     span = 6;  % calculate threshod every ? repeats;
%     slide = 1;  % slide threshod with increment of ? repeats;
%     BegTrial_shift = BegTrial;
%     EndTrial_shift = BegTrial_shift + span*one_repetition-1;
%     n=0;
%     while EndTrial_shift <= EndTrial    
%         n = n + 1;
%         select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
%         stim_type_shift = temp_stim_type( select_trials_shift );
%         mask_status_shift = temp_mask_status( select_trials_shift );
%         heading_shift = temp_heading( select_trials_shift );
%         unique_stim_type_shift = munique(stim_type_shift');
%         unique_mask_status_shift = temp_mask_status( select_trials_shift );
%         unique_heading_shift = munique(heading_shift');
%         total_trials_shift = temp_total_trials( select_trials_shift);
%         if find(unique_mask_status == 1) > 1        
%             condition_shift = mask_status_shift;
%             unique_condition_shift = unique_mask_status_shift;
%         else
%             condition_shift = stim_type_shift;
%             unique_condition_shift = unique_stim_type_shift;
%         end
%         for k = 1:length(unique_condition_shift)
%             for i = 1:length(unique_heading)
%                  trials_shift =logical( (heading_shift == unique_heading(i)) & (condition_shift == unique_condition_shift(k)) ) ;
%                  correct_trials_shift = (trials_shift & (total_trials_shift == CORRECT) );
%                  % make 'S' curve by using the rightward choice for y-axis
%                  if sum(trials_shift)>0
%                      if ( unique_heading(i) < 0 )
%                          correct_rate_shift(i) = 1 - 1*sum(correct_trials_shift) / sum(trials_shift); 
%                      else
%                          correct_rate_shift(i) = 1*sum(correct_trials_shift) / sum(trials_shift); 
%                      end    
%                  end
%                  Trials_num(i) = sum(trials_shift);
%             end
%             aa = find(correct_rate_shift >-1 );
%             for j = 1:length(aa)
%                  fit_data_psycho_cum_shift{k}(j, 1) = fit_data_psycho_cum{k}(aa(j), 1);  
%                  fit_data_psycho_cum_shift{k}(j, 2) = correct_rate_shift(aa(j));
%                  fit_data_psycho_cum_shift{k}(j, 3) = Trials_num(aa(j));
%             end
%             % this fixes a strange error: cum_gaussfit/pfit sometimes fail when pct choices are all 0's or 1's -CRF 8-13-08
%             if fit_data_psycho_cum_shift{k}(:,2)==0 | fit_data_psycho_cum_shift{k}(:,2)==1
%                 fit_data_psycho_cum_shift{k}(fit_data_psycho_cum_shift{k}==0) = 0.001;
%                 fit_data_psycho_cum_shift{k}(fit_data_psycho_cum_shift{k}==1) = 0.999;
%             end
%             [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{k}); % to save time, use a different fit method
%             psy_thresh_shift(k,n) = tt;
% %             wichman_psy = pfit(fit_data_psycho_cum_shift{k},'plot_opt','plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
% %             psy_thresh_shift(k,n) = wichman_psy.params.est(2);
%         end   
%         BegTrial_shift = BegTrial_shift + slide*one_repetition;
%         EndTrial_shift = BegTrial_shift + span*one_repetition-1;
%     end
% end



% % % % plot psychometric function here
% symbo{1,1} = 'bd';    symbo{1,2} = 'bs';    symbo{1,3} = 'bo';
% symbo{2,1} = 'rd';    symbo{2,2} = 'rs';    symbo{2,3} = 'ro';
% symbo{3,1} = 'gd';    symbo{3,2} = 'gs';    symbo{3,3} = 'go';
% fitline{1,1} = 'b-';    fitline{1,2} = 'b.';    fitline{1,3} = 'b--'; 
% fitline{2,1} = 'r-';    fitline{2,2} = 'r.';    fitline{2,3} = 'r--'; 
% fitline{3,1} = 'g-';    fitline{3,2} = 'g.';    fitline{3,3} = 'g--'; 
% 
% figure(2);
% set(2,'Position', [200,20 700,600], 'Name', 'Heading Discrimination-Vestibular');
% axes('position',[0.2,0.3, 0.6,0.5] );
% % % % fit data with cumulative gaussian and plot both raw data and fitted curve
% legend_txt = [];
% 
% xi = min(unique_heading) : 0.1 : max(unique_heading);
% for k = 1:length(unique_stim_type)
%     for j = 1:num_bins
%         plot(unique_heading, reshape(correct_rate(k, j, :), 1, 8), symbo{k, j},  xi, cum_gaussfit(psy_perf{k, j}, xi),  fitline{k, j} );
%         xlabel('Heading Angles');
%         ylim([0,1]);
%         ylabel('Rightward Choices');
%         set(gca, 'YTickMode','auto');
%         set(gca, 'xTickMode','auto');
%         hold on;
%         legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
%         legend_txt{k*2} = [''];
%     end
% end




% plot psychometric function here
% symbo{1,1} = 'bd';    symbo{1,2} = 'bs';    symbo{1,3} = 'bo';
% symbo{2,1} = 'rd';    symbo{2,2} = 'rs';    symbo{2,3} = 'ro';
% symbo{3,1} = 'gd';    symbo{3,2} = 'gs';    symbo{3,3} = 'go';
% fitline{1,1} = 'b-';    fitline{1,2} = 'b-';    fitline{1,3} = 'b-'; 
% fitline{2,1} = 'r-';    fitline{2,2} = 'r-';    fitline{2,3} = 'r-'; 
% fitline{3,1} = 'g-';    fitline{3,2} = 'g-';    fitline{3,3} = 'g-'; 


% symbo{1,1} = 'bo';    symbo{1,2} = 'bo';    symbo{1,3} = 'bo';
% symbo{2,1} = 'ro';    symbo{2,2} = 'ro';    symbo{2,3} = 'ro';
% symbo{3,1} = 'go';    symbo{3,2} = 'go';    symbo{3,3} = 'go';
% fitline{1,1} = 'b-';    fitline{1,2} = 'b-';    fitline{1,3} = 'b-'; 
% fitline{2,1} = 'r-';    fitline{2,2} = 'r-';    fitline{2,3} = 'r-'; 
% fitline{3,1} = 'g-';    fitline{3,2} = 'g-';    fitline{3,3} = 'g-'; 




symbo_NaN{1,1} = 'bx';
symbo_NaN{2,1} = 'rx';
symbo_NaN{3,1} = 'gx';

symbo{1,1} = 'bo';
symbo{2,1} = 'ro';
symbo{3,1} = 'go';

fitline{1,1} = 'b-';
fitline{2,1} = 'r-';
fitline{3,1} = 'g-';

for i = 2:num_bins
    symbo{1,i} = 'bo';
    symbo{2,i} = 'ro';
    symbo{3,i} = 'go';

    fitline{1,i} = 'b-';
    fitline{2,i} = 'r-';
    fitline{3,i} = 'g-';
end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% starting plot of figurs; This section is for figs 2 and/or 4


fig_list = cell(size(unique_motion_coherence));
fig_handle_list = cell(size(unique_motion_coherence));
if length(unique_motion_coherence) == 1 
    fig2_handle = figure(2);
    fig_handle_list{1} = 'fig2_handle';
    fig_list{1} = 'figure(2)';
elseif length(unique_motion_coherence) == 2
    fig2_handle = figure(2);
    fig_list{1} = 'figure(2)';
    fig_handle_list{1} = 'fig2_handle';
           
    fig4_handle = figure(4);
    fig_list{2} = 'figure(4)';
    fig_handle_list{2} = 'fig4_handle';
end

for y = 1:length(unique_motion_coherence)
%     fig2_handle = figure(2);
    % set(2,'Position', [1 1024*1/2 1280*1/2-10 1024*1/2-100], 'Name', 'Heading Discrimination-Vestibular');
    set(eval([fig_list{y}]),'Position', [1 54 1280 1024-150], 'Name', 'Heading Discrimination-Vestibular');
    % axes('position',[0.2,0.3, 0.6,0.5] );
    % fit data with cumulative gaussian and plot both raw data and fitted curve
    % legend_txt = [];
    xi = min(unique_heading) : 0.1 : max(unique_heading);
    axes('position',[0.2,0.95, 0.6,0.05] );
    text(0.5, 0.2, ['Psychometric Plots ' FILE  '; Coherence: ' num2str(unique_motion_coherence(y)) '%'])
    axis off
end
            

for k = 1:length(unique_stim_type)
    for j = 1:num_bins
        for n = 1:length(unique_motion_coherence)
            eval([fig_list{n}])
            subplot(floor(num_bins/2) + 1, 2, j)
            plot(unique_heading, reshape(correct_rate(k, j, :, n), 1, length(unique_heading)), symbo{unique_stim_type(k), j}, xi, cum_gaussfit(psy_perf{k, j, n}, xi),  fitline{unique_stim_type(k), j} );
            hold on
            if ~isempty(row_NaNs_all{k, j}) % if NaNs were replaced then plot X's
                for m = 1:size(row_NaNs_all{k, j}, 1)
                    plot( fit_data_psycho_cum{k, j, n}(row_NaNs_all{k, j}(m), 1), fit_data_psycho_cum{k, j, n}(row_NaNs_all{k, j}(m), column_NaNs_all{k, j}(m)), symbo_NaN{unique_stim_type(k)} )
                end
            end
            xlabel('Heading Angles');
            ylim([0,1]);
            ylabel('Rightward Choices');
            %         title(['Stimulus duration: ' num2str(bin_sections_mid(j)) ' ms']) % middle of bins
            title([ 'Stimulus duration bin range: ' num2str(bin_sections(j)) 'ms to ' num2str(bin_sections(j + 1)) 'ms' ])
            set(gca, 'YTickMode','auto');
            set(gca, 'xTickMode','auto');
            hold on;
            legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
            legend_txt{k*2} = [''];
        end
    end
end




% this is just for labelling the bins
for i = 1:length(bin_sections_mid)
    bin_sections_mid_str{i} = num2str(bin_sections_mid(i));
end
for m = 1:length(unique_motion_coherence)
    % last plot of all the trials in psychometric
    eval([fig_list{m}])
    subplot(floor(num_bins/2) + 1, 2, num_bins + 1)
    % hist(stim_durations, num_bins)
    if bin_type == 0
        hist(stim_durations, bin_sections_mid)
        title(['stimulus time distributions for seesion: ' FILE])
        xlabel('time (ms)')
    elseif bin_type == 1
        frequency_counts = histc(stim_durations, bin_sections);
        bar(frequency_counts(1:length(frequency_counts) - 1), 1)
        set(gca, 'xticklabel',bin_sections_mid_str)
        title(['Psychometric across all bins: ' FILE])
        xlabel('Heading Angles')
    end
end



    
%     subplot(224)
%     hist(stim_durations, t4_relative_bin_range)
%     xlabel('time (ms)')
% 
% 
% 
%     xlim( [0,50] );
%     ylim( [-5,9] );
%     text(0, 10, FILE);
%     text(15,10,'coherence =');
%     text(30,10,'repeats =');
%     text(45,10,'maskradius =');
%     text(25,10,num2str(unique_motion_coherence) );
%     text(40,10,num2str(repetition) );
%     text(55,10,num2str(unique_mask_radius) );
%     text(0,8, '           u              sigma                  u           sigma                     u              sigma      ');
% 
%     for k = 1:length(unique_stim_type)
%         for j = 1:num_bins
%             if k == 1
%     %             text(0,8-(k*0 + j), num2str(unique_stim_type(k)));
%                 text(0,8-(k*0 + j),num2str(Bias_psy{k, j}) );
%                 text(10,8-(k*0 + j),num2str(Thresh_psy{k, j}) );
%             elseif k == 2
%     %             text(0,8-(k*1.5 + j), num2str(unique_stim_type(k)));
%                 text(20,8-(k*1.5 + j),num2str(Bias_psy{k, j}) );
%                 text(30,8-(k*1.5 + j),num2str(Thresh_psy{k, j}) );
%             elseif k ==3
%     %             text(0,8-(k*2 + j), num2str(unique_stim_type(k)));
%                 text(40,8-(k*2 + j),num2str(Bias_psy{k, j}) );
%                 text(50,8-(k*2 + j),num2str(Thresh_psy{k, j}) );
%             end
%         end
%     end
% 
% 
% 
% 
% 
% 
%         text(0,8-(0*0 + 1), 'ves');
%         text(0,8-(0*0 + 2), 'vis');
%         text(0,8-(0*0 + 3), 'comb')
%     
%     
%     
%     for k = 1:length(unique_stim_type)
%         for j = 1:num_bins
%             if k == 1
%     %             text(0,8-(k*0 + j), num2str(unique_stim_type(k)));
%                 text(5,8-(k*0 + j),num2str(Bias_psy{k, j}) );
%                 text(13,8-(k*0 + j),num2str(Thresh_psy{k, j}) );
%             elseif k == 2
%     %             text(0,8-(k*1.5 + j), num2str(unique_stim_type(k)));
%                 text(25,8-(k*0 + j),num2str(Bias_psy{k, j}) );
%                 text(33,8-(k*0 + j),num2str(Thresh_psy{k, j}) );
%             elseif k ==3
%     %             text(0,8-(k*2 + j), num2str(unique_stim_type(k)));
%                 text(45,8-(k*0 + j),num2str(Bias_psy{k, j}) );
%                 text(53,8-(k*0 + j),num2str(Thresh_psy{k, j}) );
%             end
%         end
%     end
%     axis off;
% 
% 
% 
%         text(0,8-(0*0 + 1), 'ves');
%         text(0,8-(0*0 + 2), 'vis');
%         text(0,8-(0*0 + 3), 'comb')
% 
% 
%     for k = 1:length(unique_stim_type)
%         for j = 1:num_bins
%             if unique_stim_type(k) == 1 %vestibular
%                 subplot(floor(num_bins/2) + 1, 2, j)
%     %             text(0,8-(k*0 + j), num2str(unique_stim_type(k)));
%                 text(0.5, 0.5, [{'\color{black}mu' }]  );
%                 text(20, 0.5, [{'\color{black}Sigma'  }]  );
%                 text(0.5, 0.45, [{'\color{blue}' num2str(Bias_psy{k, j}) }]  );
%                 text(20, 0.45, [{'\color{blue}' num2str(Thresh_psy{k, j}) }]  );
%             elseif unique_stim_type(k) == 2 %visual
%     %             text(0,8-(k*1.5 + j), num2str(unique_stim_type(k)));
%                 subplot(floor(num_bins/2) + 1, 2, j)
%                 text(0.5, 0.5, [{'\color{black}mu' }]  );
%                 text(20, 0.5, [{'\color{black}Sigma'  }]  );
%                 text(0.5, 0.4, [{'\color{red}' num2str(Bias_psy{k, j}) }]  );
%                 text(20, 0.4, [{'\color{red}' num2str(Thresh_psy{k, j}) }]  );
%             elseif unique_stim_type(k) == 3 %combined
%     %             text(0,8-(k*2 + j), num2str(unique_stim_type(k)));
%                 subplot(floor(num_bins/2) + 1, 2, j)
%                 text(0.5, 0.5, [{'\color{black}mu' }]  );
%                 text(20, 0.5, [{'\color{black}Sigma'  }]  );
%                 text(0.5, 0.35, [{'\color{green}' num2str(Bias_psy{k, j}) }]  );
%                 text(20, 0.35, [{'\color{green}' num2str(Thresh_psy{k, j}) }]  );
%             end
%         end
%     end
% 
% 
%     for k = 1:length(unique_stim_type)
%         for j = 1:num_bins
%             if unique_stim_type(k) == 1 %vestibular
%                 subplot(floor(num_bins/2) + 1, 2, j)
%     %             text(0,8-(k*0 + j), num2str(unique_stim_type(k)));
%                 text(-20, 0.5, [{'\color{black}mu' }]  );
%                 text(.5, 0.5, [{'\color{black}Sigma'  }]  );
%                 text(-20, 0.45, [{'\color{blue}' num2str(Bias_psy{k, j}) }]  );
%                 text(.5, 0.45, [{'\color{blue}' num2str(Thresh_psy{k, j}) }]  );
%             elseif unique_stim_type(k) == 2 %visual
%     %             text(0,8-(k*1.5 + j), num2str(unique_stim_type(k)));
%                 subplot(floor(num_bins/2) + 1, 2, j)
%                 text(20, 0.5, [{'\color{black}mu' }]  );
%                 text(.5, 0.5, [{'\color{black}Sigma'  }]  );
%                 text(-20, 0.4, [{'\color{red}' num2str(Bias_psy{k, j}) }]  );
%                 text(.5, 0.4, [{'\color{red}' num2str(Thresh_psy{k, j}) }]  );
%             elseif unique_stim_type(k) == 3 %combined
%     %             text(0,8-(k*2 + j), num2str(unique_stim_type(k)));
%                 subplot(floor(num_bins/2) + 1, 2, j)
%                 text(-20, 0.5, [{'\color{black}mu' }]  );
%                 text(.5, 0.5, [{'\color{black}Sigma'  }]  );
%                 text(-20, 0.35, [{'\color{green}' num2str(Bias_psy{k, j}) }]  );
%                 text(.5, 0.35, [{'\color{green}' num2str(Thresh_psy{k, j}) }]  );
%             end
%         end
%     end
%     
%     
%     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


clear unique_heading_max unique_heading_min diff_max_min_heading mid_max_min_heading mid_plus_max_min_heading
%clear variables as they are defined earlier above for figure(6)
    unique_heading_max = max(unique_heading);
    unique_heading_min = min(unique_heading);
    diff_max_min_heading = unique_heading_max - unique_heading_min;
    mid_max_min_heading = diff_max_min_heading/2 + unique_heading_min;
    mid_plus_max_min_heading = (diff_max_min_heading*3/4) + unique_heading_min;


    for k = 1:length(unique_stim_type)
        for j = 1:num_bins
            for m = 1:length(unique_motion_coherence)
                eval([fig_list{m}])
                if unique_stim_type(k) == 1 %vestibular
                    subplot(floor(num_bins/2) + 1, 2, j)
                    %             text(0,8-(k*0 + j), num2str(unique_stim_type(k)));
                    text(mid_max_min_heading, 0.5, [{'\color{black}mu' }]  );
                    text(mid_plus_max_min_heading, 0.5, [{'\color{black}Sigma'  }]  );
                    text(mid_max_min_heading, 0.45, [{'\color{blue}' num2str(Bias_psy{k, j, m}) }]  );
                    text(mid_plus_max_min_heading, 0.45, [{'\color{blue}' num2str(Thresh_psy{k, j, m}) }]  );
                elseif unique_stim_type(k) == 2 %visual
                    %             text(0,8-(k*1.5 + j), num2str(unique_stim_type(k)));
                    subplot(floor(num_bins/2) + 1, 2, j)
                    text(mid_max_min_heading, 0.5, [{'\color{black}mu' }]  );
                    text(mid_plus_max_min_heading, 0.5, [{'\color{black}Sigma'  }]  );
                    text(mid_max_min_heading, 0.4, [{'\color{red}' num2str(Bias_psy{k, j, m}) }]  );
                    text(mid_plus_max_min_heading, 0.4, [{'\color{red}' num2str(Thresh_psy{k, j, m}) }]  );
                elseif unique_stim_type(k) == 3 %combined
                    %             text(0,8-(k*2 + j), num2str(unique_stim_type(k)));
                    subplot(floor(num_bins/2) + 1, 2, j)
                    text(mid_max_min_heading, 0.5, [{'\color{black}mu' }]  );
                    text(mid_plus_max_min_heading, 0.5, [{'\color{black}Sigma'  }]  );
                    text(mid_max_min_heading, 0.35, [{'\color{green}' num2str(Bias_psy{k, j, m}) }]  );
                    text(mid_plus_max_min_heading, 0.35, [{'\color{green}' num2str(Thresh_psy{k, j, m}) }]  );
                    %             text(mid_max_min_heading, 0.35, [{'\color{black}' num2str(Bias_psy{k, j}) }]  );
                    %             text(mid_plus_max_min_heading, 0.35, [{'\color{black}' num2str(Thresh_psy{k, j}) }]  );
                end
            end
        end
    end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
% end of plots of figs 2 and 4 (if using 2 coherences) OR fig 2 if using
% just one coherence.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
high_res_bins = 1:50:2000;
% fig3_handle = figure(3);
if length(unique_motion_coherence) == 1 
    fig3_handle = figure(3);
    fig_handle_list{2} = 'fig3_handle';
    fig_list{2} = 'figure(3)';
elseif length(unique_motion_coherence) == 2
    fig3_handle = figure(3);
    fig_list{3} = 'figure(3)';
    fig_handle_list{3} = 'fig3_handle';
           
end





set(fig3_handle, 'Name', 'Higher resolution bins historgram', 'Position',[1 54 1280 1024-150])
% axes('position',[0.2,0.95, 0.6,0.05] );
axes('position',[0.01,0.95, 0.99,0.05] );
txt_nums = [0.05, 0.5];
for m = 1:length(unique_motion_coherence)
    % text(0.3, 0.5, ['Low and High resolution bins, Histogram plots ' FILE, '; Coherence: ' num2str(unique_motion_coherence(m)) '%'])
    text(txt_nums(m), 0.5, ['High resolution bins, Histogram plots ' FILE, '; Coherence: ' num2str(unique_motion_coherence(m)) '%'])

%     text(0.55, 0.5, ['High resolution bins, Histogram plots ' FILE, '; Coherence: ' num2str(unique_motion_coherence(m)) '%'])
end
axis off



% dividing all stimulus durations by coherence
stim_durations_coherence = cell(size(unique_motion_coherence));
for m = 1:length(unique_motion_coherence)
    trials_select = logical( motion_coherence == unique_motion_coherence(m) );
    stim_durations_coherence{m} = nonzeros(trials_select .* stim_durations);
end

%Modified by Adhira - 01/05/11
for m = 1:length(unique_motion_coherence)
for k = 1:length(unique_stim_type)
    % % % % % All trials
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %     eval([fig_list{m + length(unique_motion_coherence)}]) %bring fig3 or fig4 to the front.
    figure(3)
    %     if m == 1
    subplot(421)
    %     elseif m == 2
    %         subplot(422)
    %     end

    hist(stim_durations_coherence{m}, high_res_bins)
    title(['Total reaction time bin : High resolution bins; Mean+/-ste = ' num2str(round(mean(stim_durations))) '+/-' num2str(round( std(stim_durations)/sqrt(length(stim_durations)) ))])
    xlabel('time (ms)')
    % xlim([100, 1100]);
    max_stim_dur = max(stim_durations_coherence{m});
    min_stim_dur = min(stim_durations_coherence{m});
    xlim([min_stim_dur, max_stim_dur])
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k')

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % % vestibular and vest, vis, combined trials
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if (length(unique_stim_type) == 1 && unique_stim_type == 1) || (length(unique_stim_type) ==2 && unique_stim_type(k) ==1) || (length(unique_stim_type) == 3)
        %         if m == 1
        subplot(423)
        %         elseif m == 2
        %             subplot(424)
        %         end
        hist(rt_vest{m}, high_res_bins)
        title(['Vestibular: High Resolutions bins; Mean+/-ste = ' num2str(round(mean(rt_vest{m}))) '+/-' num2str(round( std(rt_vest{m})/sqrt(length(rt_vest{m})) )) ])
        xlabel('time (ms)')
        %     xlim([100, 1100]);
        max_stim_dur = max(stim_durations_coherence{m});
        min_stim_dur = min(stim_durations_coherence{m});
        xlim([min_stim_dur, max_stim_dur])
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b')

        % hist(rt_comb( find(rt_comb > bin_sections(1) & rt_comb < bin_sections(2)) ), 20)
        % title([ 'Reaction time bin range: ' num2str(bin_sections(1)) 'ms to ' num2str(bin_sections(2)) 'ms' ])
    end

    % % % v
    %     if (length(unique_stim_type) == 2 && unique_stim_type == 1) || (length(unique_stim_type) == 2)
    %     if (length(unique_stim_type) == 2)
    %         subplot(423)
    %         if bin_type == 0
    %             hist(rt_vest, bin_sections_mid)
    %             xlim([100, 1100]);
    %             h = findobj(gca,'Type','patch');
    %             set(h,'FaceColor','b')
    %         elseif bin_type == 1
    %             frequency_counts = histc(rt_vest, bin_sections);
    %             bar(frequency_counts(1:length(frequency_counts) - 1), 1, 'b')
    %             set(gca, 'xticklabel',bin_sections_mid_str)
    %         end
    %         title('Vestibular: Low Resolution bins')
    %         xlabel('time (ms)')
    %
    %         subplot(424)
    %         hist(rt_vest, high_res_bins)
    %         title('Vestibular: High Resolutions bins')
    %         xlabel('time (ms)')
    %         xlim([100, 1100]);
    %         h = findobj(gca,'Type','patch');
    %         set(h,'FaceColor','b')
    %
    %         % hist(rt_comb( find(rt_comb > bin_sections(1) & rt_comb < bin_sections(2)) ), 20)
    %         % title([ 'Reaction time bin range: ' num2str(bin_sections(1)) 'ms to ' num2str(bin_sections(2)) 'ms' ])
    %     end


    if (length(unique_stim_type) == 1 && unique_stim_type == 2) || (length(unique_stim_type) ==2 && unique_stim_type(k) ==2) || (length(unique_stim_type) == 3)
        %         if m == 1
        subplot(425)
        %         elseif m == 2
        %             subplot(426)
        %         end
        hist(rt_vis{m}, high_res_bins)
        title(['Visual: High Resolution bins; Mean+/-ste = ' num2str(round(mean(rt_vis{m}))) '+/-' num2str(round( std(rt_vis{m})/sqrt(length(rt_vis{m})) ))])
        xlabel('time (ms)')
        %     xlim([100, 1100]);
        max_stim_dur = max(stim_durations_coherence{m});
        min_stim_dur = min(stim_durations_coherence{m});
        xlim([min_stim_dur, max_stim_dur])
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r')

        % hist( rt_comb( find(rt_comb > bin_sections(2) & rt_comb < bin_sections(3)) ), 20)
        % title([ 'Reaction time bin range: ' num2str(bin_sections(2)) 'ms to ' num2str(bin_sections(3)) 'ms' ])
    end


    if (length(unique_stim_type) == 1 && unique_stim_type == 3) || (length(unique_stim_type) ==2 && unique_stim_type(k) ==3) || (length(unique_stim_type) == 3)
        %         if m == 1
        subplot(427)
        %         elseif m == 2
        %             subplot(428)
        %         end

        hist(rt_comb{m}, high_res_bins)
        title(['Combined: High Resolution bins; Mean+/-ste = ' num2str(round(mean(rt_comb{m}))) '+/-' num2str(round( std(rt_comb{m})/sqrt(length(rt_comb{m})) ))])
        xlabel('time (ms)')
        %     xlim([100, 1100]);
        max_stim_dur = max(stim_durations_coherence{m});
        min_stim_dur = min(stim_durations_coherence{m});
        xlim([min_stim_dur, max_stim_dur])
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','g')
    end

end
end


% % % % % % % % % comment end for batch_gui

% % if (length(unique_stim_type) == 2 && unique_stim_type == 3) || (length(unique_stim_type) == 2)
%     subplot(427)
%     if bin_type == 0
%         hist(rt_comb, bin_sections_mid)
%         xlim([100, 1100]);
%         h = findobj(gca,'Type','patch');
%         set(h,'FaceColor','g')
%     elseif bin_type == 1
%         frequency_counts = histc(rt_comb, bin_sections);
%         bar(frequency_counts(1:length(frequency_counts) - 1), 1, 'g')
%         set(gca, 'xticklabel',bin_sections_mid_str)
%     end
%     title('Combined: Low Resolusion bins')
%     xlabel('time (ms)')
% 
%     % hist( rt_comb(  find(rt_comb > bin_sections(3) & rt_comb < bin_sections(4)) ), 20);
%     % title([ 'Reaction time bin range: ' num2str(bin_sections(3)) 'ms to ' num2str(bin_sections(4)) 'ms' ])
% 
%     subplot(428)
%     hist(rt_comb, high_res_bins)
%     title('Combined: High Resolution bins')
%     xlabel('time (ms)')
%     xlim([100, 1100]);
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','g')
% % end



dot_find = findstr('.', FILE);
filename = FILE(1: dot_find - 1);
filename_mat = strcat(filename, '.mat');



% if ~exist('RT_plots', 'dir')
%     mkdir('RT_plots')
% end


rightward_trials = (select_trials & (choice == RIGHT) ); % added this so as to get all the rightward choices the monkey makes for Jan's analysis

% cd RT_plots
% if using_batch == 0
%     if length(unique_motion_coherence) == 1
%         figure(fig2_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
%         saveas(fig2_handle, [FILE(1:(end - 4)), '_psychometric', '.bmp'])
%         figure(fig3_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
%         saveas(fig3_handle, [FILE(1:(end - 4)), '_high_res_bins', '.bmp'])
%         saveas(fig6_handle, [FILE(1:(end - 4)), '_psychometric_single', '.bmp'])
%     elseif length(unique_motion_coherence) == 2
%         figure(fig2_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
%         saveas(fig2_handle, [FILE(1:(end - 4)), '_psychometric1', '.bmp'])
%         figure(fig3_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
%         saveas(fig3_handle, [FILE(1:(end - 4)), '_high_res_bins', '.bmp'])
%         figure(fig4_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
%         saveas(fig4_handle, [FILE(1:(end - 4)), '_psychometric2', '.bmp'])
%         saveas(fig6_handle, [FILE(1:(end - 4)), '_psychometric_single', '.bmp'])
%     end
% elseif using_batch == 1
% %     figure(fig2_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
% %     saveas(fig2_handle, [FILE(1:(end - 4)), '_psychometric', '.bmp'])
%     figure(fig3_handle) %placed here to enable .bmps to be made in the corect order on Jason's machine
%     saveas(fig3_handle, [FILE(1:(end - 4)), '_high_res_bins', '.bmp'])
% end
% save(filename_mat, 'stim_durations', 'rightward_trials', 'heading', 'stim_type', 'motion_coherence', 'pre_stim_delay', 'amplitude')
% clear all

cd ..


% 
% 
% for j = 1:num_bins
%     figure(3)
%     subplot(221)
%     hist(stim_durations)
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','k')
%     subplot(222)
%     hist(rt_vest)
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','b')
%     subplot(223)
%     hist(rt_vis)
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','r')
%     subplot(224)
%     hist(rt_comb)
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','g')
% end















% for k = 1:length(unique_stim_type)
%     text(0,8-k, num2str(unique_stim_type(k))); 
%     text(10,8-k,num2str(Bias_psy{k}) );
%     text(20,8-k,num2str(Thresh_psy{k}) );
% %     text(30,8-k,num2str(correct_proportion(k)) );
% end


% axis off;

% % % plot psycho over time
% if overtimeplot ==1
%     axes('position',[0.2,0.05, 0.6,0.2] );
%     for k = 1:length(unique_stim_type)
%         plot(psy_thresh_shift(k,:), fitline{1,k});
%        % semilogy(psy_thresh_shift(k,:), f{k});
%         hold on;
%         xlabel('Repetition');  
%         ylabel('Threshold');
%         xlim([0, n]);
%       %  ylim( [min(min(psy_thresh_shift(:,:))), max(max(psy_thresh_shift(:,:)))] );   
%     end
% end








% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sprint_txt = ['%s']; 
% for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f\t'];    
% end
% %buff = sprintf(sprint_txt, FILE, fit_data_psycho_cum{1}(:,1),fit_data_psycho_cum{1}(:,2),fit_data_psycho_cum{1}(:,3),fit_data_psycho_cum{2}(:,2),fit_data_psycho_cum{2}(:,3),fit_data_psycho_cum{3}(:,2),fit_data_psycho_cum{3}(:,3) );
% %buff = sprintf(sprint_txt, FILE, unique_motion_coherence, Thresh_psy{1}, Thresh_psy{2},Thresh_psy{3} );
% if length(unique_stim_type)==3
%     buff = sprintf(sprint_txt, FILE, unique_motion_coherence, Bias_psy{1},Bias_psy{2},Bias_psy{3},Thresh_psy{1},Thresh_psy{2},Thresh_psy{3});
% elseif length(unique_stim_type)==2
%     buff = sprintf(sprint_txt, FILE, unique_motion_coherence, Bias_psy{1},Bias_psy{2},Thresh_psy{1},Thresh_psy{2});
% else
%     buff = sprintf(sprint_txt, FILE, unique_motion_coherence, Bias_psy{1},Thresh_psy{1});
% end
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Psychome_combined.dat'];
% % buff = sprintf(sprint_txt, FILE, unique_heading', fit_valid{1}(:,2), fit_valid{2}(:,2),fit_valid{3}(:,2) );
% % outfile = ['Z:\Users\Yong\Inactivationm13c0r635.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t          coherence\t  1_bias\t 2_bias\t 3_bias\t 1_thresh\t 2_thresh\t 3_thresh\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
%---------------------------------------------------------------------------------------
return;
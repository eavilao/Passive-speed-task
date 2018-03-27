%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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

temp_targlum_1 = data.targ_params(TARG_LUM_MULT,:,T1);  %Jing 01/03/2013
temp_targlum_2 = data.targ_params(TARG_LUM_MULT,:,T2);

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');
unique_heading_nonzero = unique_heading(unique_heading~=0);

if length(unique_motion_coherence)==1
   one_repetition = length(unique_heading)*length(unique_stim_type);
else
   one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence)-length(unique_heading); 
end
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
        choiceOverall(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choiceOverall(i) = LEFT;
    else
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
    
    if ((temp_targlum_1(i)==1) && (temp_targlum_2(i)==1))    
        if (sum(events == IN_T1_WIN_CD) > 0)
            choice(i) = RIGHT;
        elseif (sum(events == IN_T2_WIN_CD) > 0)
            choice(i) = LEFT;
        else
            disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
        end
    else 
        choice(i) = 0;
    end
end

correct_rate = [];
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)    
         for i = 1:length(unique_heading)
             if unique_stim_type(k) == 1 % for vestibular condition, take all the data regardless of visual coherence
                 trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & choice ~= 0 ) ;
             else 
                 trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (motion_coherence==unique_motion_coherence(c))& choice ~= 0) ;
             end
             rightward_trials = (trials_select & (choice == RIGHT) );
             rightward_rate = 1*sum(rightward_trials) / sum(trials_select);  
             fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);  
             fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
             fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select); 
         end 

         % the correct rate does not take coherence into account,temporarily 05/29/09
         trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
         trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_stim_type(k))  ) ;
         trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k)) & choice ~= 0); %exclude 0 headings
         correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);

         aa = find(fit_data_psycho_cum{c,k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
         fit_valid{c,k}(:,1) = fit_data_psycho_cum{c,k}(aa,1); 
         fit_valid{c,k}(:,2) = fit_data_psycho_cum{c,k}(aa,2);
         fit_valid{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
         
         %overall correct rate
         trials_rightward_overall = find( (heading > 0) & (choiceOverall==RIGHT) & (stim_type==unique_stim_type(k))) ;
         trials_leftward_overall  = find( (heading < 0) & (choiceOverall==LEFT) & (stim_type==unique_stim_type(k))) ;
         trials_all_overall = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k))); %exclude 0 headings
         correct_proportion_overall(k) = (length(trials_rightward_overall)+length(trials_leftward_overall))/length(trials_all_overall);  
         twoTarg_proportion(k) = length(trials_all)/length(trials_all_overall);
    end
end
%%%%%% use Wichman's MLE method to estimate threshold and bias
for c = 1:length(unique_motion_coherence) % different coherence level
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
        select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) & choice ~= 0);
        stim_type_shift = temp_stim_type( select_trials_shift );
        mask_status_shift = temp_mask_status( select_trials_shift );
        heading_shift = temp_heading( select_trials_shift );
        unique_stim_type_shift = munique(stim_type_shift');
        unique_mask_status_shift = temp_mask_status( select_trials_shift );
        unique_heading_shift = munique(heading_shift');
        total_trials_shift = temp_total_trials( select_trials_shift);
        if find(unique_mask_status == 1) > 1        
            condition_shift = mask_status_shift;
            unique_condition_shift = unique_mask_status_shift;
        else
            condition_shift = stim_type_shift;
            unique_condition_shift = unique_stim_type_shift;
        end
        for k = 1:length(unique_condition_shift)
            for i = 1:length(unique_heading)
                 trials_shift =logical( (heading_shift == unique_heading(i)) & (condition_shift == unique_condition_shift(k)) ) ;
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
                 fit_data_psycho_cum_shift{k}(j, 1) = fit_data_psycho_cum{k}(aa(j), 1);  
                 fit_data_psycho_cum_shift{k}(j, 2) = correct_rate_shift(aa(j));
                 fit_data_psycho_cum_shift{k}(j, 3) = Trials_num(aa(j));
            end
            % this fixes a strange error: cum_gaussfit/pfit sometimes fail when pct choices are all 0's or 1's -CRF 8-13-08
            if fit_data_psycho_cum_shift{k}(:,2)==0 | fit_data_psycho_cum_shift{k}(:,2)==1
                fit_data_psycho_cum_shift{k}(fit_data_psycho_cum_shift{k}==0) = 0.001;
                fit_data_psycho_cum_shift{k}(fit_data_psycho_cum_shift{k}==1) = 0.999;
            end
            [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{k}); % to save time, use a different fit method
            psy_thresh_shift(k,n) = tt;
            psy_bias_shift(k,n) = bb;  % added Bias, CRF 11-5-09
             if n == 11
                 xxx = 0;
             end


        end   
        BegTrial_shift = BegTrial_shift + slide*one_repetition;
        EndTrial_shift = BegTrial_shift + span*one_repetition-1;
    end
end




% plot psychometric function here
symbo{1,1} = 'bo';    symbo{1,2} = 'ro';    symbo{1,3} = 'go'; 
symbo{2,1} = 'b*';    symbo{2,2} = 'm*';    symbo{2,3} = 'g*'; 
fitline{1,1} = 'b-';    fitline{1,2} = 'r-';    fitline{1,3} = 'g-'; 
fitline{2,1} = 'b--';    fitline{2,2} = 'r--';    fitline{2,3} = 'g--'; 

figure(2);
set(2,'Position', [200,100 700,800], 'Name', 'Heading Discrimination-Vestibular');
axes('position',[0.2,0.47, 0.6,0.4] );
% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];
xi = min(unique_heading) : 0.1 : max(unique_heading);
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)      
        %plot(unique_heading, fit_valid{c,k}(:,2), symbo{c,k},  xi, cum_gaussfit(psy_perf{c,k}, xi),  fitline{c,k} ); % Jing 01/18/2012
        plot(fit_valid{c,k}(:,1), fit_valid{c,k}(:,2), symbo{c,k},  xi, cum_gaussfit(psy_perf{c,k}, xi),  fitline{c,k} );
        xlabel('Heading Angles');   
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        hold on;
        legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
        legend_txt{k*2} = [''];
    end
end
% output some text of basic parameters in the figure
axes('position',[0.2,0.83, 0.6,0.15] );
xlim( [0,50] );
ylim( [2,10] );
text(0, 10, FILE);
text(15,10,'coherence =');
text(30,10,'repeats =');
%text(45,10,'maskradius =');
text(25,10,num2str(unique_motion_coherence) );
text(40,10,num2str(repetition) );
%text(55,10,num2str(unique_mask_radius) );
text(10,8, 'u               sigma          2targ-corr rate        corr-rate      2targ-trial%');

for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        text(0,8-k-(c-1)*3, num2str(unique_stim_type(k)));  % non-microstim
        text(10,8-k-(c-1)*3,num2str(Bias_psy{c,k}) );
        text(20,8-k-(c-1)*3,num2str(Thresh_psy{c,k}) );
        text(30,8-k-(c-1)*3,num2str(correct_proportion(k)) );
        text(40,8-k-(c-1)*3,num2str(correct_proportion_overall(k)) );
        text(50,8-k-(c-1)*3,num2str(twoTarg_proportion(k)) );
    end
end

axis off;

xx = [fit_valid{c,k}(:,1), fit_valid{c,k}(:,2)]
yy =[ xi', cum_gaussfit(psy_perf{c,k}, xi)']

% % plot psycho over time
if overtimeplot ==1
    axes('position',[0.2,0.26, 0.6,0.16] );
    for k = 1:length(unique_stim_type)
%         plot(3:n+2,psy_thresh_shift(k,:), fitline{1,k});
%        semilogy(psy_thresh_shift(k,:), f{k});
       
       %%%%%%%%%%%%%%%%%%%%%%%add by 6/18/2010
       m_threshold = psy_thresh_shift(k,:);
       index = find(m_threshold > 150);
       m_threshold(index) = 150;
       plot(3:n+2,m_threshold, fitline{1,k});
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        xlabel('Repetition');  
        ylabel('Threshold');
        xlim([0+3, n+2]);   % Jing 01/18/2012
      %  ylim( [min(min(psy_thresh_shift(:,:))), max(max(psy_thresh_shift(:,:)))] );   
    end
    % added Bias, CRF 11-5-09
    axes('position',[0.2,0.05, 0.6,0.16] );
    for k = 1:length(unique_stim_type)
%         plot(3:n+2,psy_bias_shift(k,:), fitline{1,k});
       % semilogy(psy_thresh_shift(k,:), f{k});
        %%%%%%%%%%%%%%%%%%%%%%%add by 6/18/2010
        m_bias = psy_bias_shift(k,:);
        index = find(m_bias > 100)
        m_bias(index) = 100;
        index = find(m_bias < -100)
        m_bias(index) = -100;
        plot(3:n+2,m_bias, fitline{1,k});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        xlabel('Repetition');  
        ylabel('Bias');
        xlim([0+3, n+2]); % Jing 01/18/2012
      %  ylim( [min(min(psy_thresh_shift(:,:))), max(max(psy_thresh_shift(:,:)))] );   
    end
end
orient tall;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprint_txt = ['%s']; 
for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f\t'];    
end

if length(unique_stim_type)==3
    buff = sprintf(sprint_txt, FILE, unique_motion_coherence, unique_stim_type, Thresh_psy{1,1},Thresh_psy{1,2},Thresh_psy{1,3});
elseif length(unique_stim_type)==2
    buff = sprintf(sprint_txt, FILE, unique_motion_coherence, unique_stim_type,Thresh_psy{1,1},Thresh_psy{1,2});
else
    buff = sprintf(sprint_txt, FILE, unique_motion_coherence, unique_stim_type,Thresh_psy{1,1});
end

outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Psychome_combined.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t          coherence\t  1_bias\t 2_bias\t 3_bias\t 1_thresh\t 2_thresh\t 3_thresh\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
return;
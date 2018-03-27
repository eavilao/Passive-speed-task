%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_FineRotDiscrm_2targ_deg1(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, order);

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
temp_rotamplitude = data.moog_params(ROT_AMPLITUDE,:,MOOG);

temp_targlum_1 = data.targ_params(TARG_LUM_MULT,:,T1);  %Jing 12/17/2012
temp_targlum_2 = data.targ_params(TARG_LUM_MULT,:,T2);

%convert Rotation amp to vel. 
for i=1:length(temp_heading)
    %[velmax maxAcc1s]=calpeakVel(1, 4, abs(temp_rotamplitude(i)));
    velmax = 2*abs(temp_rotamplitude(i))/(0.7+0.5);  % Trapzoid profile. Jing 09/06/2013
 
    temp_heading(i) = velmax *sin(pi*temp_heading(i)/180);
end

trials = 1:length(temp_heading);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial));

stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);
mask_status= temp_mask_status( select_trials );
mask_radius= temp_mask_radius( select_trials );
microstim = temp_microstim(select_trials );

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
                 trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (motion_coherence==unique_motion_coherence(c)) & choice ~= 0 ) ;
             end
             rightward_trials = (trials_select & (choice == RIGHT) );
             rightward_rate = 1*sum(rightward_trials) / sum(trials_select);  
             fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);  
             fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
             fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select); 
         end 

         % the correct rate does not take coherence into account,temporarily 05/29/09
         % correct rate for two target trial
         trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_stim_type(k))) ;
         trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_stim_type(k))) ;
         trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k))& choice ~= 0 ); %exclude 0 headings
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
        %wichman_psy = pfit(fit_valid{c,k}(1:end,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','true','verbose','true');  
        wichman_psy = pfit(fit_valid{c,k}(1:end,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');
        Thresh_psy{c,k} = wichman_psy.params.est(2);
        Bias_psy{c,k} = wichman_psy.params.est(1);
        psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
        
        try
        wichman_psy = pfit(fit_data_psycho_cum{c,k}(2:end-1, :),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');  
        Thresh_psy_dec = wichman_psy.params.est(2);
        Bias_psy_dec = wichman_psy.params.est(1);
        catch,
        Thresh_psy_dec = 0;
        Bias_psy_dec = 0;
        end
        
        Rtotal =  sum((fit_valid{c,k}(:,2) -  mean( fit_valid{c,k}(:,2) ) ).*(fit_valid{c,k}(:,2) -  mean( fit_valid{c,k}(:,2) ) ));
        Rerror =  sum(( fit_valid{c,k}(:,2) -  cum_gaussfit(psy_perf{c,k}, fit_valid{c,k}(:,1)) ).* ( fit_valid{c,k}(:,2) -  cum_gaussfit(psy_perf{c,k}, fit_valid{c,k}(:,1)) ));
        m_Rsquar{c,k} =  1 - Rerror/Rtotal;
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
        %plot(unique_heading, fit_valid{c,k}(:,2), symbo{c,k},  xi, cum_gaussfit(psy_perf{c,k}, xi),  fitline{c,k} );
        plot(fit_valid{c,k}(:,1), fit_valid{c,k}(:,2), symbo{c,k},  xi, cum_gaussfit(psy_perf{c,k}, xi),  fitline{c,k} );
        %xlabel('Rotation Angles');   
        xlabel('Rotation Vel(deg/s)');   
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        hold on;
        legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
        legend_txt{k*2} = [''];
    end
end

xx = [fit_valid{c,k}(:,1), fit_valid{c,k}(:,2)]
yy =[ xi', cum_gaussfit(psy_perf{c,k}, xi)']
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
text(10,8, 'u                   sigma         2targ-corr rate        corr-rate      2targ-trial%      R square');

text(10,6,num2str(Bias_psy_dec) );
text(20,6,num2str(Thresh_psy_dec) );

for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        text(0,8-k-(c-1)*3, num2str(unique_stim_type(k)));  % non-microstim
        text(10,8-k-(c-1)*3,num2str(Bias_psy{c,k}) );
        text(20,8-k-(c-1)*3,num2str(Thresh_psy{c,k}) );
        text(30,8-k-(c-1)*3,num2str(correct_proportion(k)) );
        text(40,8-k-(c-1)*3,num2str(correct_proportion_overall(k)) );
        text(50,8-k-(c-1)*3,num2str(twoTarg_proportion(k)) );
        text(60,8-k-(c-1)*3,num2str(m_Rsquar{c,k}));
    end
end

axis off;

dt = [  Bias_psy{1,1}   Thresh_psy{1,1}  m_Rsquar{1,1}   Bias_psy_dec  Thresh_psy_dec  repetition   1];
% % plot psycho over time


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flnm = [FILE(1:end-4)  '_'  num2str(order)];

addpath('Z:\Users\sheng\program\newtool');
writeexceldt(['Z:\Users\Courtney\yawtime\behavior11' ], flnm, dt);

saveas(gcf,  ['Z:\Users\Courtney\yawtime\'  flnm],  'png');
close(gcf);


%---------------------------------------------------------------------------------------
return;

%convert Rotation amp to vel. 08/15/2013
function [velmax maxAcc1s]=calpeakVel(dur, sig, mag)
%%% rot: dur=1,sig=4,
%%% trans: dur=1,sig=4,

vel = GenGaussian(dur, sig);
dist = cumtrapz(vel)/60;
for i = 1 : length(mag)
    dist = dist/max(dist)*mag(i);
    vel = diff(dist)*60;    
    accel= diff(vel)*60;   
end
velmax = abs(max(vel)) ;       %%
maxAcc1s = abs(max(accel))*100;%% cm/s/s

function [data] = GenGaussian(dur, sig)
t0 = dur/2;
t = 0:1/60:dur;

% Generate the Gaussian.
data = exp(-sqrt(2)* ((t-t0) / (dur/sig)).^2);
%end 08/15/2013

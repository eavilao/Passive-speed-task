%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingRaw_cum(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(SpikeChan,:);   % spike rasters

temp_ROT_VELOCITY = data.moog_params(ROT_VELOCITY,:,MOOG);
temp_Duration = data.moog_params(DURATION,:,MOOG);


unique_duration = munique(temp_Duration');
%convert Rotation amp to vel. 08/15/2013
for i=1:length(temp_heading)
    %
    if  unique_duration ==  1000
        [vel Acc]=CalVel(1, 4, abs(temp_heading(i)), 0.5);
         velmax = vel; 
    else
        velmax = 2*abs(temp_heading(i))/(0.7+0.5);  % Trapzoid profile. Jing 09/05/2013
    end

     
    if(temp_heading(i)>=0)
        temp_heading(i)=velmax;
    else
        temp_heading(i)=-velmax;
    end
end
%end 08/15/2013
%%%%%%%%%%revised by xjy
% for i=1:length(temp_heading)
%    [velmax maxAcc1s]=calpeakVel(1, 4,abs(temp_heading(i)));
%    
%   if(temp_heading(i)>=0)
%       temp_heading(i)=velmax;
%   else
%       temp_heading(i)=-velmax;
%   end
% end
      

%%%%%%%%%%5

%%%%the following code is to count CV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spkdataforcv= squeeze(data.spike_data(SpikeChan, StartEventBin(1)-300:StartEventBin(1),:));
ind = 1;
for trialnum=1:length(StartEventBin)
    currtspk=spkdataforcv(:,trialnum);
    currtspktime=find( currtspk==1);
    curretisi=diff(currtspktime);
            curtmisi=mean(curretisi);
           curtsdisi=std(curretisi);  
           if curtmisi == 0
               cv_false = 1;
               continue;
           end
           currentcv(ind)=curtsdisi/curtmisi ;
           ind = ind + 1;
           clear currtspk  currtspktime curretisi  curtmisi curtsdisi
    
end
CV=nanmean( currentcv);
clear  currentcv

%%%


temp_targlum_1 = data.targ_params(TARG_LUM_MULT,:,T1);  %Jing 01/03/2013
temp_targlum_2 = data.targ_params(TARG_LUM_MULT,:,T2);

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
spike_rates = temp_spike_rates( select_trials);
total_trials = temp_total_trials( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
disc_heading = unique_heading( floor(length(unique_heading)/2)+1 : end );

StartEventBin = StartEventBin +115;
StopEventBin  = StopEventBin +115;

 StartOffset=200;
 StopOffset=-200;


   for jj = 1 : length(StartEventBin)
        start_spikebin = StartEventBin(1) + StartOffset;
        stop_spikebin  = StopEventBin(1) + StopOffset;
        spike_rates(jj) = mean(double(data.spike_data(SpikeChan, start_spikebin:stop_spikebin,jj)))*1000;
   end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part deals with for whatever reason, some trials are absolutely
% wrong, like >1000 spikes/s, makes no sense. tick this bad trial out,
% replace with one of the other data within same condition, this is similar
% to bootstrap. But in order to keep every running the same, use the mean
% of the rest instead. this happen only very rarely!!!
bad_tri = find(temp_spike_rates > 1000 );
if length(bad_tri) > 0
    for i=1:length(bad_tri)
		same_tri = find( (stim_type == stim_type(bad_tri(i))) & (heading == heading(bad_tri(i))) & (motion_coherence == motion_coherence(bad_tri(i))) ); % find other trials within same condition
		rest_tri = same_tri(same_tri ~= bad_tri(i));
		spike_rates(bad_tri(i)) = mean( spike_rates(rest_tri) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%neurometric dataset and calculate ROC, Choice Probability(CP)
%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
LEFT = 1;
RIGHT = 2;
for i= 1 : length(spike_rates) 
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
%     if (sum(events == IN_T1_WIN_CD) > 0)
%         choice(i) = RIGHT;
%     elseif (sum(events == IN_T2_WIN_CD) > 0)
%         choice(i) = LEFT;
%     else
%         disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
%     end
    
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
%choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for this trial
% psychometric dataset
psycho_correct = [];
fit_data_psycho = [];
N_obs = [];
for k = 1:length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)    
        if unique_stim_type(k) ==1
            c=1;
        end
         for i = 1:length(unique_heading)            
             trials_p =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) ) ;
             % make 'S' curve by using the rightward choice for y-axis
             correct_trials = (trials_p & (choice == RIGHT) );
             psycho_correct{c}(k,i) = 1*sum(correct_trials) / sum(trials_p); 
             fit_data_psycho_cum{c,k}(i, 1) = unique_heading( i );  
             fit_data_psycho_cum{c,k}(i, 2) = psycho_correct{c}(k,i);
             fit_data_psycho_cum{c,k}(i, 3) = sum(trials_p); 
         end
    end
end

% this part needs work later, for vestibular condition, coherence does not duplicate
one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence); 
repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition
resp_heading = [];
% Z_Spikes = spike_rates;
% % z-score data for later cp analysis across headings
% for k = 1:length(unique_stim_type)
%     for c = 1:length(unique_motion_coherence)
%         if unique_stim_type(k) ==1
%             c=1;
%         end
%         for i = 1:length(unique_heading)
%             select = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) ) ;  
%             z_dist = spike_rates(select);
%             z_dist = (z_dist - mean(z_dist))/std(z_dist);
%             Z_Spikes(select) = z_dist;
%         end
%     end
% end


Z_Spikes = spike_rates;
% z-score data for later cp analysis across headings
for k = 1:length(unique_stim_type)
    for i = 1:length(unique_heading)
       select1 = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & choice == LEFT   );
       select2 = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & choice == RIGHT   );
       if sum(select1) ==  0  ||  sum(select2) ==  0 
           select = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) );
           z_dist = spike_rates(select);
           z_dist = (z_dist - mean(z_dist))/std(z_dist);
           Z_Spikes(select) = z_dist;
       else
         select = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) );
         resp_left = spike_rates(select1);
         resp_right = spike_rates(select2);
         nn = length(resp_left) + length(resp_right);
         f1 = length(resp_left)/nn;
         f2 = length(resp_right)/nn;
         mu1 = mean(resp_left);
         var1 = var(resp_left);
         mu2 = mean(resp_right);
         var2 = var(resp_right);
         mu = f1*mu1 + f2*mu2;
         vv = f1*var1 + f2*var2 +f1*f2*(mu1-mu2)^2;
         z_dist = spike_rates(select);
         z_dist = (z_dist - mu)/sqrt(vv);
         Z_Spikes(select) =  z_dist;
       end        
    end
end
Z_Spikes_Ori = Z_Spikes; % keep a Z_Spikes unchanged for later use


% now group neuronal data into two groups according to monkey's choice
for k = 1:length(unique_stim_type)    % notice that the stim_type is double than disc_heading 
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) ) ;  
            resp{c,k,i} = spike_rates(select);   

            resp_mat{c,k}(i) = mean(resp{c,k,i});  % the mean firing rate for each heading 
            resp_mat_std{c,k}(i)= std(resp{c,k,i});
            resp_mat_err{c,k}(i) = std(resp{c,k,i}) / sqrt(repetition);
            % calculate CP, group data based on monkey's choice 
            resp_left_choose{c,k,i} = spike_rates(select & (choice == LEFT) );
            resp_right_choose{c,k,i} = spike_rates(select & (choice == RIGHT) );
            if (length(resp_left_choose{c,k,i}) <= 3) | (length(resp_right_choose{c,k,i}) <= 3)   % make sure each stim_type has at least 3 data values
          %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)  
                Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
            end 
        end  
        % now across all data 
        resp_left_all{c,k} = Z_Spikes( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_right_all{c,k} = Z_Spikes( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
        resp_all{c,k} = Z_Spikes( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) ); 
    end
end

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %----------- plot neuronal response aross time, see whether it is stable-------------------------------------------------------------------
% figure(9);
% kk{1}='b.-'; kk{2}='r.-'; kk{3}='k.-'; kk{4}='g.-'; kk{5}='y.'; kk{6}='y.-'; 
% nu{1}='b.--'; nu{2}='r.--'; nu{3}='k.--'; nu{4}='g.--'; nu{5}='y.--'; nu{6}='y.--'; 
% plot(resp{1,2,1},kk{1});
% hold on;
% plot(resp{1,2,5},kk{4});
% hold on;
% plot(resp{1,2,9}, kk{2});
% hold on;
% ylim([0 100]);
% set(gca, 'ytick',[0:20:100]);
% set(gca, 'xtick',[0 repetition]);

for k=1:length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        % decide whether ves and vis is congruent tuning. Fit line by linear
        % regression first and compare the sign of each stim_type to decide whether
        % congruent or opposite, this is used to check whether congruent cells lead
        % to better neuronal performance in combined stim_type, and vice versa
        [rr,pp] = corrcoef(unique_heading, resp_mat{c,k}(:));
        line_re{c,k} = rr(1,2);
        line_p{c,k} = pp(1,2);
    end
end
   
hold off;
% %------------------------------------------------------------------------

% now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
% neurothreshold 
fit_data_neuro = [];
fit_data_neuro_cut = [];
for k = 1 : length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        for i = 1 : length(unique_heading)   % subtract the 0 heading               
            trials_n =logical( (motion_coherence == unique_motion_coherence(c)) & (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
            fit_data_neuro_cum{c,k}(i,3) = sum(trials_n);  % for later function fit use
            Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i+1},resp{c,k,i},100 );
           
     %        if  resp_mat{k}(1) < resp_mat{k}(end)
             if line_re{c,k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
                 % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
                 % here we asume if left response larger than right response then asign the left to be preferred direction
                 Neuro_correct{c,k}(i) = 1 - Neuro_correct{c,k}(i);            
             end  
        end
%          for i = 1 : length(unique_heading)-2   % for symetric headings/no 0 heading             
%             trials_n =logical( (motion_coherence == unique_motion_coherence(c)) & (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
%             fit_data_neuro_cum{c,k}(i,3) = sum(trials_n);  % for later function fit use
%             resp_comparison = [resp{c,k,length(unique_heading)/2},resp{c,k,1+length(unique_heading)/2} ];
%             if i < length(unique_heading)/2
%                  % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
%                  Neuro_correct{c,k}(i) =  rocN( resp_comparison,resp{c,k,i},100 ); % compare to the 0 heading stim_type, which is straght ahead
%             % use anti-neuron model instead, so compare each plus and minus headings 
%            %      Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i+1},resp{c,k,i},100 );
%             else
%                  Neuro_correct{c,k}(i) =  rocN( resp_comparison, resp{c,k,(i+2)},100 ); % compare to the 0 heading stim_type, which is straght ahead
%         %         Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i}, resp{c,k,(i+1)},100 );
%             end
%      %        if  resp_mat{k}(1) < resp_mat{k}(end)
%              if line_re{c,k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
%                  % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
%                  % here we asume if left response larger than right response then asign the left to be preferred direction
%                  Neuro_correct{c,k}(i) = 1 - Neuro_correct{c,k}(i);            
%              end  
%         end
    end
end

%choice probability
for k = 1 : length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
           c=1;
        end
        for i = 1 : length(unique_heading)  
            if (length(resp_left_choose{c,k,i}) > 3) & (length(resp_right_choose{c,k,i}) > 3)
               CP{c,k}(i) = rocN( resp_left_choose{c,k,i},resp_right_choose{c,k,i},100 );
            else
                CP{c,k}(i) = NaN;
            end
            if (length(resp_left_all{c,k}) > 3) & (length(resp_right_all{c,k}) > 3)
                CP_all{c,k} = rocN( resp_left_all{c,k},resp_right_all{c,k},100 );
            else
                CP_all{c,k} = NaN; 
            end
            if  line_re{c,k} > 0
                CP{c,k}(i) = 1 - CP{c,k}(i);
                CP_all{c,k} = 1 - CP_all{c,k};
            end
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% use Wichman's MLE method to estimate threshold and bias


for k = 1:length(unique_stim_type)
     for c = 1:length(unique_motion_coherence)
       wichman_psy = pfit(fit_data_psycho_cum{c,k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');  
        Thresh_psy{c,k} = wichman_psy.params.est(2);
        Bias_psy{c,k} = wichman_psy.params.est(1);
        psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
         fit_data_neuro_cum{c,k}(:,1) = unique_heading(unique_heading~=0);
%        fit_data_neuro_cum{c,k}(:,1) = unique_heading( unique_heading~=unique_heading(length(unique_heading)/2) & unique_heading~=unique_heading(1+length(unique_heading)/2) );
         fit_data_neuro_cum{c,k}(:,2) = Neuro_correct{c,k}(:);
       wichman_neu = pfit(fit_data_neuro_cum{c,k}(:,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');
       Thresh_neu{c,k} = wichman_neu.params.est(2);
%         [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_cum{c,k});
%         if ttt>300
%             ttt=300;
%         end
%         Thresh_neu{c,k} = ttt;
        % negative and positive infinite value means flat tuning
        if Thresh_neu{c,k}<0 | Thresh_neu{c,k}> 300
            Thresh_neu{c,k} = 300;
            wichman_neu.params.est(2) = 300;
        end
        Bias_neu{c,k} = wichman_neu.params.est(1);
        neu_perf{c,k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%do permutation to test the significance of CP_all{k}, re-calculate CP 2000 times
perm_num = 1000;
Z_Spikes_perm = Z_Spikes;
bin = 0.005;
x_bin = 0 : bin : 1;
for k = 1:length(unique_stim_type) 
    for c = 1:length(unique_motion_coherence)
       
        if unique_stim_type(k) ==1
           c=1;
        end
        for n = 1 : perm_num
            % temperarilly only use near-threshold heading angles where monkey make a guess mainly
            select = logical( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
            Z_Spikes_con{c,k} = Z_Spikes_perm( select );
            Z_Spikes_con{c,k} = Z_Spikes_con{c,k}(randperm(length(Z_Spikes_con{c,k})));   % permute spike_rates
            Z_Spikes_perm(select) = Z_Spikes_con{c,k};    % now in spike_rates, the corresponding data were permuted already

            resp_left_all_perm{c,k} = Z_Spikes_perm( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
            resp_right_all_perm{c,k} = Z_Spikes_perm( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 

            if  (length(resp_left_all{c,k}) > 3) & (length(resp_right_all{c,k}) > 3) 
                CP_all_perm{c,k}(n) = rocN( resp_left_all_perm{c,k}, resp_right_all_perm{c,k},100 );
            else
                CP_all_perm{c,k}(n) = NaN; 
            end
            if  line_re{c,k} > 0  
                CP_all_perm{c,k}(n) = 1 - CP_all_perm{c,k}(n);             
            end  

            resp_left_choose_perm = Z_Spikes_perm((heading == unique_heading(5)) & (stim_type == unique_stim_type(1)) & (motion_coherence == unique_motion_coherence(1)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
            resp_right_choose_perm = Z_Spikes_perm((heading == unique_heading(5)) & (stim_type == unique_stim_type(1)) & (motion_coherence == unique_motion_coherence(1)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 

            if  (length(resp_left_choose{1,1,5}) > 3) & (length(resp_right_choose{1,1,5}) > 3)
                CP_perm(n) = rocN( resp_left_choose_perm, resp_right_choose_perm,100 );
            else
                CP_perm(n) = NaN; 
            end
            if  line_re{1} > 0  
                CP_perm(n) = 1 - CP_perm(n);             
            end  

        end
        % now calculate p value or significant test
        if (length(resp_left_all{c,k}) > 3) & (length(resp_right_all{c,k}) > 3) 
            hist_perm = [];
            hist_perm = hist( CP_all_perm{c,k}(:), x_bin );  % for permutation
            bin_sum = 0;
            n = 0;
            while ( n < (CP_all{c,k}/bin) )
                 n = n+1;
                 bin_sum = bin_sum + hist_perm(n);
                 if CP_all{c,k} > 0.5                  % note it's two tail test
                    p{c,k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
                 else
                    p{c,k} = 2* bin_sum / perm_num;
                 end
            end
        else
            p{c,k} = NaN;
        end 

        % calculate p value for CP during straight ahead motion
        if (length(resp_left_choose{1,round(length(unique_heading)/2)}) > 3) & (length(resp_right_choose{1,1,round(length(unique_heading)/2)}) > 3)  
            hist_perm = [];
            hist_perm = hist( CP_perm(:), x_bin );  % for permutation
            bin_sum = 0;
            n = 0;
            while ( n < (CP{1}(round(length(unique_heading)/2))/bin) )
                 n = n+1;
                 bin_sum = bin_sum + hist_perm(n);
                 if CP{c,k}(round(length(unique_heading)/2)) > 0.5    % note it's two tail test
                    pp = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
                 else
                    pp = 2* bin_sum / perm_num;
                 end
            end
        else
            pp = NaN;
        end  
        p_a(c,k) = pp;
    end
end

% % use self-written maximum liklihood method, this should run faster than
% above
% for c = 1:length(unique_motion_coherence)
    % for k = 1:length(unique_stim_type)    
    %     [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum{c,k});
    %     Bias_psy{c,k} = bb;
    %     Thresh_psy{c,k} = tt;
    %     psy_perf{c,k}=[bb, tt];
    %     [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_cum{c,k});
    %     Bias_neu{c,k} = bbb;
    %     Thresh_neu{c,k} = ttt;
    %     neu_perf{c,k}=[bb, tt];
    % end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot psychometric and neurometric function here
h{1} = 'bo';  f{1} = 'b-';  g{1} = 'bo-';
h{2} = 'rd';  f{2} = 'r-';  g{2} = 'rd-';
h{3} = 'gs';  f{3} = 'g-';  g{3} = 'gs-';
figure(10);
set(10,'Position', [5,25, 980,650], 'Name', 'psycho_neurometic function');
orient landscape;
%plot psychometric function
axes('position',[0.05 0.3, 0.26 0.45]);
title('psychometric');
legend_txt = [];
for k = 1:length(unique_stim_type)
    xi = min(unique_heading) : 0.1 : max(unique_heading);   
    beta = [0, 1.0];
 %   plot data in logarithmical space instead of linspace
    plot(unique_heading, psycho_correct{1}(k,:), h{k}, xi, cum_gaussfit(psy_perf{1,k}, xi),  f{k} );
    xlabel('Rotation Peak vel(deg/s)');   
    ylim([0,1]);
    ylabel('Rightward Choices');
    hold on;
    legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
    legend_txt{k*2} = [''];
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot neurometric function
axes('position',[0.36 0.3, 0.26 0.45]);
title('neuroometric');
for k = 1:length(unique_stim_type)
    neu_heading = unique_heading(unique_heading~=0);
    xi = min(unique_heading) : 0.1 : max(unique_heading); 
%    xi = unique_heading(2):0.1:unique_heading(8);
    plot(neu_heading, Neuro_correct{1,k}, h{k},  xi, cum_gaussfit(neu_perf{1,k}, xi),  f{k} );
   % plot(neu_heading(2:7), Neuro_correct{k}(2:7), h{k},  xi, cum_gaussfit(neu_perf{k}, xi),  f{k} );
    xlabel('Rotation Peak vel(deg/s)');   
    ylim([0,1]);
    hold on;
%    betafit_ne_cutt(k)=betafit_ne_cut{k}(2);
%    also fit data with weibull function
%     [neuro_alpha(k) neuro_beta(k)]= weibull_fit(fit_data_neuro{k});
%     [neuro_alpha_cut(k) neuro_beta_cut(k)]= weibull_fit(fit_data_neuro_cut{k});  % tick the most outside datapoint out
end

%%%%%%  neurological raw data based on firing rate instead of ROC
axes('position',[0.7 0.3, 0.26 0.45]);
title('firing rate');
for k = 1:length(unique_stim_type)
    errorbar(unique_heading, resp_mat{1,k}(:), resp_mat_err{1,k}(:),g{k} );
    xlabel('Rotation Peak vel(deg/s)');
    ylabel('Firing rate(spikes/s)');   
    xlim([min(unique_heading),max(unique_heading)]);
    hold on;
end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output some text of basic parameters in the figure
axes('position',[0.05,0.8, 0.9,0.15] );
xlim( [0,100] );
ylim( [2,10] );
text(0, 10, FILE);
text(20,10,'coherence =');
text(30,10,num2str(unique_motion_coherence) );
text(45,10,'repetition =');
text(55,10,num2str(repetition) ); 
text(5,8, 'Psy: u      threshold                   Neu:u         threshold               CP           p');
text(0,8, 'stim');
for k = 1:length(unique_stim_type)
    text(0,8-k, num2str(unique_stim_type(k)));
    text(5,8-k,num2str(Bias_psy{c,k} ));
    text(12,8-k,num2str(Thresh_psy{c,k} ));

    text(27,8-k,num2str(Bias_neu{1,k} ));
    text(35,8-k,num2str(Thresh_neu{1,k} ));

     text(47,8-k,num2str(CP_all{1,k}') ); 
     text(55,8-k,num2str(p{1,k}') );   

end
axis off;


xx = [unique_heading, psycho_correct{1}(k,:)']
yy =[ xi', cum_gaussfit(psy_perf{1,k}, xi)']

aa = [neu_heading, Neuro_correct{1,k}']
bb =[ xi', cum_gaussfit(neu_perf{1,k}, xi)']

xlswrite(['Z:\Users\Courtney\figure\ ' FILE(1:end-4)],  xx);
xlswrite(['Z:\Users\Courtney\figure\ ' FILE(1:end-4)],  yy,  ['C1:D' num2str(size(yy,1)) ] );
xlswrite(['Z:\Users\Courtney\figure\ ' FILE(1:end-4)],  aa,  ['E1:F' num2str(size(aa,1)) ] );
xlswrite(['Z:\Users\Courtney\figure\ ' FILE(1:end-4)],  bb,  ['G1:H' num2str(size(bb,1)) ] );

% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
%buff = sprintf(sprint_txt, FILE, unique_motion_coherence, repetition, Bias_psy{:}, Thresh_psy{:}, Bias_neu{:}, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:}, line_p{:}, line2_re, line2_p );
if length(unique_stim_type)==1 
   buff = sprintf(sprint_txt, FILE, 2, Thresh_neu{length(unique_motion_coherence),1} );   % visual 100% coherence
elseif length(unique_stim_type)>1
    buff = sprintf(sprint_txt, FILE, 1,CP{1,1}(5),CP{1,2}(5),CP{1,3}(5) );   % vestibular
end
%buff = sprintf(sprint_txt, resp_mat{1}, resp_mat_std{1}.^2, resp_mat{2}, resp_mat_std{2}.^2, resp_mat{3}, resp_mat_std{3}.^2  );
%buff = sprintf(sprint_txt, FILE, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:} );
%buff = sprintf(sprint_txt, FILE, Thresh_neu{:},line_re{1}, line_re{2}, line_p{1}, line_p{2} );
%buff = sprintf(sprint_txt, FILE, resp_mat{1}(:),resp_mat{2}(:),resp_mat{3}(:) );
%buff = sprintf(sprint_txt, FILE, fit_data_psycho_cum{1}(:, 2),fit_data_psycho_cum{2}(:, 2),fit_data_psycho_cum{3}(:, 2),Neuro_correct{1}(:),Neuro_correct{2}(:),Neuro_correct{3}(:) );
outfile = ['Z:\Users\Yong\CP.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         Coher\t repet\t Ve_P_u\t Vi_P_u\t co_P_u\t Ve_P_thr\t vi_P_th\t Co_P_thr\t Ve_N_u\t Vi_N_u\t co_N_u\t Ve_N_thr\t vi_N_th\t Co_N_thr\t Ves_CP\t Vis_CP\t Com_CP\t Ves_p\t Vis_p\t Com_p\t sign\t vesMax\t visMax\t coMax\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%--------------------------------------------------------------------------
function [velmax maxAcc1s]=calpeakVel(dur, sig, mag)


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
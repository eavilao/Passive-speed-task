%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cum(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
spike_rates = temp_spike_rates(select_trials);
total_trials = temp_total_trials( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
disc_heading = unique_heading(floor(length(unique_heading)/2)+1 : end );

fname = sprintf('%s_%d', FILE(1:end-4), SpikeChan); 
save (fname, 'trials', 'stim_type', 'heading', 'amplitude', 'motion_coherence', 'spike_rates', 'unique_stim_type', 'unique_heading'); 
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
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
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
% Z_Spikes_Ori = Z_Spikes; % keep a Z_Spikes unchanged for later use
%% balanced CP
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


%%

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
        for i = 1 : length(unique_heading)-1   % subtract the 0 heading               
            trials_n =logical( (motion_coherence == unique_motion_coherence(c)) & (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
            fit_data_neuro_cum{c,k}(i,3) = sum(trials_n);  % for later function fit use
            if i < (1+length(unique_heading))/2
                 % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
                 Neuro_correct{c,k}(i) =  rocN( resp{c,k,(1+length(unique_heading))/2},resp{c,k,i},100 ); % compare to the 0 heading stim_type, which is straght ahead
            % use anti-neuron model instead, so compare each plus and minus headings 
       %          Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i+1},resp{c,k,i},100 );
            else
                Neuro_correct{c,k}(i) =  rocN( resp{c,k,(1+length(unique_heading))/2}, resp{c,k,(i+1)},100 ); % compare to the 0 heading stim_type, which is straght ahead
      %           Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i}, resp{c,k,(i+1)},100 );
            end
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
        wichman_psy = pfit(fit_data_psycho_cum{c,k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        Thresh_psy{c,k} = wichman_psy.params.est(2);
        Bias_psy{c,k} = wichman_psy.params.est(1);
        psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
        fit_data_neuro_cum{c,k}(:,1) = unique_heading(unique_heading~=0);
        %        fit_data_neuro_cum{c,k}(:,1) = unique_heading( unique_heading~=unique_heading(length(unique_heading)/2) & unique_heading~=unique_heading(1+length(unique_heading)/2) );
        fit_data_neuro_cum{c,k}(:,2) = Neuro_correct{c,k}(:);
        wichman_neu = pfit(fit_data_neuro_cum{c,k}(:,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        Thresh_neu{c,k} = wichman_neu.params.est(2);
       
        
        %          [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_cum{c,k});
        %          if ttt>300
        %              ttt=300;
        %          end
        %          Thresh_neu{c,k} = ttt;
        % negative and positive infinite value means flat tuning
        
        if Thresh_neu{c,k}<0 | Thresh_neu{c,k}> 300
            Thresh_neu{c,k} = 300;
            wichman_neu.params.est(2) = 300;
        end
        
        % For normal CPs
        Bias_neu{c,k} = wichman_neu.params.est(1);
        neu_perf{c,k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
        
        % For abnormal Bias/Thr
        %           Bias_neu{c,k} = wichman_neu.params.est(1);
        %         neu_perf{c,k} = [bbb,ttt];
        
        
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
    xlabel('Heading Angles');   
    ylim([0,1]);
    ylabel('Rightward Choices');
    hold on;
    legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
    legend_txt{k*2} = [''];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psychometric, neurometric, CP over time
% run the slide threshold over time, see whether performance fluctuate across time
span = 5;  % calculate threshod every ? repeats;
slide = 1;  % slide threshod with increment of ? repeats;
BegTrial_shift = BegTrial;
EndTrial_shift = BegTrial_shift   + span*one_repetition-1;
n=0;
while EndTrial_shift <= EndTrial
    n = n + 1;
    select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
    stim_type_shift = temp_stim_type( select_trials_shift );
    heading_shift = temp_heading( select_trials_shift );
    unique_stim_type_shift = munique(stim_type_shift');
    unique_heading_shift = munique(heading_shift');
    total_trials_shift = temp_total_trials( select_trials_shift);
    spike_rates_shift = temp_spike_rates( select_trials_shift );
    Z_Spikes_shift = Z_Spikes_Ori;
    for k = 1:length(unique_stim_type)
        for i = 1:length(unique_heading)
             trials_shift =logical( (heading_shift == unique_heading(i)) & (stim_type_shift == unique_stim_type(k)) ) ;
             trials_shift2 = logical( (trials >= BegTrial_shift-BegTrial) & (trials <= EndTrial_shift-BegTrial) );
             trials_shift_CP = logical( (trials >= BegTrial) & (trials <= EndTrial_shift-BegTrial) );
             trials_shift3 = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
             correct_trials_shift = (trials_shift & (total_trials_shift == CORRECT) );
             % neural
             resp_heading_shift{k,i} = spike_rates_shift(trials_shift );
             % choice probability
             resp_left_choose_shift{k,i} = spike_rates(trials_shift3 & trials_shift_CP & (choice == LEFT) );
             resp_right_choose_shift{k,i} = spike_rates(trials_shift3 & trials_shift_CP & (choice == RIGHT) );
             if (length(resp_left_choose_shift{k,i}) <= span*0.25) | (length(resp_right_choose_shift{k,i}) <= span*0.25)
                 Z_Spikes_shift(trials_shift3 &trials_shift_CP) = 9999;
             end
             % make 'S' curve by using the rightward choice for y-axis
             if ( unique_heading(i) < 0 )
                 correct_rate_shift(i) = 1 - 1*sum(correct_trials_shift) / sum(trials_shift); 
             else
                 correct_rate_shift(i) = 1*sum(correct_trials_shift) / sum(trials_shift); 
             end         
         end
         fit_data_psycho_cum_shift{k}(:, 1) = fit_data_psycho_cum{k}(:, 1);  
         fit_data_psycho_cum_shift{k}(:, 2) = correct_rate_shift(:);
         fit_data_psycho_cum_shift{k}(:, 3) = span;
         [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{k});
         psy_thresh_shift(k,n) = tt;
         % for neuronal performence over time
         for i = 1 : length(unique_heading)-1   % subtract the 0 heading
             if i < (1+length(unique_heading))/2
                 Neuro_correct_shift{k}(i) =  rocN( resp_heading_shift{k,(1+length(unique_heading))/2},resp_heading_shift{k,i},100 ); % compare to the 0 heading stim_type, which is straght ahead
             else
                 Neuro_correct_shift{k}(i) =  rocN( resp_heading_shift{k,(1+length(unique_heading))/2},resp_heading_shift{k,i+1},100 ); % compare to the 0 heading stim_type, which is straght ahead
             end
             if  resp_mat{k}(1) < resp_mat{k}(end)  
                 Neuro_correct_shift{k}(i) = 1 - Neuro_correct_shift{k}(i);            
             end  
         end
         fit_data_neu_cum_shift{k}(:, 1) = unique_heading(unique_heading~=0);  
         fit_data_neu_cum_shift{k}(:, 2) = Neuro_correct_shift{k}(:);
         fit_data_neu_cum_shift{k}(:, 3) = span;
         [bbb,ttt] = cum_gaussfit_max1(fit_data_neu_cum_shift{k});
         neu_thresh_shift(k,n) = ttt; 
         % choice probability
         resp_left_all_shift{k} = Z_Spikes_shift( trials_shift_CP & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes_shift~=9999) ); 
         resp_right_all_shift{k} = Z_Spikes_shift( trials_shift_CP & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes_shift~=9999) );
         if (length(resp_left_all_shift{k}) <= span*0.25) | (length(resp_right_all_shift{k}) <= span*0.25)
             CP_all_shift{k}(n) = NaN;
         else
             CP_all_shift{k}(n) = rocN( resp_left_all_shift{k},resp_right_all_shift{k},100 );
         end
         if  resp_mat{k}(1) < resp_mat{k}(end)  
              CP_all_shift{k}(n) = 1 - CP_all_shift{k}(n);
         end   
    end   
    BegTrial_shift = BegTrial_shift + slide*one_repetition;
    EndTrial_shift = EndTrial_shift + slide*one_repetition; 
end
% plot psycho
axes('position',[0.05,0.05, 0.26,0.15] );
for k = 1:length(unique_stim_type)
    plot(psy_thresh_shift(k,:), f{k});
    hold on;
    xlabel('Repetition');  
    ylabel('Threshold');
    xlim([1, n]);
    ylim( [min(min(psy_thresh_shift(:,:))), max(max(psy_thresh_shift(:,:)))] );   
end
% plot neuro
axes('position',[0.36,0.05, 0.26,0.15] );
for k = 1:length(unique_stim_type)  
    plot(neu_thresh_shift(k,:), f{k});
    hold on;
    xlabel('Repetition');  
    ylabel('Threshold');
    xlim([1, n]);
%    ylim( [min(min(neu_thresh_shift(:,:))), 100] );   % just cut off at 100 deg, don't trust too big numbers
    ylim( [min(min(neu_thresh_shift(:,:))), max(max(neu_thresh_shift(:,:)))] );   % just cut off at 100 deg, don't trust too big numbers
end
% plot Choice Probability
axes('position',[0.7,0.05, 0.26,0.15] );
for k = 1:length(unique_stim_type)    
    plot(CP_all_shift{k}(:), f{k});
    hold on;
    xlabel('Repetition');  
    ylabel('CP');
    xlim([1, n]);
    ylim( [0, 1] );   
    grid on;
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
    xlabel('Heading Angles');   
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
    xlabel('Heading Angle (deg)');
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
text(65,10,'SpikeChan =' );
text(75,10,num2str(SpikeChan) );

text(5,8, 'Psy: u      threshold        err           Neu:u         threshold         err              CP        p');
text(0,8, 'stim');
for k = 1:length(unique_stim_type)
    text(0,8-k, num2str(unique_stim_type(k)));
    text(5,8-k,num2str(Bias_psy{1,k} ));
    text(12,8-k,num2str(Thresh_psy{1,k} ));

    text(30,8-k,num2str(Bias_neu{1,k} ));
    text(40,8-k,num2str(Thresh_neu{1,k} ));

     text(50,8-k,num2str(CP_all{1,k}') ); 
     text(60,8-k,num2str(p{1,k}') );   

end
axis off;
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
%buff = sprintf(sprint_txt, FILE, unique_motion_coherence, repetition, Bias_psy{:}, Thresh_psy{:}, Bias_neu{:}, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:}, line_p{:}, line2_re, line2_p );
buff = sprintf(sprint_txt, FILE, SpikeChan, unique_motion_coherence, repetition, Bias_psy{:}, Thresh_psy{:}, Bias_neu{:}, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:}, line_p{:});

% if length(unique_stim_type)==1 
%    buff = sprintf(sprint_txt, FILE, 2, Thresh_neu{length(unique_motion_coherence),1} );   % visual 100% coherence
% elseif length(unique_stim_type)>1
%     buff = sprintf(sprint_txt, FILE, 1,CP{1,1}(5),CP{1,2}(5),CP{1,3}(5) );   % vestibular
% end
%buff = sprintf(sprint_txt, resp_mat{1}, resp_mat_std{1}.^2, resp_mat{2}, resp_mat_std{2}.^2, resp_mat{3}, resp_mat_std{3}.^2  );
%buff = sprintf(sprint_txt, FILE, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:} );
%buff = sprintf(sprint_txt, FILE, Thresh_neu{:},line_re{1}, line_re{2}, line_p{1}, line_p{2} );
%buff = sprintf(sprint_txt, FILE, resp_mat{1}(:),resp_mat{2}(:),resp_mat{3}(:) );
%buff = sprintf(sprint_txt, FILE, fit_data_psycho_cum{1}(:, 2),fit_data_psycho_cum{2}(:, 2),fit_data_psycho_cum{3}(:, 2),Neuro_correct{1}(:),Neuro_correct{2}(:),Neuro_correct{3}(:) );
%outfile = ['Z:\Users\Yong\CP.dat'];
outfile = ['C:\Users\erico\Documents\MATLAB\HD_Analysis\All_HD_Analysis.dat'];  
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t SpikeChan\t         Coher\t repet\t Ve_P_bias\t Vi_P_bias\t co_P_bias\t Ve_P_thr\t vi_P_th\t Co_P_thr\t Ve_N_bias\t Vi_N_bias\t co_N_bias\t Ve_N_thr\t vi_N_th\t Co_N_thr\t Ves_CP\t Vis_CP\t Com_CP\t Ves_p\t Vis_p\t Com_p\t sign\t vesMax\t visMax\t coMax\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fname = sprintf('%s_%d', FILE(1:end-4), SpikeChan); 
%print (fname, '-dpsc2')
fclose(fid);
%saveas(gcf, [ 'c:\'  FILE(1:end-4) ],  'png');
saveas(gcf, [ 'C:\Users\erico\Documents\MATLAB\HD_Analysis\'  fname],  'png');

close(gcf);

%---------------------------------------------------------------------------------------

% figure(10);saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '.png'],'png')
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '.fig'])
% close (9);close (10);
% [StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
% HeadingDis_cum_PSTH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% 
% figure(2);orient landscape
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '_PSTH.png'],'png')
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '_PSTH.fig'])
% close(2);



return;
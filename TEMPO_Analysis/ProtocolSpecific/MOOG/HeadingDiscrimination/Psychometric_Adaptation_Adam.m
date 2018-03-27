%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for heading discrimination task
%--	07/16/04 GY, edited by AZ 2011
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_Adaptation_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
clear mex %need to clear mex functions otherwise the plexon crashes after a limited number of calls.
PLOT=0;
overtimeplot = 1;  % whether to calcuate performance over time
span = 5;  % for overtimeplot: calculate over 'span; repeats;
slide = 1;  % for overtimeplot: slide with increment of 'slide' repeats;
global monkey_no monkey prePST postPST preSACC postSACC postMVMT head_cut FR_del segval PCAcentered base_use
RESULTS_DIR = 'C:\Zaidel\Monkey\BCM\Visual-vestibular adaptation\results\';
TEMPO_Defs;
Path_Defs; %AZ: takes too much time (and I don't like that it places Z:\ stuff above my work stuff) - replaced by a dummy in my path
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
channels=mod(SpikeChan,100);
cells=0; cellID=[]; %initialize
sigma=20; half_guas=sigma*3.5;
dist=(-(half_guas):(half_guas))'; gaus=1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))*1000; %sigma in ms, 3.5 X sigma. Gaussian for smoothing the PSTH (to form the SDF)
use_tempo_spikes=0;

%AZ this is in order to get PRE data e.g. pre-thresholds in order to use a 'prior' for calculating POSTadapt thresh (if there is only one datapoint POST)
load(strcat(RESULTS_DIR,'monkey_XLSlist.mat'))
monk=find((str2num(FILE((strfind(FILE,'m')+1):(strfind(FILE,'c')-1))))==monkey_no); %I added 1000 to the monkey number if it had labyrinthectomy
list=cell2mat(eval(sprintf('%s_list',upper(monkey{monk}))));
id=find(all(cell2mat(arrayfun(@(x)(strcmp(x.htb_file,FILE(1:end-4))),list,'uniformoutput',false)'),2)); %AZ find which file in the list we are currently testing
precellID=[];
if list(id).block==3 || list(id).block==5 %AZ if this is the POST adaptation block
    PRIORid=find(arrayfun(@(x)(x.session==list(id).session & x.block==1),list)); %AZ find PRE block of the session
    load(strcat(RESULTS_DIR,'Psychometric_adapt_',FILE(strfind(FILE,'m'):strfind(FILE,'c')-1),'.mat'),list(PRIORid).htb_file) %AZ I assume that the PRE data is calculated/saved already
    try PREunique_stim_type=eval(sprintf('%s.unique_stim_type',list(PRIORid).htb_file)); end %incase we don't have the first block
    try PRIORthresh95CI=eval(sprintf('%s.Thresh95CI_psy',list(PRIORid).htb_file)); end %incase we don't have the first block
    %get pre PCA data
    %TBD for block 5 I want the shift from block3?
    try precellID=eval(sprintf('%s.cellID',list(PRIORid).htb_file)); end %incase we don't have the first block
    try prePCAcoeffVes=eval(sprintf('%s.PCAcoeffVes',list(PRIORid).htb_file)); end %try incase we don't have the first block
    try prePCAcoeffVis=eval(sprintf('%s.PCAcoeffVis',list(PRIORid).htb_file)); end %try incase we don't have the first block
    try prePCAsdf_meanVes=eval(sprintf('%s.PCAsdf_meanVes',list(PRIORid).htb_file)); end %try incase we don't have the first block
    try prePCAsdf_meanVis=eval(sprintf('%s.PCAsdf_meanVis',list(PRIORid).htb_file)); end %try incase we don't have the first block
    clear(list(PRIORid).htb_file)
end

%get the column of values for azimuth and elevation and stim_type
temp_azimuthM = data.moog_params(AZIMUTH,:,MOOG);
temp_azimuthC = data.moog_params(AZIMUTH,:,CAMERAS); %added by AZ (to deduce selection of visual/vest from tempo)
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_duration = data.moog_params(DURATION,:,MOOG); %added by AZ 2/21/11
if ~all(temp_duration==1000) %added by AZ
    error('AZ: Trial duration was not 1s! Exclude experiment!');
end
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);
temp_mask_status = data.moog_params(MASK_STATUS,:,MOOG);
temp_mask_radius = data.moog_params(MASK_RADIUS,:,MOOG);
temp_microstim = data.moog_params(MICROSTIM,:,MOOG);
temp_delta = data.moog_params(VESTIB_HEADING_OFFSET,:,MOOG); %added by AZ

trials = 1:length(temp_heading);		% a vector of trial indices.
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials);
amplitude= temp_amplitude( select_trials );
duration= temp_duration( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);
mask_status= temp_mask_status( select_trials );
mask_radius= temp_mask_radius( select_trials );
microstim = temp_microstim(select_trials );
delta = temp_delta(select_trials ); %added by AZ
azimuthM = temp_azimuthM(select_trials ); %added by AZ
azimuthC = temp_azimuthC(select_trials ); %added by AZ

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_duration = munique(duration');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');
unique_heading_nonzero = unique_heading(unique_heading~=0);
unique_delta = munique(delta'); %added by AZ
unique_azimuthM_comb = munique(azimuthM(stim_type==3)'); %added by AZ (only looking at combined trials)
unique_azimuthC_comb = munique(azimuthC(stim_type==3)'); %added by AZ (only looking at combined trials)
if sum(select_trials)/(length(unique_heading)*length(unique_stim_type))~=round(sum(select_trials)/(length(unique_heading)*length(unique_stim_type)))
    error(sprintf('AZ: trial count (%u) is not a multiple of unique_headingXunique_stim_type (choose trials 1:%u)',sum(select_trials),floor(sum(select_trials)/(length(unique_heading)*length(unique_stim_type)))*length(unique_heading)*length(unique_stim_type)))
end
%AZ: find meanDelta_OnOff (used the mean of the visual/vestibular stimuli as heading) and VestVisual (=0 rewarded by vestibular cue; =1 rewarded by visual cue)
%AZ: the test below should use '==', but there are sometimes some small (0.05) differences due to the way it is saved (only allows 5 digits o/w truncates), hence the '<' comparison.
if all(abs(sort(unique_azimuthM_comb)-sort(90-unique_heading-unique_delta))<0.1) & all(abs(sort(unique_azimuthC_comb)-sort(90-unique_heading))<0.1)
    meanDelta_OnOff=0; VestVisual=1;
elseif all(abs(sort(unique_azimuthM_comb)-sort(90-unique_heading))<0.1) & all(abs(sort(unique_azimuthC_comb)-sort(90-unique_heading+unique_delta))<0.1)
    meanDelta_OnOff=0; VestVisual=0;
elseif all(abs(sort(unique_azimuthM_comb)-sort(90-unique_heading-0.5*unique_delta))<0.1) & all(abs(sort(unique_azimuthC_comb)-sort(90-unique_heading+0.5*unique_delta))<0.1)
    meanDelta_OnOff=1; VestVisual=2; %fixation condition
else error('AZ: unable to determine whether rewarded by vestibular/visual cue.')
end

% 'one repetition' is different when values are added (i.e., extra stim types to increase ratio of combined to vestibular trials)  [** currently only works for 2 stim types **]  -CRF 12/2009
num_extras = 0;
if length(unique_stim_type) == 2 && sum(stim_type == unique_stim_type(1)) ~= sum(stim_type == unique_stim_type(2))
    if sum(stim_type == unique_stim_type(1)) > sum(stim_type == unique_stim_type(2))
        duplicated_stim = 1;
        other_stim = 2;
    else
        duplicated_stim = 2;
        other_stim = 1;
    end
    num_extras = round(sum(stim_type == unique_stim_type(duplicated_stim))/sum(stim_type == unique_stim_type(other_stim))) - 1;
end
one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence);
repetition = floor( length(heading)/one_repetition ); % take minimum repetition

if ~isempty(channels) %get plexon event data - added by AZ 4/12/11
    plex_PATH=strcat(PATH(1:end-4),'Plexon\sorted\');
    plex_FILE=strcat(plex_PATH,FILE(1:end-4),'-01.plx'); %sorted filename convention
    [Sync_n, Sync_ts, Sync_sv] = plx_event_ts(plex_FILE, 2); %sync pulse
    try
        [Event_n, Event_ts, Event_sv] = plx_event_ts(plex_FILE, 257); %plexon strobed channel is 257
        Succ_idx=find(Event_sv==12); %the trials that finished successfully
        Succ_idx=Succ_idx(select_trials); %of those that finished successfully, take only the trials selected in the XLS database
        for ss=1:length(Succ_idx) %find the stop and start timestamps that came just before successful trial completion
            Stop_diff=Succ_idx(ss)-find(Event_sv==StopCode); Stop_idx(ss)=Succ_idx(ss)-min(Stop_diff(Stop_diff>0)); %min distance of stop code from successful code
            Start_diff=Succ_idx(ss)-find(Event_sv==StartCode); Start_idx(ss)=Succ_idx(ss)-min(Start_diff(Start_diff>0)); %min distance of start code from successful code
            %currently I do not support saccade PSTHs for fixation experiments (VestVisual==2) since single cue trials have saccades and combined cue trials don't...
            if VestVisual<2, Sacc_diff=Succ_idx(ss)-find(Event_sv==7); Sacc_idx(ss)=Succ_idx(ss)-min(Sacc_diff(Sacc_diff>0)); end %min distance of Saccade code (7) from successful code
            STARTSync_diff=Sync_ts-Event_ts(Start_idx(ss)); STARTSync_shift(ss)=min(STARTSync_diff(STARTSync_diff>0))+1/60; %difference between start-code and SECOND sync-pulse (that is why I add 1/60 = one syncpulse) which signifies the actual stim start
            Start_ts(ss)=Event_ts(Start_idx(ss))+STARTSync_shift(ss); %NB: this does not represent the beginning of movement, rather stim on (since I can have a delay at the beginning).
            Stop_ts(ss)=Start_ts(ss)+round(10*(Event_ts(Stop_idx(ss))-Event_ts(Start_idx(ss))))/10; %duration of the trial (from start to stop, rounded off to nearest 100ms since I see that a 1sec stim is 1.015ms long in tempo!)
            if VestVisual<2, Sacc_del(ss)=Event_ts(Sacc_idx(ss))-Stop_ts(ss); end %delay between the end of movement and choice selection (saccade). NB I use syncpulse to correct stim start/stop but the saccade event is realtime (no syncpulse correction) since it is input to tempo
        end
        if all((Stop_ts-Start_ts)<1.02 & (Stop_ts-Start_ts)>0.98), vis_del=0; elseif all((Stop_ts-Start_ts)<1.32 & (Stop_ts-Start_ts)>1.28), vis_del=315; end %vis_delay is meant to be 300, but is ~315 in practice. 'vis_del' is only used for calculating the tempo delay (delay between code 04 and actual movement)
    catch %the strobed channel was not saved. maybe try to recover the data in plexon from the sync-pulse...? Would need to align the sync pulses? (not trivial and not done yet)...Meantime here I use the tempo spikes rather
        use_tempo_spikes=1;
        channel=1;
        step_idx=find(diff(Sync_ts)>1); %Deduce trial length. NB I skip the first trial (negligible for the purpose of deducing trial length)
        max_len=max((Sync_ts(step_idx(2:end))-Sync_ts(step_idx(1:end-1)+1)));
        if (max_len>1.27 & max_len<1.32), vis_del=315; elseif (max_len>0.98 & max_len<1.02), vis_del=0; end %vis_delay is meant to be 300, but is ~315 in practice. 'vis_del' is only used for calculating the tempo delay (delay between code 04 and actual movement)
    end
    tempo_sync=squeeze(data.spike_data(2,:,:)); %sync pulse as saved in tempo
    for i=1:size(tempo_sync,2), Tdel(i)=max(find(tempo_sync(:,i),2))-1000; end %find the time after 1000s that the trial actually began, and add it when extracting spikes
    tempo_del=round(mean(Tdel)); %average delay between tempo and system
    del=tempo_del+vis_del; %[ms]. Average tempo delay (unwanted delay between tempo go signal and actual movement) + intended visual-delay (delay between vis onset and movement). The latter manifests as a 315ms delay (intended 300ms)
    accelFB = mean(data.eye_data(6,[(1000/5+1 - prePST*200):(2000/5 + postMVMT*200)]+round(del/5),select_trials & (temp_stim_type==1 | temp_stim_type==3)),3); %added by AZ - calculate the mean accelerometer values (at 200Hz)
    %confirm that BG_delay was selected (if not then remove the delay added above)
    [Cmax,Imax] = max(accelFB); [Cmin,Imin] = min(accelFB);
    Tmid=mean([Imax Imin])*0.005-prePST; %middle of the gaussian profile
    if (Cmax-Cmin) < 0.1
        warning('AZ: Accelerometer not working')
    elseif Tmid<0.35 % really before 0.5 (the expected middle)
        warning('AZ: BG_delay was not selected. There was no delay before movement')
        vis_del=0; del=tempo_del+vis_del; %corrected from above
        if use_tempo_spikes==0, Stop_ts=Stop_ts-0.3; end %since the movement ended 300ms before the 05 code (b/c there was no delay added to the beginning)
        accelFB = mean(data.eye_data(6,[(1000/5+1 - prePST*200):(2000/5 + postMVMT*200)]+round(del/5),select_trials & (temp_stim_type==1 | temp_stim_type==3)),3); %added by AZ - calculate the mean accelerometer values (at 200Hz)
    end
end
if use_tempo_spikes
    temp_spike_data = data.spike_data(channel,:,select_trials);
    PST{1}=squeeze(data.spike_data(channel,[(1001 - half_guas - prePST*1000):(2000 + half_guas + postPST*1000)]+del,select_trials))'; %tempo starts saving 1s before code-04, and there is a 300ms delay after code04 before movement
    spike_rates{1}= sum(temp_spike_data(channel,[(1001+FR_del):2000]+del,:),2)./((unique_duration-FR_del)*1e-3); %tempo code-04 to code-05 (excluding the first FR_del)
    for r=segval
        eval(sprintf('spike_rates%u{1}= sum(temp_spike_data(channel,[1001:1200]+r-100+del,:),2)./((200)*1e-3);',r)) %tempo code-04 to code-05 (excluding the first FR_del)
    end
    base_rates250{1}= sum(temp_spike_data(channel,1000+(-249:0)+del,:),2)./(250*1e-3); %find the baseline FR (average of the 250ms before mvmt)
    base_rates500{1}= sum(temp_spike_data(channel,1000+(-499:0)+del,:),2)./(500*1e-3); %find the baseline FR (average of the 500ms before mvmt)
    base_rates1000{1}= sum(temp_spike_data(channel,1000+(-999:0)+del,:),2)./(1000*1e-3); %find the baseline FR (average of the 1s before mvmt)
    cells=1; cellID(1)=101;
    if VestVisual<2, Sacc_del=(find(squeeze(data.event_data(channel,:,select_trials))==7)-find(squeeze(data.event_data(channel,:,select_trials))==5))/1000; end %delay between end of movement and choice saccade (for PSTH around saccade)
else %get plexon spike data
    [channels IX]=sort(channels); %sort according to channel #
    if ~isempty(channels), [~, ~, waveformFreq] = plx_information(plex_FILE); end
    for ch=1:length(channels)
        channel=channels(ch);
        units{ch}=find([mod(floor(SpikeChan(IX(ch))/1e2),10) mod(floor(SpikeChan(IX(ch))/1e3),10) mod(floor(SpikeChan(IX(ch))/1e4),10) mod(floor(SpikeChan(IX(ch))/1e5),10)]); %unit numbers
        for u=1:length(units{ch})
            cells=cells+1;
            unit=units{ch}(u);
            cellID(cells)=channel+ 10^(unit+1);
            [n_temp p_temp ts_temp wav_temp]=plx_waves_v(plex_FILE, channel, unit); %spike waveform
            wav=wav_temp(ts_temp>min(Start_ts) & ts_temp<max(Stop_ts),:); %select waveforms within the trials selected
            waveform{cells}=wav(randperm(length(wav),100),:); %randomly selec 100 spikes
            [Spike_n{cells}, Spike_ts{cells}] = plx_ts(plex_FILE, channel, unit);
            PST{cells}(1:sum(select_trials),1:(half_guas + prePST*1000 + unique_duration + postPST*1000 + half_guas +1))=0; %INITIALIZE PST=peri-stimulus (spike) times.  (+1 for spike at 0s)
            for ss = 1:sum(select_trials) % ss marks the index of trial
                % Start_ts is when the visual stimulus appears (not when mvmnt begins). Since I keep the stim static before mvmt (Jing changed tempo) I define the PSTH epoch rather by Stop_ts.
                epoch_start=Stop_ts(ss) - unique_duration*1e-3 - (prePST + half_guas/1000); %time epoch to save for PSTH in seconds
                epoch_end=Stop_ts(ss) + postPST + half_guas/1000;
                IN_epoch=(epoch_start < Spike_ts{cells}) & (Spike_ts{cells} < epoch_end);
                epoch_spike_time=ceil((Spike_ts{cells}(IN_epoch)-epoch_start)*1000);
                PST{cells}(ss,epoch_spike_time)=1; %take all the spike timestamps that lie around the trial. NB PST has half_guas + prePST*1000 before the actual response
            end
            %calculate the spike rate for each interval as the sum of spikes WITHIN the PST interval from (after prePST+0.5*guass time : before postPST+0.5*guass time)
            spike_rates{cells}= sum(PST{cells}(:,prePST*1000+half_guas+(1+FR_del:unique_duration)),2)./((unique_duration-FR_del)*1e-3);
            for r=segval
                eval(sprintf('spike_rates%u{cells}= sum(PST{cells}(:,prePST*1000+half_guas+r-100+(1:200)),2)./((200)*1e-3);',r))
            end
            base_rates250{cells}= sum(PST{cells}(:,prePST*1000+half_guas+(-249:0)),2)./(250*1e-3); %find the baseline FR (average of the 250ms before mvmt)
            base_rates500{cells}= sum(PST{cells}(:,prePST*1000+half_guas+(-499:0)),2)./(500*1e-3); %find the baseline FR (average of the 500ms before mvmt)
            base_rates1000{cells}= sum(PST{cells}(:,prePST*1000+half_guas+(-999:0)),2)./(1000*1e-3); %find the baseline FR (average of the 1s before mvmt)
        end
    end
    if ~isempty(channels) & isempty(strfind(plex_FILE,'m26c130')) %get the LFP data (except for one recording session which had a bug: fn and ts are extremely long? adfreq wrong?)
        [n,samplecounts] = plx_adchan_samplecounts(plex_FILE);
        [n,names]=plx_adchan_names(plex_FILE);
        [n,adchans] = plx_ad_chanmap(plex_FILE);
        ch_names=str2num(names(samplecounts>0,3:4));
        LFPch{1} = adchans(ch_names(ch_names<=16)); %analog channel numbers for location1 LFPs (AD01..AD16 are ADchannel # 0-15)
        LFPdata{1}=zeros(sum(select_trials),length(LFPch{1}),(prePST*1000 + unique_duration + postPST*1000)/10); %100Hz
        LFPch{2} = adchans(ch_names(ch_names>=17 & ch_names<=32)); %analog channel numbers for location2 LFPs (AD17..AD32 are ADchannel # 16-31))
        LFPdata{2}=zeros(sum(select_trials),length(LFPch{2}),(prePST*1000 + unique_duration + postPST*1000)/10); %100Hz
        for loc=1:2 %two possible electrode locations
            for i=1:length(LFPch{loc})
                [adfreq, n, ts, fn, ad] = plx_ad_v(plex_FILE, LFPch{loc}(i));
                if length(ts)>2, error(strcat(plex_FILE,' paused multiple times - is the file corrupt?')); end
                ad_1kHz=downsample(ad(1:fn(1)),adfreq/1000); %1st segment (incase the plexon file was paused). Also, downsample the LFP to 1kHz
                for segment=2:length(ts) %add the additional segments if plexon was paused
                    ad_1kHz(end:round((ts(segment)-ts(1))*1000))=NaN; %NaN pad the time that the system was paused. subtract ts(1) since ad_1kHz starts at ts(1) and it is added later
                    ad_1kHz=[ad_1kHz; downsample(ad(fn(segment-1)+(1:fn(segment))),adfreq/1000)]; %next segment (incase the plexon file was paused). Also, downsample the LFP to 1kHz
                end
                for ss =   1:sum(select_trials) % ss marks the index of trial
                    % Start_ts is when the visual stimulus appears (not when mvmnt begins). Since I keep the stim static before mvmt (Jing changed tempo) I define the PSTH epoch rather by Stop_ts.
                    LFPdata{loc}(ss,i,:)=decimate(ad_1kHz(round(1000*(ts(1)+Stop_ts(ss) - unique_duration*1e-3 - prePST))+(0:(prePST*1000 + unique_duration + postPST*1000)-1)),10); %downsample to 100Hz
                end
            end
        end
    end
end

%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)
LEFT = 1;
RIGHT = 2;
fixation=0; %AZ: for fixation blocks, the monkey makes no choice (only for combined visual/vestibular stimuli). Rather he is rewarded for fixation alone.
for i= 1 : length(total_trials)
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
        choice(i) = 0;
        %disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
        fixation=1; %a completed trial w/o a choice means that the monkey only needed to fixate
    end
end

correct_rate = [];
Z_Spikes_new=ones(cells,length(total_trials))*9999;
Z_Spikes_old=ones(cells,length(total_trials))*9999;
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
%         if unique_stim_type(k)==3 %combined condition
%             fixation = ~sum(choice(stim_type==3)); %AZ: if monkey made no choice fixation=1. i.e. rewarded for fixation alone
%         end
        for i = 1:length(unique_heading)
            if unique_stim_type(k) == 1 % for vestibular condition, take all the data regardless of visual coherence
                trials_select(i,:) =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) ) ;
            else
                trials_select(i,:) =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (motion_coherence==unique_motion_coherence(c)) ) ;
            end
            rightward_trials = (trials_select(i,:) & (choice == RIGHT) );
            leftward_trials = (trials_select(i,:) & (choice == LEFT) );
            rightward_rate = 1*sum(rightward_trials) / sum(trials_select(i,:));
            fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);
            fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
            fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select(i,:));
            
            for cell=1:cells
                resp{c,k,cell}(i,:) = spike_rates{cell}(trials_select(i,:)); %the response is the spike rate per trial (for heading i, stim k)
                resp_choice{c,k,cell}(i,:)= [mean(spike_rates{cell}(leftward_trials)) mean(spike_rates{cell}(rightward_trials))]; %two columns - av FR for R/L choice
                for r=segval
                    eval(sprintf('resp%u{c,k,cell}(i,:) = spike_rates%u{cell}(trials_select(i,:));',r,r))
                    eval(sprintf('resp_choice%u{c,k,cell}(i,:) = [mean(spike_rates%u{cell}(leftward_trials)) mean(spike_rates%u{cell}(rightward_trials))];',r,r,r))
                end
                base250{c,k,cell}(i,:) = base_rates250{cell}(trials_select(i,:));
                base500{c,k,cell}(i,:) = base_rates500{cell}(trials_select(i,:));
                base1000{c,k,cell}(i,:) = base_rates1000{cell}(trials_select(i,:));
                resp_mat{c,k,cell}(i) = mean(resp{c,k,cell}(i,:));  % the mean firing rate for each heading
                resp_mat_std{c,k,cell}(i)= std(resp{c,k,cell}(i,:));
                resp_mat_err{c,k,cell}(i) = std(resp{c,k,cell}(i,:)) / sqrt(repetition);
                raster{c,k,cell}(i,:,:) = PST{cell}(trials_select(i,:),:); % raster matrix (headings,trials,spikes). NB PST has half_guas + prePST*1000 before the actual response
                if VestVisual<2, SACCdel{c,k,cell}(i,:)=Sacc_del(trials_select(i,:)); 
                else SACCdel{c,k,cell}(i,:)=NaN; end
                
                %CP - group data based on monkey's choice. Zscore for pooling (below) and calculate CP per heading.
                resp_left_choose{c,k,cell,i} = spike_rates{cell}(trials_select(i,:) & (choice == LEFT)); %the response is the spike rate per trial (for heading i, stim k)
                resp_right_choose{c,k,cell,i} = spike_rates{cell}(trials_select(i,:) & (choice == RIGHT)); %the response is the spike rate per trial (for heading i, stim k)
                if (length(resp_left_choose{c,k,cell,i}) <= 3) || (length(resp_right_choose{c,k,cell,i}) <= 3)   % make sure each stim_type has at least 3 data values
                    CP{c,k,cell}(i) = NaN;
                    Z_Spikes_old(cell,trials_select(i,:)) = 9999;   % similar to NaN, just make a mark
                    Z_Spikes_new(cell,trials_select(i,:)) = 9999;   % similar to NaN, just make a mark
                else
                    CP{c,k,cell}(i) = rocN( resp_left_choose{c,k,cell,i},resp_right_choose{c,k,cell,i},100 );
                    dist = spike_rates{cell}(trials_select(i,:)); %resp{c,k,cell}(i,:); same thing
                    Z_Spikes_old(cell,trials_select(i,:)) = (dist - mean(dist))/std(dist); %original version.
                    %new, based on Kang and Maunsell JNP 2012
                    mu=(mean(resp_left_choose{c,k,cell,i}) + mean(resp_right_choose{c,k,cell,i}))/2;
                    sd=sqrt((std(resp_left_choose{c,k,cell,i})^2 + std(resp_right_choose{c,k,cell,i})^2)/2 + ((mean(resp_left_choose{c,k,cell,i}) - mean(resp_right_choose{c,k,cell,i}))^2)/4);
                    Z_Spikes_new(cell,trials_select(i,:)) = (dist - mu)/sd;
                end  %flipped below (1-CP) according to corr
            end
            %TBD LFP is only aligned to mvmt - maybe align also to Saccade and plot around choice saccade (like spikes)
            if exist('LFPdata')==1
                for loc=1:length(LFPdata)
                    LFP_MVMT{loc,k}(i,:,:)=mean(LFPdata{loc}(trials_select(i,:),:,1:(prePST*1000 + unique_duration + postMVMT*1000)/10));
                    % LFP_SACC{loc,k}(i,:,:)=TBD needs to loop through the trials like sdf
                end
            end %mean LFP matrix (headings,channel,LFPdata) TBD: currently I save headings seperately. is this necessary?
        end
        for cell=1:cells %calculate a regression of the FRs vs. heading (only for headings <head_cut degrees)
            [~,line_coef{k,cell},corr_coeff{k,cell},p_value{k,cell}]= xy_regress(repmat(unique_heading(abs(unique_heading)<head_cut),[repetition 1]),reshape(resp{c,k,cell}(abs(unique_heading)<head_cut,:),[],1));
            %calculate the center (i.e. stright ahead), and other parameters, of sigmoid fit
            if p_value{k,cell}<0.05
                init_param=[1 1 0 0];
                out = fminsearch(@(params) acotfun(params,unique_heading(abs(unique_heading)<head_cut), resp_mat{c,k,cell}(abs(unique_heading)<head_cut)),init_param) ;
                sigm_center{k,cell}=-out(3)/out(2);
                %figure; hold on; plot(unique_heading(abs(unique_heading)<head_cut),resp_mat{c,k,cell}(abs(unique_heading)<head_cut)); plot([-20:20],out(1) .* mod(acot(out(2) .* [-20:20] + out(3))+pi,pi) + out(4))
            else sigm_center{k,cell}=NaN; end
            %CP, across all data. TBD - do I really want to remove Z_Spikes = 9999 (i.e. headings with <=3?
            %original version
            resp_left_all_old{c,k,cell} = Z_Spikes_old(cell, (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes_old(cell,:)~=9999) );
            resp_right_all_old{c,k,cell} = Z_Spikes_old(cell, (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes_old(cell,:)~=9999) );
            resp_all_old{c,k,cell} = Z_Spikes_old(cell, (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (Z_Spikes_old(cell,:)~=9999) );
            if (length(resp_left_all_old{c,k,cell}) > 3) & (length(resp_right_all_old{c,k,cell}) > 3)
                CP_all_old{c,k,cell} = rocN( resp_left_all_old{c,k,cell},resp_right_all_old{c,k,cell},100 );
            else CP_all_old{c,k,cell} = NaN; end
            %new version (only diff here in name, the changes are above)
            resp_left_all_new{c,k,cell} = Z_Spikes_new(cell, (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes_new(cell,:)~=9999) );
            resp_right_all_new{c,k,cell} = Z_Spikes_new(cell, (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes_new(cell,:)~=9999) );
            resp_all_new{c,k,cell} = Z_Spikes_new(cell, (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (Z_Spikes_new(cell,:)~=9999) );
            if (length(resp_left_all_new{c,k,cell}) > 3) & (length(resp_right_all_new{c,k,cell}) > 3)
                CP_all_new{c,k,cell} = rocN( resp_left_all_new{c,k,cell},resp_right_all_new{c,k,cell},100 );
            else CP_all_new{c,k,cell} = NaN; end
            if  corr_coeff{k,cell} > 0 %if the FR regression is +ve
                CP_all_old{c,k,cell} = 1 - CP_all_old{c,k,cell};
                CP_all_new{c,k,cell} = 1 - CP_all_new{c,k,cell};
                CP{c,k,cell} = 1 - CP{c,k,cell}; %individually for the headings               
            end
        end
        % the correct rate does not take coherence into account,temporarily 05/29/09
        trials_rightward = find( (heading > 0) & (choice==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
        trials_leftward  = find( (heading < 0) & (choice==LEFT) & (stim_type==unique_stim_type(k))  ) ;
        trials_non0 = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k)) ); %exclude 0 headings
        correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_non0);
        
        aa = find(fit_data_psycho_cum{c,k}(:,2)>-99); % b/c choice rate could be NaN if a heading condition is absent.
        fit_valid{c,k}(:,1) = fit_data_psycho_cum{c,k}(aa,1);
        fit_valid{c,k}(:,2) = fit_data_psycho_cum{c,k}(aa,2);
        fit_valid{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
        bb = abs(fit_valid{c,k}(:,1))<head_cut; % only use headings around straight ahead for psychometric curve
        fit_use{c,k}=fit_valid{c,k}(bb,:);
    end
end

% If combined trials are fixation-only, remove that stim type for all subsequent analyses/plots -- CRF 1-28-10
if unique_stim_type(end) == 3
    index = find(unique_stim_type==3);
    if sum(fit_valid{1,index}(:,2)) == 0
        unique_stim_type(end) = [];
    end
end

%%%%%% use Wichman's MLE method to estimate threshold and bias
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        valid_thresh{c,k}=sum(~(fit_use{c,k}(1:end,2)==0 | fit_use{c,k}(1:end,2)==1))>1; %AZ the threshold calculated from pfit is only valid if there are >1 data points that aren't 1 or 0
        if unique_stim_type(k)==3 | valid_thresh{c,k} | ~exist('PRIORthresh95CI') | ~any(PREunique_stim_type==unique_stim_type(k)) | any(isnan(PRIORthresh95CI{PREunique_stim_type==unique_stim_type(k)})), thresh_prior=''; else thresh_prior=sprintf('-cosine %f %f',PRIORthresh95CI{PREunique_stim_type==unique_stim_type(k)}(1),PRIORthresh95CI{PREunique_stim_type==unique_stim_type(k)}(4)); end %AZ if the fit will be based on one data point, use PREthresh as a 'prior'
        %next line differs from Psychometric_Adam
        if sum(fit_use{c,k}(1:end,2)~=0 & fit_use{c,k}(1:end,2)~=1)>1 || ~strcmp(thresh_prior,'') %if there are at least two numbers that are not zero OR one or there is a prior
            fit_use{c,k}(fit_use{c,k}(1:end,2)==0,2)=1e-32; %AZ pfit program doesn't like pefect 0s (performs badly)
            fit_use{c,k}(fit_use{c,k}(1:end,2)==1,2)=1-(1e-32); %AZ pfit program doesn't like perfect 1s
            [wichman_psy wichman_Full] = pfit(fit_use{c,k}(1:end,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'BETA_PRIOR',thresh_prior,'FIX_LAMBDA',1e-3,'LAMBDA_EQUALS_GAMMA',1,'sens',0,'compute_stats','false','verbose','false','CONF',[0.025 0.159 0.841 0.975]); %95 percent and 1SD equivalent
            Thresh_psy{c,k} = wichman_psy.params.est(2);
            Thresh95CI_psy{c,k} = prctile(wichman_Full.params.sim(:,2),[2.5 15.9 84.1 97.5]'); %AZ I calculate percentile to prevent getting NaNs (vs. wichman_psy.params.lims(:,2);)
            Bias_psy{c,k} = wichman_psy.params.est(1);
            Bias95CI_psy{c,k} = wichman_psy.params.lims(:,1);
            under=prctile(wichman_Full.params.sim(:,1),0.25); %remove the most extreme 0.5 percent incase the bootstrap gave far outliers (due to the limited # of datapoints) which would bias the SD
            over=prctile(wichman_Full.params.sim(:,1),99.75);
            BiasSD_psy{c,k}=std(wichman_Full.params.sim(wichman_Full.params.sim(:,1)>=under & wichman_Full.params.sim(:,1)<=over,1));
            psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
            gamma{c,k} = wichman_psy.params.est(3);
            lambda{c,k} = wichman_psy.params.est(4);
        else %calculate the bias only using a simpler function
            [bb,tt] = cum_gaussfit_max1(fit_use{c,k});
            Thresh_psy{c,k} = NaN; Thresh95CI_psy{c,k}=[NaN NaN NaN NaN]';
            Bias_psy{c,k} = bb; Bias95CI_psy{c,k}=[NaN NaN NaN NaN]';
            BiasSD_psy{c,k}=NaN;
            psy_perf{c,k} =[bb,tt];
            gamma{c,k}=0; lambda{c,k}=0;
        end
        if Thresh_psy{c,k} <0 %invalid psychometric function
            Thresh_psy{c,k} = NaN; Thresh95CI_psy{c,k}=[NaN NaN NaN NaN]';
            Bias_psy{c,k} = NaN; Bias95CI_psy{c,k}=[NaN NaN NaN NaN]';
            BiasSD_psy{c,k}=NaN;
            psy_perf{c,k} =[NaN,NaN];
            warning('AZ: negative threshold value. Psychometric invalid')
        end
        %AZ Goodness-of-fit
        pfitcurve_dir_only = cum_gaussfit(psy_perf{c,k} , fit_use{c,k}(:,1)); %the model values only at the actual directions used (used to calc. the R-square)
        SStot = sum((fit_use{c,k}(:,2)-mean(fit_use{c,k}(:,2))).^2);
        SSerr = sum((pfitcurve_dir_only(:) - fit_use{c,k}(:,2)).^2);
        Rsquare{c,k} = 1-SSerr/SStot;
    end
end
% this is the output, you can use it for plot of example cells
for k = 1:length(unique_stim_type)
    step=(max(fit_use{c,k}(:,1))-min(fit_use{c,k}(:,1)))/99; %AZ to consistently have 100 points
    xi{k} = min(fit_use{c,k}(:,1)) : step : max(fit_use{c,k}(:,1));
    yi{k} = cum_gaussfit(psy_perf{k}, xi{k});
end
if length(unique_stim_type) == 3
    Thresh_pred = sqrt( Thresh_psy{1}^2*Thresh_psy{2}^2/(Thresh_psy{1}^2+Thresh_psy{2}^2) ); % now this is the prediction when there are three stimuli conditions
    % if length(unique_stim_type) == 3 && skip_combined == 0
    yi_pred = cum_gaussfit([Bias_psy{3},Thresh_pred], xi{3}); % smoothed line for prediction with PSE at actual combined condition
end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot psychometric, neurometric, CP over time
% % run the slide threshold over time, see whether performance fluctuate across time
% not work for coherence, temporarily 05/29/09
if overtimeplot == 1
    BegTrial_shift = BegTrial;
    EndTrial_shift = BegTrial_shift + span*one_repetition-1;
    n=1; % trial count during adaptation
    if EndTrial_shift > EndTrial %shift not calclated at all
        psy_thresh_shift=cell2mat(Thresh_psy');     psy_bias_shift=cell2mat(Bias_psy');
    else
        while EndTrial_shift <= EndTrial
            select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
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
                    % if sum(trials_shift)~=span, error('AZ wrong number of trials (compare=%u to span=%u)',sum(trials_shift),span); end %removed to get at least one value for psy_thresh_shift and psy_bias_shift
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
                    fit_data_psycho_cum_shift{k}(fit_data_psycho_cum_shift{k}(:,2)==0,2) = 1e-32; %AZ - updated a little
                    fit_data_psycho_cum_shift{k}(fit_data_psycho_cum_shift{k}(:,2)==1,2) = 1-(1e-32);
                end
                [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{k}); % to save time, use a different fit method
                psy_thresh_shift(k,n:(n+slide*length(unique_heading_shift)-1)) = tt;
                psy_bias_shift(k,n:(n+slide*length(unique_heading_shift)-1)) = bb;  % give the 'bias' value to "slide" trials used for the calculation. Later this is centered in the middle of "span" (all trials used for the calculation)
                %             wichman_psy = pfit(fit_data_psycho_cum_shift{k},'plot_opt','plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',1e-5,'sens',0,'compute_stats','false','verbose','false');
                %             psy_thresh_shift(k,n:(n+slide*one_repetition-1)) = wichman_psy.params.est(2);
                %             psy_bias_shift(k,n:(n+slide*one_repetition-1)) = wichman_psy.params.est(1);
            end
            BegTrial_shift = BegTrial_shift + slide*one_repetition;
            EndTrial_shift = BegTrial_shift + span*one_repetition-1;
            n = n + slide*length(unique_heading_shift); %AZ 1/11/12 changed here and above to be a function of trials and not reps (some experiments may have more headings per rep)
        end
        padding=nan(length(unique_condition_shift),length(unique_heading_shift)*(span-slide)/2);
        psy_thresh_shift=[padding psy_thresh_shift padding];%in order to center it correctly
        psy_bias_shift=[padding psy_bias_shift padding];%in order to center it correctly
    end
end

%AZ: I save a variable named with the file name (e.g. m26c0r185) since it is unique.
outfile = [RESULTS_DIR 'Psychometric_adapt_' FILE(strfind(FILE,'m'):strfind(FILE,'c')-1) '.mat']; %save each monkey to his own file
var=FILE(1:end-4);
stim={'Ves' 'Vis' 'Comb'};
for k = 1:length(unique_stim_type) %NB assuming one level of coherence
    eval(sprintf('%s.dir%s=fit_valid{1,k}(:,1);',var,stim{unique_stim_type(k)}))
    eval(sprintf('%s.right%s=fit_valid{1,k}(:,2);',var,stim{unique_stim_type(k)}))
    eval(sprintf('%s.count%s=fit_valid{1,k}(:,3);',var,stim{unique_stim_type(k)})) %added - # of trials for each heading
    eval(sprintf('%s.Xi%s=xi{k};',var,stim{unique_stim_type(k)}))
    eval(sprintf('%s.pfitcurve%s=cum_gaussfit(psy_perf{1,k}, xi{k});',var,stim{unique_stim_type(k)}))
    if exist('LFP_MVMT')==1
        for loc=1:size(LFP_MVMT,1)
            eval(sprintf('%s.LFP_%s{loc}=(LFP_MVMT{loc,k});',var,stim{unique_stim_type(k)}));
        end
    end %mean LFP matrix (headings,channel,LFPdata)
    for cell=1:cells
        eval(sprintf('%s.CP%s{cell}=CP{1,k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.CP_all_old%s{cell}=CP_all_old{1,k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.CP_all_new%s{cell}=CP_all_new{1,k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.FR%s{cell}=resp_mat{1,k,cell};',var,stim{unique_stim_type(k)}))%FR response of the neuron
        eval(sprintf('%s.FR_err%s{cell}=resp_mat_err{1,k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.line_coef%s{cell}=line_coef{k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.corr_coeff%s{cell}=corr_coeff{k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.p_value%s{cell}=p_value{k,cell};',var,stim{unique_stim_type(k)}))
        eval(sprintf('%s.sigm_center%s{cell}=sigm_center{k,cell};',var,stim{unique_stim_type(k)}))
        %         eval(sprintf('%s.raster_idx%s{cell}=uint32(find(temp_raster));',var,stim{unique_stim_type(k)})); %save raster as indices rather than 0s and 1s (save space). NB raster has half_guas + prePST*1000 before the actual response
        %         eval(sprintf('%s.raster_size%s{cell}=size(temp_raster);',var,stim{unique_stim_type(k)})); %save raster size so that able to decipher raster indices. NB raster has half_guas + prePST*1000 before the actual response
        if max(resp{1,k,cell}(:))<1 %the cell did not fire. 2011_12_07 changed from 'if max(resp_mat{1,k,cell})<1' because could have low FR here but may respond to other modality (not corrected here, just made less likely)
            eval(sprintf('%s.MVMTsdf%s{cell}(1:%u,1:(prePST*100 + unique_duration/10 + 100*postMVMT + 1))=NaN;',var,stim{unique_stim_type(k)},length(unique_heading))); %length = time of sdf
            eval(sprintf('%s.PEAKsdf%s{cell}(1:%u,1:size(resp{1,k,cell},2))=NaN;',var,stim{unique_stim_type(k)},length(unique_heading))); %length = reps
            eval(sprintf('%s.RMSsdf%s{cell}(1:%u,1:size(resp{1,k,cell},2))=NaN;',var,stim{unique_stim_type(k)},length(unique_heading))); %length = reps
            eval(sprintf('%s.SACCsdf%s{cell}(1:%u,:)=NaN;',var,stim{unique_stim_type(k)},length(unique_heading)));
            eval(sprintf('%s.FR_per_trial%s{cell}(1:%u,1:size(resp{1,k,cell},2))=NaN;',var,stim{unique_stim_type(k)},length(unique_heading)))
            eval(sprintf('%s.FR_per_choice%s{cell}(1:%u,1:size(resp_choice{1,k,cell},2))=NaN;',var,stim{unique_stim_type(k)},length(unique_heading)))
            for r=segval
                eval(sprintf('%s.FR_per_trial%u%s{cell}(1:%u,:)=NaN;',var,r,stim{unique_stim_type(k)},length(unique_heading)))
                eval(sprintf('%s.FR_per_choice%u%s{cell}(1:%u,:)=NaN;',var,stim{unique_stim_type(k)},length(unique_heading)))
            end
            eval(sprintf('%s.Baseline%s{cell}=[NaN NaN NaN NaN];',var,stim{unique_stim_type(k)})) %save different possible baselines
            eval(sprintf('%s.BaselineSD%s{cell}=[NaN NaN NaN NaN];',var,stim{unique_stim_type(k)}))
        else % calculated the FR per trial and spike density function SDF per heading (and peak SDF per trial)
            eval(sprintf('%s.FR_per_trial%s{cell}=resp{1,k,cell};',var,stim{unique_stim_type(k)}))%FR for each trial
            eval(sprintf('%s.FR_per_choice%s{cell}=resp_choice{1,k,cell};',var,stim{unique_stim_type(k)}))%FR for each choice (R/L)
            for r=segval
                eval(sprintf('%s.FR_per_trial%u%s{cell}=resp%u{1,k,cell};',var,r,stim{unique_stim_type(k)},r))%FR for each trial
                eval(sprintf('%s.FR_per_choice%u%s{cell}=resp_choice%u{1,k,cell};',var,r,stim{unique_stim_type(k)},r))%FR for each trial
            end
            for i = 1:length(unique_heading)
                temp_sdf{i}=downsample(convn(squeeze(raster{1,k,cell}(i,:,:))',gaus,'valid'),10); %'valid' - cuts off the extra 2X half_gaus that I added to the signal (to take care of edge effects). %downsample from 1000Hz to 100Hz
                pca_sdf(i,:,:)=temp_sdf{i}(prePST*100+(1:unique_duration/10),:); %100 and /10 because downsampled to 100Hz
                eval(sprintf('%s.MVMTsdf%s{cell}(i,:)=mean(temp_sdf{i}(1:(prePST*100 + unique_duration/10 + 100*postMVMT + 1),:),2);',var,stim{unique_stim_type(k)})); %movement SDF
                for tr=1:length(SACCdel{1,k,cell}(i,:)) %for this heading, find the sdf around the saccade (at different times for each trial)
                    times=prePST*100 + unique_duration/10 + round(100*SACCdel{1,k,cell}(i,tr)) + ((-100*preSACC):(100*postSACC)); %100 and /10 because downsampled to 100Hz
                    if max(times)<=size(temp_sdf{i},1)
                        temp_SACCsdf(:,tr)=temp_sdf{i}(times,tr);
                    else
                        temp_SACCsdf(:,tr)=NaN;
                        warning('Saccade time greater than sdf matrix limits')
                    end
                end
                eval(sprintf('%s.SACCsdf%s{cell}(i,:)=mean(temp_SACCsdf,2);',var,stim{unique_stim_type(k)})); %saccade SDF
            end
            %calculate baseline
            eval(sprintf('prePSTmedian=median(mean(%s.MVMTsdf%s{cell}(:,1:prePST*100)));',var,stim{unique_stim_type(k)})); %calculate the median of the prePST region of the SDF to be used as a baseline
            eval(sprintf('%s.Baseline%s{cell}=[mean(base250{1,k,cell}(:)) mean(base500{1,k,cell}(:)) mean(base1000{1,k,cell}(:)) prePSTmedian];',var,stim{unique_stim_type(k)})) %save many possible baselines
            eval(sprintf('%s.BaselineSD%s{cell}=[std(base250{1,k,cell}(:)) std(base500{1,k,cell}(:)) std(base1000{1,k,cell}(:)) NaN];',var,stim{unique_stim_type(k)}))
            eval(sprintf('%s.Base_per_trial%s{cell}=base1000{1,k,cell};',var,stim{unique_stim_type(k)}))%Baseline for each trial
            %ALT values
            for i = 1:length(unique_heading)
                eval(sprintf('%s.PEAKsdf%s{cell}(i,:)=max(temp_sdf{i}(prePST*100 + [(1+FR_del/10):(unique_duration/10)],:));',var,stim{unique_stim_type(k)}));
                eval(sprintf('%s.RMSsdf%s{cell}(i,:)=sqrt(mean((temp_sdf{i}(prePST*100 + [(1+FR_del/10):(unique_duration/10)],:)-%s.Baseline%s{cell}(base_use)).^2));',var,stim{unique_stim_type(k)},var,stim{unique_stim_type(k)}));
            end
            %PCA
            eval(sprintf('pca_sdf=pca_sdf-%s.Baseline%s{cell}(base_use);',var,stim{unique_stim_type(k)})) %subtract baseline for PCA (only relevant if UNcentered, o/w meaningless)
            collapse_sdf=reshape(shiftdim(pca_sdf,1),unique_duration/10,length(unique_heading)*repetition); %100 and /10 because downsampled to 100Hz
            PCAsdf_mean=mean(collapse_sdf');
            eval(sprintf('%s.PCAsdf_mean%s{cell}=PCAsdf_mean;',var,stim{unique_stim_type(k)}))
            sdf_heading_mean=mean(pca_sdf,3);
            %calculate the PCA across (averaged) headings
            [PCAcoeff,PCAscore,PCAlatent,PCAtsquared,PCAexplained] = pca(sdf_heading_mean,'NumComponents',3,'centered',PCAcentered);
            eval(sprintf('%s.PCAcoeff%s{cell}=PCAcoeff;',var,stim{unique_stim_type(k)}))
            eval(sprintf('%s.PCAscore%s{cell}=PCAscore;',var,stim{unique_stim_type(k)}))
            eval(sprintf('%s.PCAlatent%s{cell}=PCAlatent;',var,stim{unique_stim_type(k)}))
            eval(sprintf('%s.PCAtsquared%s{cell}=PCAtsquared;',var,stim{unique_stim_type(k)}))
            eval(sprintf('%s.PCAexplained%s{cell}=PCAexplained;',var,stim{unique_stim_type(k)}))
            eval(sprintf('[~,~,~,%s.PCAp_value%s{cell}]= xy_regress(unique_heading(abs(unique_heading)<head_cut),%s.PCAscore%s{cell}(abs(unique_heading)<head_cut,1));',var,stim{unique_stim_type(k)},var,stim{unique_stim_type(k)}))
            %using the PCA calculated (FROM BLOCK 1!) across averaged headings, calculate the projection for each trial (for the first PCA only)
            if ~isempty(precellID) & unique_stim_type(k)~=3 & any(precellID==cellID(cell)) %if this cell exists in the pre block. Do only for vis/vest
                eval(sprintf('prePCAcoeff=prePCAcoeff%s{precellID==cellID(cell)};',stim{unique_stim_type(k)}))
                eval(sprintf('prePCAsdf_mean=prePCAsdf_mean%s{precellID==cellID(cell)};',stim{unique_stim_type(k)}))
                if strcmp(PCAcentered,'on'),
                    X=collapse_sdf'-repmat(prePCAsdf_mean,size(collapse_sdf',1),1); %centered on the pre data
                    %TBD try this: X=collapse_sdf'-repmat(PCAsdf_mean,size(collapse_sdf',1),1); %centered on THIS (post) data
                else X=collapse_sdf'; end
                PCAscore_perTrial=reshape(X*prePCAcoeff(:,1),repetition,length(unique_heading))';
                eval(sprintf('%s.PCAscore_perTrial%s{cell}=PCAscore_perTrial;',var,stim{unique_stim_type(k)}))
                eval(sprintf('[~,~,~,%s.PCAp_value_perTrial%s{cell}]= xy_regress(repmat(unique_heading(abs(unique_heading)<head_cut),[repetition 1]),reshape(%s.PCAscore_perTrial%s{cell}(abs(unique_heading)<head_cut,:),[],1));',var,stim{unique_stim_type(k)},var,stim{unique_stim_type(k)}))
            else
                eval(sprintf('%s.PCAscore_perTrial%s{cell}=NaN;',var,stim{unique_stim_type(k)}))
                eval(sprintf('%s.PCAp_value_perTrial%s{cell}=NaN;',var,stim{unique_stim_type(k)}))
            end
            %caluculate PEAK and RMS p_values (for plotting purposes only)
            eval(sprintf('[~,~,~,%s.PEAKp_value_perTrial%s{cell}]= xy_regress(repmat(unique_heading(abs(unique_heading)<head_cut),[repetition 1]),reshape(%s.PEAKsdf%s{cell}(abs(unique_heading)<head_cut,:),[],1));',var,stim{unique_stim_type(k)},var,stim{unique_stim_type(k)}))
            eval(sprintf('[~,~,~,%s.RMSp_value_perTrial%s{cell}]= xy_regress(repmat(unique_heading(abs(unique_heading)<head_cut),[repetition 1]),reshape(%s.RMSsdf%s{cell}(abs(unique_heading)<head_cut,:),[],1));',var,stim{unique_stim_type(k)},var,stim{unique_stim_type(k)}))
        end
    end
end
eval(sprintf('%s.use_tempo_spikes=use_tempo_spikes;',var))
eval(sprintf('%s.cellID=cellID;',var))
eval(sprintf('%s.Bias_psy=Bias_psy;',var))
eval(sprintf('%s.Bias95CI_psy=Bias95CI_psy;',var))
eval(sprintf('%s.BiasSD_psy=BiasSD_psy;',var))
eval(sprintf('%s.Thresh_psy=Thresh_psy;',var))
eval(sprintf('%s.Thresh95CI_psy=Thresh95CI_psy;',var))
eval(sprintf('%s.gamma=gamma;',var))
eval(sprintf('%s.valid_thresh=valid_thresh;',var))
eval(sprintf('%s.Rsquare=Rsquare;',var))
eval(sprintf('%s.coherence=unique_motion_coherence;',var))
eval(sprintf('%s.unique_stim_type=unique_stim_type;',var))
eval(sprintf('%s.unique_heading=unique_heading;',var))
eval(sprintf('%s.psy_thresh_shift=psy_thresh_shift;',var))
eval(sprintf('%s.psy_bias_shift=psy_bias_shift;',var))
eval(sprintf('%s.delta=delta;',var))
eval(sprintf('%s.unique_delta=unique_delta;',var))
eval(sprintf('%s.num_trials=length(total_trials);',var))
eval(sprintf('%s.meanDelta_OnOff=meanDelta_OnOff;',var))
eval(sprintf('%s.VestVisual=VestVisual;',var))
eval(sprintf('%s.fixation=fixation;',var))
try eval(sprintf('%s.accelFB=accelFB;',var)); end %I only save this for neuronal data
try eval(sprintf('%s.vis_del=vis_del;',var)); end %I only save this for neuronal data
try eval(sprintf('%s.waveform=waveform;',var)); end %(a cellarray of length cells)I only save this for neuronal data
try eval(sprintf('%s.waveformFreq=waveformFreq;',var)); end %(a cellarray of length cells)I only save this for neuronal data

if (exist(outfile, 'file') == 0), save(outfile,var)    %file does not yet exist
else, save(outfile,var,'-append'); end

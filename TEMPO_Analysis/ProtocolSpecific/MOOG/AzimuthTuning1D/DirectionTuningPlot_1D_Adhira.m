% isolate more than 2 units from single unit recording data by offline
% spikesorting, analyze spike timing correlation among different
% units --SAM, 05/08
% %-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_1D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE)
%--------------------------------------------------------------------------
% process data into blocks
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_interocular_dist = data.moog_params(INTEROCULAR_DIST,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);

%get indices of any NULL conditions (for measuring spontaneous activity)
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);
interocular_dist = temp_interocular_dist(~null_trials & select_trials);
num_sigmas = temp_num_sigmas(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_motion_coherence = munique(motion_coherence');
unique_interocular_dist = munique(interocular_dist');
unique_num_sigmas = munique(num_sigmas');

% extract channel information
channelnum_temp = size(data.spike_rates);
channelnum = channelnum_temp(1,1); % how many channels
channelcount = 0;
for c = 1: channelnum
    spikesum = sum(data.spike_rates(c,:));
    if spikesum>0
        channelcount = channelcount+1;
    end
end

for c = 1 : 2 %only use last 2 channels
    temp_spike_rates = data.spike_rates(c+channelnum-2, :); 
    spike_rates = temp_spike_rates(~null_trials & select_trials);
    spike_rates_channel(c,:) = spike_rates;
    
    % check whether there is a slow drift of responsiveness over time
    [r,p] = corrcoef(1:length(temp_spike_rates), temp_spike_rates);
    R_slowdrift(c) = r(1,2);
    p_slowdrift(c) = p(1,2);     
    %now remove the slow drift responsiveness effet, comment out this part
    %if you don't want to, this is for z-scored only
    xx = 1:length(temp_spike_rates);
    sl = polyfit(xx, temp_spike_rates,1);
    temp_spike_rates_z = temp_spike_rates - (xx*sl(1)+sl(2)); % remove trend
    minf = min(temp_spike_rates_z);
    temp_spike_rates_z = temp_spike_rates_z + minf(1); % make firing rates >=0
    spike_rates_z = temp_spike_rates_z(~null_trials & select_trials);

    repetition = floor( length(spike_rates) / (1+length(unique_azimuth)*length(unique_stim_type)) ); % take minimum repetition
         
    % creat basic matrix represents each response vector
    resp = [];
    for k=1:length(unique_stim_type)        
        for i=1:length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i))  & (stim_type==unique_stim_type(k)) );
            for j = 1 : repetition; 
                spike_temp = spike_rates(select);   
                resp_trial{k}(j, i) = spike_temp( j );  
                resp_trial_plot{k}(j,i) = resp_trial{k}(j, i);
                resp_trial_channel{k,c}(j, i) = resp_trial{k}(j, i);
            end
            resp(i, k) = mean(spike_rates(select));        
            resp_std(i,k) = std(spike_rates(select));        
            resp_err(i,k) = std(spike_rates(select)) / sqrt(repetition); 
            
            % z-score data for spike count correlation analysis
            z_dist = spike_rates_z(select);
            if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
               z_dist = (z_dist - mean(z_dist))/std(z_dist);
            else
                z_dist = 0;
            end
            Z_Spikes(select) = z_dist;            
        end        
    end    
    Z_Spikes_channel(c,:) = Z_Spikes;

    % vectorsum and calculate preferred direction
    % vectors must be symetric, otherwise there will be a bias both on
    % preferred direction and the length of resultant vector
    % the following is to get rid off non-symetric data, hard code temporally
    if length(unique_azimuth) >8
        resp_s(1,:) = resp(1,:);
        resp_s(2,:) = resp(2,:);
        resp_s(3,:) = resp(4,:);
        resp_s(4,:) = resp(6,:);
        resp_s(5,:) = resp(7,:);
        resp_s(6,:) = resp(8,:);
        resp_s(7,:) = resp(9,:);
        resp_s(8,:) = resp(10,:);
    else
        resp_s(:,:) = resp(:,:);
    end
    unique_azimuth_s(1:8) = [0,45,90,135,180,225,270,315];
    unique_elevation_s(1:8) = 0;  
    resp_pair{c}(:,:) = resp(:,:);
    resp_err_pair{c}(:,:) = resp_err(:,:);
    
    for k = 1: length(unique_stim_type)
        [az(c,k), el(c,k), amp(c,k)] = vectorsumAngle(resp_s(:,k), unique_azimuth_s, unique_elevation_s);
        p_1D(c,k) = anova1(resp_trial{k},'','off');  
%         p_1D(c,k) = kruskalwallis(resp_trial{k});
    end  
end

for k=1:length(unique_stim_type) % ananlyze noise correlation in different conditions, if find no difference, combine later
    select_stim = logical( stim_type==unique_stim_type(k) );
    % noise correlation with stim type seperated
    [rr,pp] = corrcoef(Z_Spikes_channel(1,select_stim),Z_Spikes_channel(2,select_stim)); % temporarily hard coded
    noise_r_stim(k) = rr(1,2);
    noise_p_stim(k) = pp(1,2);     
    % this is only the regular correlation between two tuning curves
    [rr,pp] = corrcoef(resp_pair{1}(:,k),resp_pair{2}(:,k));
    corrcoef_r_unit(k) = rr(1,2);
    corrcoef_p_unit(k) = pp(1,2);
end

% Define figure
xoffset=0;
yoffset=0;
figure(2);
set(2,'Position', [5,15 980,650], 'Name', '1D Direction Tuning');
orient landscape;
set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'auto');

% temporarily hard coded, will be probematic if there are more than 12
% repetitions -GY
% Data only calculated for last 2 channels - AS edit
f{1}='bo-'; f{2}='ro-'; f{3}='go-'; f{4}='ko'; f{5}='co'; f{6}='mo';
for k=1: length(unique_stim_type)
    axes('position',[0.05+(k-1)*0.35 0.6 0.25 0.2]);
    for c = 1: 2 
%         if c~=2 % don't plot the second channel which is the synpulse signal
        errorbar(unique_azimuth, resp_pair{c}(:,k), resp_err_pair{c}(:,k), f{c} );
%         end
    hold on;
    end
    ylabel('spikes/s');
    xlabel('azimuth');
    xlim( [0, 315] );
    title(num2str(unique_stim_type(k)));
    set(gca, 'xtick',[unique_azimuth]);  
    
    % noise correlation
    axes('position',[0.1+(k-1)*0.35 0.38 0.15 0.18]);
    select_stim = logical( stim_type==unique_stim_type(k) );
    plot(Z_Spikes_channel(1,select_stim),Z_Spikes_channel(2,select_stim), 'o');
    xmin = min([Z_Spikes_channel(1,select_stim),Z_Spikes_channel(2,select_stim)]);
    xmax = max([Z_Spikes_channel(1,select_stim),Z_Spikes_channel(2,select_stim)]);
    xlim([xmin, xmax]);
    ylim([xmin, xmax]);
    hold on;
    plot([xmin,xmax], [xmin, xmax],'-');
end

%show file name and some values in text
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_stim_type)*(channelcount-1)] );
text(0, length(unique_stim_type)*(channelcount-1), FILE);
text(15,length(unique_stim_type)*(channelcount-1),'preferred           p                  noisecorrelation');
count = 0;
for k=1:length(unique_stim_type) 
    for c=1:2
        count = count+1;
%         if c>1 %skip second channel
%             cc = c+1;
%         else
%             cc = c;
%         end
        text(0,length(unique_stim_type)*(channelcount-1)-count, num2str(unique_stim_type(k)));
        text(15,length(unique_stim_type)*(channelcount-1)-count, num2str(az(c,k)) );
        text(25,length(unique_stim_type)*(channelcount-1)-count, num2str(p_1D(c,k)) );
        text(35,length(unique_stim_type)*(channelcount-1)-count, num2str(noise_r_stim(k)) );       
    end
end
axis off;
% ---------------------------------------------------------------------------------------
% % Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s\t'];
for i = 1 : 1000
      sprint_txt = [sprint_txt, ' %1.4f\t'];
end
% buff = sprintf(sprint_txt); %xcorr_stimtype(1,:), xcorr_stimtype(2,:) );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\PlexonPairwiseunitsSUSU_noisecorr.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t        neuro1 \t neuro2 \t stim \t pref1 \t p1 \t pref2 \t p2 \t DiffPref \t noisecorrelation');
    fprintf(fid, '\r\n');
end
% for pair=1:numpairs(1)
    for k=1:length(unique_stim_type)
        for c=1:2
%             if c==allpairs(pair,1)%finds pair data to output
                preferred(1) = az(c,k);
                pvalue(1)=p_1D(c,k);
%             elseif c==allpairs(pair,2)
                preferred(2) = az(c,k);
                pvalue(2)=p_1D(c,k);
%             end
        end
        prefdiff=abs(preferred(1)-preferred(2)); 
        buff = sprintf(sprint_txt, FILE, allpairs(pair,1), allpairs(pair,2), unique_stim_type(k), preferred(1), pvalue(1), preferred(2), pvalue(2), prefdiff , noise_r_stim(k), corrcoef_r_unit(k));
        fprintf(fid, buff);
        fprintf(fid, '\r\n');
    end
% end
fclose(fid);
return;
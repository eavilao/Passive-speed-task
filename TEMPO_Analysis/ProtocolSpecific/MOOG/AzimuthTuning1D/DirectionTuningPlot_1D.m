% % AZIMUTH_TUNING_1D.m -- Plots response as a function of azimuth for
% % heading task
% %--	YG, 07/12/04
% %-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_1D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);   

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
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

% descide whether loop is stim_type or disparity (when there is vary of
% disparity, it's the visual only condition

if length(unique_num_sigmas)>1
    condition_num = num_sigmas;
else
    condition_num = stim_type;
end
unique_condition_num = munique(condition_num');

% calculate spontaneous firing rate 
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found))
repetition = floor( length(spike_rates) / (length(unique_azimuth)*length(unique_stim_type)*length(unique_motion_coherence)) ); % take minimum repetition

% creat basic matrix represents each response vector
resp = [];
for c=1:length(unique_motion_coherence)  % motion coherence
    for k=1:length(unique_stim_type)
        for i=1:length(unique_azimuth)
            select = logical( azimuth==unique_azimuth(i) & stim_type==unique_stim_type(k) & motion_coherence==unique_motion_coherence(c) );
            if (sum(select) > 0)  % there are situations where visual/combined has >2 coherence and vestibular only has one coherence
                resp{c}(i, k) = mean(spike_rates(select));        
                resp_std{c}(i,k) = std(spike_rates(select));        
                resp_err{c}(i,k) = std(spike_rates(select)) / sqrt(repetition);
                
                spike_temp = spike_rates(select);
                resp_trial{c,k}(1:repetition,i) = spike_temp(1:repetition);
            else
                resp{c}(i, k) = resp{1}(i, k);    % actually duplicate vestibular condition    
                resp_std{c}(i,k) = resp_std{1}(i,k);        
                resp_err{c}(i,k) = resp_err{1}(i,k);
                
                resp_trial{c,k}(1:repetition,i) = resp_trial{1,k}(1:repetition,i);
            end
        end
    end
    % vectorsum and calculate preferred direction
    % vectors must be symetric, otherwise there will be a bias both on
    % preferred direction and the length of resultant vector
    % the following is to get rid off non-symetric data, hard code temporally
    if length(unique_azimuth) >8
        resp_s{c}(1,:) = resp{c}(1,:);
        resp_s{c}(2,:) = resp{c}(2,:);
        resp_s{c}(3,:) = resp{c}(4,:);
        resp_s{c}(4,:) = resp{c}(6,:);
        resp_s{c}(5,:) = resp{c}(7,:);
        resp_s{c}(6,:) = resp{c}(8,:);
        resp_s{c}(7,:) = resp{c}(9,:);
        resp_s{c}(8,:) = resp{c}(10,:);
    else
        resp_s{c}(:,:) = resp{c}(:,:);
    end
end

unique_azimuth_s(1:8) = [0,45,90,135,180,225,270,315];
unique_elevation_s(1:8) = 0;
unique_azimuth_plot = [unique_azimuth',360];

for c=1:length(unique_motion_coherence)
    for k = 1: length(unique_stim_type)
 %       if length(resp_s) >= length(unique_azimuth_s)
%             [az(c,k), el(c,k), amp(c,k)] = vectorsumAngle(resp_s{c}(:,k), unique_azimuth_s, unique_elevation_s);
        az(c,k) = 999; el(c,k) = 999; amp(c,k) = 999;
 %       else
 %           az(c,k) = NaN;   % this hard-coded method cannot handle < 8 azims  -CRF
 %       end
    %    Modulation Index
        DDI(c,k) = ( max(resp{c}(:,k))-min(resp{c}(:,k)) ) / ( max(resp{c}(:,k))-min(resp{c}(:,k))+mean(resp_std{c}(:,k)) );
        index_90 = find(unique_azimuth == 90);    
        Dprime(c,k) = (resp{c}(index_90+1,k)-resp{c}(index_90-1,k)) / sqrt( (resp_std{c}(index_90+1,k)^2+resp_std{c}(index_90-1,k)^2)/2 );
        p_1D(c,k) = anova1(resp_trial{c,k},'','off');
    end
end

% % ------------------------------------------------------------------
% Define figure
xoffset=0;
yoffset=0;
figure(3);
set(3,'Position', [5,15 980,650], 'Name', '1D Direction Tuning');
orient landscape;
set(0, 'DefaultAxesXTickMode', 'manual', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'manual');

spon_azimuth = min(unique_azimuth_plot) : 1 : max(unique_azimuth_plot);
% temporarily hard coded, will be probematic if there are more than 3*3 conditions
% repetitions -GY
f{1,1}='bo-'; f{1,2}='bo-'; f{1,3}='bo-'; 
f{2,1}='ro-'; f{2,2}='ro-'; f{2,3}='ro-'; 
f{3,1}='go-'; f{3,2}='go-'; f{3,3}='go-'; 

for k=1: length(unique_stim_type)     
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.45;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.64+yoffset 0.32 0.3]);
    for c=1:length(unique_motion_coherence)
        errorbar(unique_azimuth, resp{c}(:,k), resp_err{c}(:,k), f{c,k} );
        hold on;
    end
    plot(spon_azimuth, spon_resp, 'k--');

    ylabel('spikes/s');
    xlabel('azimuth');
    xlim( [min(unique_azimuth), max(unique_azimuth)] );
    title(num2str(unique_stim_type(k)));
    set(gca, 'xtick', unique_azimuth);

    xoffset=xoffset+0.48;    
end

%show file name and some values in text
axes('position',[0.5,0.1, 0.4,0.4] );
xlim( [0,100] );
ylim( [0,length(unique_stim_type)*length(unique_motion_coherence)+2] );
text(0, length(unique_stim_type)*length(unique_motion_coherence)+2, FILE);
text(40, length(unique_stim_type)*length(unique_motion_coherence)+2, 'SpikeChan=');
text(60, length(unique_stim_type)*length(unique_motion_coherence)+2, num2str(SpikeChan));
text(0,length(unique_stim_type)*length(unique_motion_coherence)+1,'stim coherence    prefer               DDI                  p');
count=0;
for k=1:length(unique_stim_type)  
    for c=1:length(unique_motion_coherence)     
        count=count+1;
        text(0,length(unique_stim_type)*length(unique_motion_coherence)-(count-1),num2str(unique_stim_type(k)));
        text(10,length(unique_stim_type)*length(unique_motion_coherence)-(count-1), num2str(unique_motion_coherence(c)) );              
        text(20,length(unique_stim_type)*length(unique_motion_coherence)-(count-1), num2str(az(c,k)) );
        text(50,length(unique_stim_type)*length(unique_motion_coherence)-(count-1), num2str(DDI(c,k)) );
        text(70,length(unique_stim_type)*length(unique_motion_coherence)-(count-1), num2str(p_1D(c,k)) );
    end
end
axis off;

%% ---------------------------------------------------------------------------------------
% Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s\t'];
for i = 1 : 200
     sprint_txt = [sprint_txt, ' %1.3f\t'];    
end
% if length(unique_stim_type)~=1
    buff = sprintf(sprint_txt, FILE, spon_resp, az, p_1D, Dprime, DDI );
%    buff = sprintf(sprint_txt, FILE, unique_stim_type, unique_motion_coherence, DDI );
    outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\DirectionTuning1D.dat'];
% else    
%     buff = sprintf(sprint_txt, FILE, spon_resp, az(:), amp(:), p_1D{:}, DI);
%     outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\DirectionTuning1D_Hui.dat'];
% end

% buff = sprintf(sprint_txt, FILE, az(:), p_1D{:},congruency, Z_Spikes );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\DirectionTuning1D_Zscore.dat'];
    
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
	if length(unique_stim_type)~=1
        %     fprintf(fid, 'FILE\t Spon\t HTIve\t HTIvi\t VesP\t VisP\t VesPre\t VisPre\t VesSlo\t Vis\Slo\t VesMin\t VisMin\t VesMax\t VisMax\t');
        fprintf(fid, 'FILE\t Spon\t VesPref\t VisPref\t VesP\t VisP\t VesDDI\t VisDDI\t Congruency');
	else
        fprintf(fid, 'FILE\t Spon\t Az_6\t Az_7_5\t Az_9\t Amp_6\t Amp_7_5\t Amp_9\t P_6\t P_7.5\t P_9\t DI_6\t DI_7.5\t DI_9\t');
	end
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
% SaveTrials(FILE,BegTrial,EndTrial,p_1D)
return;


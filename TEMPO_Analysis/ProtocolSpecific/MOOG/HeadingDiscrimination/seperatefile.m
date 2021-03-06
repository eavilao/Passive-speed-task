%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function seperatefile(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP


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

Nsep  =  floor(repetition/10);


for ii = 1 : Nsep
    clear eye_data  event_data moog_params   targ_params  misc_params  eye_positions
    if ii < Nsep
        sel = (ii-1)*10*one_repetition + 1 : ii*10*one_repetition; 
        eye_data = data.eye_data(:,:,sel);
        data1.eye_data = eye_data;
        event_data = data.event_data(:,:,sel);
        data1.event_data = event_data;
        targ_params = data.targ_params(:,sel,:);
        data1.targ_params = targ_params;
         misc_params = data.misc_params(:,sel);
        data1.misc_params = misc_params;
         moog_params = data.moog_params(:,sel,:);
        data1.moog_params = moog_params;
        eye_positions = data.eye_positions(:,sel);
        data1.eye_positions =eye_positions;
        Psychometric_FineRotDiscrm_2targ_deg1(data1, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, ii);
    else
        sel = (ii-1)*10*one_repetition + 1 : length(heading);
         eye_data = data.eye_data(:,:,sel);
        data1.eye_data = eye_data;
        event_data = data.event_data(:,:,sel);
        data1.event_data = event_data;
        targ_params = data.targ_params(:,sel,:);
        data1.targ_params = targ_params;
         misc_params = data.misc_params(:,sel);
        data1.misc_params = misc_params;
         moog_params = data.moog_params(:,sel,:);
        data1.moog_params = moog_params;
        eye_positions = data.eye_positions(:,sel);
        data1.eye_positions =eye_positions;
       Psychometric_FineRotDiscrm_2targ_deg1(data1, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, ii);
    end
end
    
 
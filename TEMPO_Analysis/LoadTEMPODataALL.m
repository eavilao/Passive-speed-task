%----------------------------------------------------------------------------------------------------------
%-- LoadTEMPOData.m: This function opens the TEMPO HTB file and reads in all of the data and header information
%--	contained therein.  All data is returned in the large structure arrays, good_data and bad_data.  By doing it this way,
%-- 	I am keeping every scrap of data about every trial all together in these big structures.  Note that various
% %--	defines in TEMPO_Defs.m are needed to access the data from good/bad_data.  GCD, starting 12/25/99
%----------------------------------------------------------------------------------------------------------
function	[return_value, good_data, bad_data, newChans,protocol] =LoadTEMPODataALL(PATH,FILE)

global binned_ced_spikes good_data listtext num_recorded_spike_channels
TEMPO_Defs;		%reads in some definitions that we need
% Path_Defs;		%some path definitions

good_data = [];
bad_data = [];
 
l = length(FILE);
if (FILE(l-3:l) == '.htb')	% .htb extension already there
    filename = [PATH FILE];   %the HTB data file
    logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
else	%no extension in FILE, add extensions
    filename = [PATH FILE '.htb'];   %the HTB data file
    logfile = [PATH FILE '.log'];   %the TEMPO log file
end
 
%disp('opening  HTB file...');
fid = htbOpen(filename);            % Open the HTB file
if (fid == -1)		%file could not be opened
    return_value = -1;
    return;
end

%disp('counting databases...');
ndbs = htbCount(fid);               % Get number of databases in it

%read in data and headers from the TEMPO .htb file
%identify type of each database and store it in appropriate array
temp1 = []; temp2 = []; temp3 = []; temp4 = [];
for database = 1:ndbs               % Loop though each one...
    %store the database headers in good_data and in bad_data
    temp_hd = htbGetHd(fid, database);

    if strcmp(temp_hd.title, 'Eye Traces')
        good_data.htb_header{EYE_DB} = temp_hd;
        bad_data.htb_header{EYE_DB} = temp_hd;
        %disp('reading eye data...');
        temp1 = htbGetDa(good_data.htb_header{EYE_DB});  %gets data for *all* epochs into a single matrix with dimentions: (sample#*trial#, channel#)
        temp1 = reshape(temp1', [good_data.htb_header{EYE_DB}.nchannels, good_data.htb_header{EYE_DB}.period, good_data.htb_header{EYE_DB}.sweep]);
        %*NOTE*:after this manipulation, eye_data has the following indices: (channel #, sample #, trial #).  I am doing this because it is 10x faster
        %than reading in the data for each trial using htbGetEp()
    end
    
    if strcmp(temp_hd.title, 'Spikes')
        good_data.htb_header{SPIKE_DB} = temp_hd;
        bad_data.htb_header{SPIKE_DB} = temp_hd;
        %disp('reading spike data...');
        temp2 = htbGetDa(good_data.htb_header{SPIKE_DB});
%                 temp2 = zeros(1000, 5);
        if size(temp2) == [good_data.htb_header{SPIKE_DB}.period * good_data.htb_header{SPIKE_DB}.sweep good_data.htb_header{SPIKE_DB}.nchannels]
            temp2 = reshape(temp2', [good_data.htb_header{SPIKE_DB}.nchannels, good_data.htb_header{SPIKE_DB}.period, good_data.htb_header{SPIKE_DB}.sweep]);
        else              % *********************************************************************************************************
            temp2 = [];   % CRF -- 3/8/06 -- TEMPORARY -- error in Tempo screwed up the spike DB in some of my data files (discrim_2I)
        end               % *********************************************************************************************************
        num_recorded_spike_channels = size(temp2, 1);
    end
    
    if strcmp(temp_hd.title, 'Events')
        good_data.htb_header{EVENT_DB} = temp_hd;
        bad_data.htb_header{EVENT_DB} = temp_hd;
        %disp('reading event data...');
        temp3 = htbGetDa(good_data.htb_header{EVENT_DB});
        temp3 = reshape(temp3', [good_data.htb_header{EVENT_DB}.nchannels, good_data.htb_header{EVENT_DB}.period, good_data.htb_header{EVENT_DB}.sweep]);
    end
    
    if ( strcmp(temp_hd.title, 'Local Field Potentials') | strcmp(temp_hd.title, 'LFP') )
        good_data.htb_header{LFP_DB} = temp_hd;
        bad_data.htb_header{LFP_DB} = temp_hd;
        disp('reading LFP data...');
        temp4 = htbGetDa(good_data.htb_header{LFP_DB});
        temp4 = reshape(temp4', [good_data.htb_header{LFP_DB}.nchannels, good_data.htb_header{LFP_DB}.period, good_data.htb_header{LFP_DB}.sweep]);
    end
    
    
end

err = htbClose(fid);                % Close HTB file

%determine number of trials so that we can declare arrays in ReadTEMPOLog
num_trials = good_data.htb_header{EVENT_DB}.sweep
%now, read in parameters from the TEMPO log file, which contains all experimental parameters for each trial
[revcorr_params, dots_params, moog_params, targ_params, cue_params, misc_params, one_time_params, eye_calib_params, obj_params, bar_params, bkgnd_params, neuron_params] = ReadTEMPOLog(logfile,num_trials);		
%dots_params[] has indices (param#, trial#, patch#)
%moog_params[] has indices (param#, trial#, item#)
%targ_params[] has indices (param#, trial#, targ#)
%cue_params[] has indices (param#, trial#, cue#)
%misc_params[] has indices (param#, trial#)
%one_time_params[] has index (param#)
%eye_calib_params[] has indices (param#, value#)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% added by BJP 1/3/01 to incorporate binding data
%obj_params[] has indices (param#, trial#, obj#)
%bar_params[] has indices (param#, trial#, obj#, bar#)
%bkgnd_params[] has indices (param#, trial#)
%neuron_params[] has indices (param#, #spike channels)

%now, remove any trials that were not completed successfully as we store the data in good_data
%and store data for trials not completed in bad_data
if isempty(temp2)    
    all_trials = 1:size(temp3,3);	%list of indices for all trials; CRF, 2/21/06 -- added IF statement because some of my behavior-only (discrim) datasets have no spike DB (and thus temp2 cannot be used to index all_trials)
else
    all_trials = 1:size(temp2,3);	%list of indices for all trials; GCD changed this 5/18/04 to deal with a strange data file that had extra trial in event channel; should be OK this way
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%changed by Xiongjie Yu on May 19, 2010
% %%%the following is th original code
% good_trials = find(temp3(:,:,:) == SUCCESS_CD);
%%the following is the changed code
[keys, data] = textread(logfile, '%s %[^\n]', 'bufsize', 300000);
num_lines_read = size(keys, 1);
kum=1;
protoc=sscanf(data{1},'%d');
protocol=protoc(1); 

%%%%%read sorted data and event data from 
spikefile = [PATH  'ALLspikes.txt']  
parafile = [PATH  'parameter.txt'] 


tf=exist(parafile);
if  exist(parafile)  &&  exist(spikefile) 
    allspikes=dlmread(spikefile);
    try
         parameter=dlmread(parafile,'\t',1,0);
    catch
        parameter=dlmread(parafile);
    end

    if exist([PATH  'add_spk.txt' ])
        addspk = dlmread([PATH  'add_spk.txt' ]);
%         for ii = 1 : length(addspk)
%             index = find( allspikes > addspk(ii));
%             if ~isempty(index)
%                 allspikes(index(1)+1: end+1) = allspikes(index(1): end);
%                 allspikes(index(1)) = addspk(ii);
%             end
%         end
        [temp11, temp22] = size(addspk);
        addspk = reshape(addspk, temp11*temp22,1);
        allspikes =  [allspikes;  addspk];
        allspikes = sort(allspikes);

    end
         trialstarttime=parameter(:,4);
         alldata=convert(allspikes,trialstarttime);
         spkdt = convert_1(allspikes,trialstarttime);
else
   
   importdata(PATH,FILE)
   allspikes=dlmread(spikefile);
    try
         parameter=dlmread(parafile,'\t',1,0);
    catch
        parameter=dlmread(parafile);
    end

    if exist([PATH  'add_spk.txt' ])
        addspk = dlmread([PATH  'add_spk.txt' ]);
        [temp11, temp22] = size(addspk);
        addspk = reshape(addspk, temp11*temp22,1);
        allspikes =  [allspikes;  addspk];
        allspikes = sort(allspikes);
    end
      trialstarttime=parameter(:,4);
      alldata=convert(allspikes,trialstarttime);
      spkdt =convert_1(allspikes,trialstarttime);
end


 

if(protocol==101)|(protocol==134)|(protocol==139)|(protocol==141)

good_trialsindex1=find(temp3(:,:,:) == IN_T1_WIN_CD);
good_trialsindex2=find(temp3(:,:,:) == IN_T2_WIN_CD);
good_trialsindex3=find(temp3(:,:,:) == SUCCESS_CD);
good_trialsnumber1 = ceil(good_trialsindex1/good_data.htb_header{EVENT_DB}.period);
good_trialsnumber2 = ceil(good_trialsindex2/good_data.htb_header{EVENT_DB}.period);
good_trialsnumber3 = ceil(good_trialsindex3/good_data.htb_header{EVENT_DB}.period);
good_trial1=false(num_trials,1);
good_trial2=false(num_trials,1); 
good_trial3=false(num_trials,1); 
good_trial1(good_trialsnumber1)=true;
good_trial2(good_trialsnumber2)=true;
good_trial3(good_trialsnumber3)=true;
good_trialsdex12 =xor(good_trial1,good_trial2);
good_trialsde=good_trialsdex12|good_trial3;

if(tf)
candidatetrial=parameter(:,6);
good_trialswithsorting=candidatetrial>0;
good_trialsdex=good_trialswithsorting&good_trialsde;
good_trials=find(good_trialsdex==true);
else
    good_trials=find(good_trialsde==true);
end



else
    
candidatetrial=parameter(:,6);
good_trials=find(candidatetrial>0);
end

% good_trialsfile=[PATH  'goodtrial_use.txt'] ;
% good_trials=dlmread(good_trialsfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% good_trials = ceil(good_trials/good_data.htb_header{EVENT_DB}.period);  %these are now trial indices
%NOTE: don't use htb_header.sweep hereafter for the number of trials, since this includes trials that were
%not completed successfully (e.g. breaks of fixation).  GCD, 12/28/99
all_trials(good_trials) = NaN;	%mark the good trials with NaNs
bad_trials = all_trials(~isnan(all_trials));  %and then the bad trials are the members of all_trials ~= NaN


%now, I am going to stash the data into the appropriate
%parts of my big data structures, called good_data and bad_data
if isempty(temp1)
    good_data.eye_data = [];
    bad_data.eye_data = [];
else
    %Convert the eye_data here from A/D units to degrees of visual angle
    %This is a change in the approach implemented to make software calibration easier.  GCD 12/28/00
    temp1(LEYE_H, :, :) = temp1(LEYE_H, :, :) .* (one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
    temp1(LEYE_V, :, :) = temp1(LEYE_V, :, :) .* (one_time_params(Y_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
    temp1(REYE_H, :, :) = temp1(REYE_H, :, :) .* (one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
    temp1(REYE_V, :, :) = temp1(REYE_V, :, :) .* (one_time_params(Y_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
    
    if (size(temp1,1) > 4)  % if there are D/A channels in addition to eye channels
   	 temp1(DA_H, :, :) = temp1(DA_H, :, :) .* (one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
   	 temp1(DA_V, :, :) = temp1(DA_V, :, :) .* (one_time_params(Y_DEG_FULL_SCALE)/one_time_params(AD_RANGE));       
    end
    
    %Do the software calibration of eye position signals here, if software calibration was used in collection
    %Here, I operate directly on the temp1 array, replacing the previous values
    %Added by GCD, 12/22/00
    if (one_time_params(SOFTWARE_CALIB_STATUS) == 1)
        junk = temp1;
        if isnan(eye_calib_params(SOFT_CAL_LEYE_HORIZ,4)) %if only 3 calibration parameters
            temp1(LEYE_H, :, :) = eye_calib_params(SOFT_CAL_LEYE_HORIZ,1) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,2).*junk(LEYE_H,:,:) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,3).*junk(LEYE_V,:,:);
            temp1(LEYE_V, :, :) = eye_calib_params(SOFT_CAL_LEYE_VERT,1) + eye_calib_params(SOFT_CAL_LEYE_VERT,2).*junk(LEYE_V,:,:) + eye_calib_params(SOFT_CAL_LEYE_VERT,3).*junk(LEYE_H,:,:);
            temp1(REYE_H, :, :) = eye_calib_params(SOFT_CAL_REYE_HORIZ,1) + eye_calib_params(SOFT_CAL_REYE_HORIZ,2).*junk(REYE_H,:,:) + eye_calib_params(SOFT_CAL_REYE_HORIZ,3).*junk(REYE_V,:,:);
            temp1(REYE_V, :, :) = eye_calib_params(SOFT_CAL_REYE_VERT,1) + eye_calib_params(SOFT_CAL_REYE_VERT,2).*junk(REYE_V,:,:) + eye_calib_params(SOFT_CAL_REYE_VERT,3).*junk(REYE_H,:,:);
        else %if all 4 parameters are used
            temp1(LEYE_H, :, :) = eye_calib_params(SOFT_CAL_LEYE_HORIZ,1) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,2).*junk(LEYE_H,:,:) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,3).*junk(LEYE_V,:,:) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,4).*junk(LEYE_V,:,:).*junk(LEYE_H,:,:);
            temp1(LEYE_V, :, :) = eye_calib_params(SOFT_CAL_LEYE_VERT,1) + eye_calib_params(SOFT_CAL_LEYE_VERT,2).*junk(LEYE_V,:,:) + eye_calib_params(SOFT_CAL_LEYE_VERT,3).*junk(LEYE_H,:,:) + eye_calib_params(SOFT_CAL_LEYE_VERT,4).*junk(LEYE_H,:,:).*junk(LEYE_V,:,:);
            temp1(REYE_H, :, :) = eye_calib_params(SOFT_CAL_REYE_HORIZ,1) + eye_calib_params(SOFT_CAL_REYE_HORIZ,2).*junk(REYE_H,:,:) + eye_calib_params(SOFT_CAL_REYE_HORIZ,3).*junk(REYE_V,:,:) + eye_calib_params(SOFT_CAL_REYE_HORIZ,4).*junk(REYE_V,:,:).*junk(REYE_H,:,:);
            temp1(REYE_V, :, :) = eye_calib_params(SOFT_CAL_REYE_VERT,1) + eye_calib_params(SOFT_CAL_REYE_VERT,2).*junk(REYE_V,:,:) + eye_calib_params(SOFT_CAL_REYE_VERT,3).*junk(REYE_H,:,:) + eye_calib_params(SOFT_CAL_REYE_VERT,4).*junk(REYE_H,:,:).*junk(REYE_V,:,:);
        end
        ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
        listtext{length(listtext)+1} = 'Using software-calibrated eye data.';
        set(ListHandle, 'String', listtext);
    end
    good_data.eye_data = temp1(:,:,good_trials);
    bad_data.eye_data = temp1(:,:,bad_trials);
end
clear temp1;

if isempty(temp2)
    good_data.spike_data = [];
    bad_data.spike_data = [];
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%      his is comented by xjy and changed to the below code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%T
%     good_data.spike_data = temp2(:,:,good_trials);
%     good_data.spike_data = cat(1,temp2(:,:,good_trials),alldata(:,:,good_trials) );
if(tf)
    good_data.spike_data = cat(1,temp2(:,:,good_trials),alldata(:,:,good_trials) );
else
 good_data.spike_data = temp2(:,:,temp2(:,:,good_trials));
end

for ii = 1 : length(good_trials)
    spkdts{ii} =  spkdt{good_trials(ii)};
end

good_data.spkdt = spkdts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
bad_data.spike_data = temp2(:,:,bad_trials);    
end
clear temp2;

if isempty(temp3)
    good_data.event_data = [];
    bad_data.event_data = [];
else
    good_data.event_data = temp3(:,:,good_trials);
    bad_data.event_data = temp3(:,:,bad_trials);
end
clear temp3;

if isempty(temp4)
    good_data.lfp_data = [];
    bad_data.lfp_data = [];
else
    for i = 1:size(temp4,1)
        temp4(i, :, :) = temp4(i, :, :)/one_time_params(AD_RANGE);
    end
    good_data.lfp_data = temp4(:,:,good_trials);
    bad_data.lfp_data = temp4(:,:,bad_trials);
end
clear temp4;

good_data.targ_params = targ_params(:,good_trials,:);
good_data.misc_params = misc_params(:,good_trials,:);
good_data.one_time_params = one_time_params;
good_data.eye_calib_params = eye_calib_params;
good_data.neuron_params = neuron_params;		%neuron params array is used for both binding and dots protocols

if (~isempty(revcorr_params))
    good_data.revcorr_params = revcorr_params(:, good_trials, :);
end
if (~isempty(moog_params))
    good_data.moog_params = moog_params(:, good_trials, :);
end
if ~isempty(dots_params)
    good_data.dots_params = dots_params(:,good_trials,:);
end
if ~isempty(cue_params)
    good_data.cue_params = cue_params(:,good_trials,:);
end
if (~isempty(obj_params))
    good_data.obj_params = obj_params(:,good_trials,:);
end
if (~isempty(bar_params))
   good_data.bar_params = bar_params(:,good_trials,:,:);
end
if (~isempty(bkgnd_params))
    good_data.bkgnd_params = bkgnd_params(:,good_trials,:);
end
%else%if ~isempty(moog_params)
    %good_data.moog_params = moog_params(:,good_trials,:);
    %else
    %binding protocol values
%    good_data.obj_params = obj_params(:,good_trials,:);
%    good_data.bar_params = bar_params(:,good_trials,:,:);
%   good_data.bkgnd_params = bkgnd_params(:,good_trials,:);
%end
 
bad_data.targ_params = targ_params(:,bad_trials,:);
bad_data.misc_params = misc_params(:,bad_trials,:);
bad_data.one_time_params = one_time_params;
bad_data.eye_calib_params = eye_calib_params;
bad_data.neuron_params = neuron_params;

if (~isempty(revcorr_params))
	bad_data.revcorr_params = revcorr_params(:, bad_trials, :);
end
if ~isempty(cue_params)
    bad_data.cue_params = cue_params(:,bad_trials,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% added by BJP 1/3/01 for binding stuff
if ~isempty(dots_params)
    bad_data.dots_params = dots_params(:,bad_trials,:);
elseif ~isempty(moog_params)
    bad_data.moog_params = moog_params(:,bad_trials,:);
else
    bad_data.obj_params = obj_params(:,bad_trials,:);
    bad_data.bar_params = bar_params(:,bad_trials,:,:);
    bad_data.bkgnd_params = bkgnd_params(:,bad_trials,:);
end

%%%%the follwing is commented by xjy
% % Check to see if the htb file was in a standard location.
% PATH2 = upper(PATH);
% wasinstdloc = findstr(PATH2, 'Z:\DATA\TEMPO');
% 
% % This adds the good data from spikesort2.
% fa = find(PATH == '\');
% 
% % Prevent and indexing error if the htb file was in a nonstandard location.
% if (~isempty(wasinstdloc))    
%         fName = ['Z:\Data\CED', PATH(fa(3):fa(4)), 'SortedSpikes\', FILE(1:length(FILE) - 4)];
% else
%     if findstr(PATH, 'Z:\Data\MOOG')%Add by AHC
%         fName = ['Z:\Data\MOOG', PATH(fa(3):fa(4)), 'CED\SortedSpikes\', FILE(1:length(FILE) - 4)]; %Add by AHC
%         good_data.fName=fName;
%     else      
%         fName = [];
%     end
% end
%  fName = ['C:\data\sorted\', FILE(1:length(FILE) - 4)]
% 
newChans=0;
  newChans = 0;
% check_spont = 0;
% global LOAD_MU_DATA;
% if (LOAD_MU_DATA & exist([fName, '.mat'], 'file'))
%      disp('Running Packaroni, Booyah!!!');
%      [good_data winDiscrim] = Packaroni(good_data, fName, 1, 0);  %get window discriminator data
%      [good_data shiftValue]= SpikeCorr(good_data, 1, 1, winDiscrim);
%      if  ( (shiftValue == -50) & (sum( sum( good_data.spike_data(winDiscrim,:,:) ) == 0) ) )
%          % remove window discrimination data if not being used
%          good_data.spike_data = good_data.spike_data(1:size(good_data.spike_data,1) - 1,:,:);
%          shiftValue = 0;
%      end
%      [good_data newChans] = Packaroni(good_data, fName, 0, shiftValue);  %get rest of data
%      for i = 1:length(newChans)
%         [good_data shiftValue]= SpikeCorr(good_data, 0, 1, newChans(i));
%         if check_spont == 1
%             num_temp = size(good_data.event_data);
%             num_t = num_temp(3);
%             spont = 0;
%             for j=1:num_t
%                 start_code_bin = find(good_data.event_data(1, :, j) == 3);
%                 stop_code_bin = find(good_data.event_data(1, :, j) == 4);
%                 spikes = sum(good_data.spike_data(newChans(i), start_code_bin:stop_code_bin, j));
%                 time = stop_code_bin - start_code_bin;
%                 spont = spikes/time + spont;
%             end
%             spont_avg = spont/j * 1000
%         end
%     end    
% end

%%%%the upper is commented by xjy

%add by AHC 02-22-06
% sName = ['C:\data\sorted\', FILE(1:length(FILE) - 4),'.mat']
% if (exist(sName,'file'))
%     s=1
%    disp('Loading Spike2 sorted data!')
%    [good_data SortedChannel]=PackData(good_data,sName);%get the Spike2 sorted data 
% end

%%%%%%%%%%%%%%%%%%%xyj commented%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added2loadtempodata
%%%%%%%%%%%%%%%%%%%xyj%%%%%%%%%%%%%%%%%%%commented%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return_value = 1;		%indicates completed OK
 
return;   
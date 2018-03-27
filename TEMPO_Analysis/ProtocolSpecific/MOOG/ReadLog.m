%----------------------------------------------------------------------------------------------------------
%-- ReadTEMPOLog.m: This function reads in all the data from the TEMPO *.log file, which contains a list
%--	of experiment parameters for each trial.	xjy 022810
%----------------------------------------------------------------------------------------------------------
function trialpara = ReadLog(logfile);

TEMPO_Defs;	%some defines that we'll need

%% main bottleneck is reading and parsing data out into these two cell arnrays - BJP 5/2/01
logfile
[keys, data] = textread(logfile, '%s %[^\n]', 'bufsize', 300000);
num_lines_read = size(keys, 1);
kum=1;
protoc=sscanf(data{1},'%d');
protocol=protoc(1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(protocol==112)
for i=1:num_lines_read
    
    
    switch keys{i}
    case 'TRIAL#'
        curr_trial( kum) = sscanf(data{i},'%d');  
          
    case 'OUTCOME'
        curr_out( kum) = sscanf(data{i},'%d');  
 
    case 'ROT_ELEVATION'
        t=sscanf(data{i},'%d');
         curr_elev(kum)=t(1);
         
             case 'ROT_AZIMUTH'
        t=sscanf(data{i},'%d');
        curr_azim(kum)=t(1);
        kum=kum+1;
 
    end
   
end
 
trialpara=[curr_trial',curr_elev',curr_azim',curr_out'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(protocol==100)
for i=1:num_lines_read
    
    
    switch keys{i}
    case 'TRIAL#'
        curr_trial( kum) = sscanf(data{i},'%d');  
          

    case 'OUTCOME'
        curr_out( kum) = sscanf(data{i},'%d');  
 
 
    case 'ELEVATION'
        t=sscanf(data{i},'%d');
         curr_elev(kum)=t(1);
         kum=kum+1;
         
             case 'AZIMUTH'
        t=sscanf(data{i},'%d');
        curr_azim(kum)=t(1);
        
 
    end
   
end

trialpara=[curr_trial',curr_elev',curr_azim',curr_out'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(protocol==122|protocol==123|protocol==124|protocol==125|protocol==126|protocol==127|protocol==128|protocol==129|protocol==130)  %SINUSOID_YAW_P5HZ
    
for i=1:num_lines_read
    
    
    switch keys{i}
    case 'TRIAL#'
        curr_trial( kum) = sscanf(data{i},'%d');  
          

    case 'OUTCOME'
        curr_out( kum) = sscanf(data{i},'%d');  
 curr_out( kum)=protocol;
 
    case 'SIN_ROT_AMPLITUDE'
        t=sscanf(data{i},'%f');
         curr_elev(kum)=t(1);
         kum=kum+1;
         
             case 'SIN_FREQUENCY'
        t=sscanf(data{i},'%f');
        curr_azim(kum)=t(1);
        
 
    end
   
end

trialpara=[curr_trial',curr_elev',curr_azim',curr_out'];
end
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(protocol==133|protocol==132|protocol==131)  %Tanslatoin threshold detection
    
for i=1:num_lines_read
    
    
    switch keys{i}
    case 'TRIAL#'
        curr_trial( kum) = sscanf(data{i},'%d');  
          

    case 'OUTCOME'
        curr_out( kum) = sscanf(data{i},'%d');  
 curr_out( kum)=protocol;
 
    case 'SIN_TRANS_AMPLITUDE'
        t=sscanf(data{i},'%f');
         curr_elev(kum)=t(1);
         kum=kum+1;
         
             case 'SIN_FREQUENCY'
        t=sscanf(data{i},'%f');
        curr_azim(kum)=t(1);
        
 
    end
   
end

trialpara=[curr_trial',curr_elev',curr_azim',curr_out'];
end
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(protocol==101)
     for i=1:num_lines_read
  
    
       switch keys{i}
      case 'TRIAL#'
        curr_trial( kum) = sscanf(data{i},'%d');  

    case 'OUTCOME'
        curr_out( kum) = sscanf(data{i},'%d');  
 

        case 'HEADING'
        t=sscanf(data{i},'%f') ;
        curr_headintangle(kum)=t(1) ;
         kum=kum+1;
 
        end
   
     end
     
curr_headintangle
 l = length(logfile);
htbfile = [logfile(1:l-4) '.htb']; 
%disp('opening  HTB file...');
fid = htbOpen(htbfile);            % Open the HTB file
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

  TEMPO_Defs;
    if strcmp(temp_hd.title, 'Events')
        good_data.htb_header{EVENT_DB} = temp_hd;
        bad_data.htb_header{EVENT_DB} = temp_hd;
        %disp('reading event data...');
        temp3 = htbGetDa(good_data.htb_header{EVENT_DB});
        temp3 = reshape(temp3', [good_data.htb_header{EVENT_DB}.nchannels, good_data.htb_header{EVENT_DB}.period, good_data.htb_header{EVENT_DB}.sweep]);
    end
    
    
    
end

err = htbClose(fid);                % Close HTB file

%determine number of trials so that we can declare arrays in ReadTEMPOLog
num_trials = good_data.htb_header{EVENT_DB}.sweep

%%%%%%changed by Xiongjie Yu on May 19, 2010
good_trialsindex1=find(temp3(:,:,:) == IN_T1_WIN_CD);
good_trialsindex2=find(temp3(:,:,:) == IN_T2_WIN_CD);
good_trialsindex3=find(temp3(:,:,:) == SUCCESS_CD);
good_trialsnumber1 = ceil(good_trialsindex1/good_data.htb_header{EVENT_DB}.period);
good_trialsnumber2 = ceil(good_trialsindex2/good_data.htb_header{EVENT_DB}.period);
good_trial1=false(num_trials,1);
good_trial2=false(num_trials,1); 
good_trial1(good_trialsnumber1)=true;
good_trial2(good_trialsnumber2)=true;
good_trialsdex12 =xor(good_trial1,good_trial2);
good_trialsdex=good_trialsdex12 ;
good_trials=find(good_trialsdex==true);

%%%%%%%%%%add the choice that monkey made
choice=nan(num_trials,1);
choice(good_trial1)=IN_T1_WIN_CD;
choice(good_trial2)=IN_T2_WIN_CD;
choice(good_trial1&good_trial2)=nan; %%%nan is the indicator that show the trial is not correct.

  trialpara=[curr_trial',choice,curr_headintangle',curr_out'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(protocol==134)|(protocol==141)
     for i=1:num_lines_read
  
    
       switch keys{i}
      case 'TRIAL#'
        curr_trial( kum) = sscanf(data{i},'%d');  

    case 'OUTCOME'
        curr_out( kum) = sscanf(data{i},'%d');  
 

        case 'HEADING'
        t=sscanf(data{i},'%f') ;
        curr_headintangle(kum)=t(1) ;
         kum=kum+1;
 
        end
   
     end
     
curr_headintangle
 l = length(logfile);
htbfile = [logfile(1:l-4) '.htb']; 
%disp('opening  HTB file...');
fid = htbOpen(htbfile);            % Open the HTB file
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

  TEMPO_Defs;
    if strcmp(temp_hd.title, 'Events')
        good_data.htb_header{EVENT_DB} = temp_hd;
        bad_data.htb_header{EVENT_DB} = temp_hd;
        %disp('reading event data...');
        temp3 = htbGetDa(good_data.htb_header{EVENT_DB});
        temp3 = reshape(temp3', [good_data.htb_header{EVENT_DB}.nchannels, good_data.htb_header{EVENT_DB}.period, good_data.htb_header{EVENT_DB}.sweep]);
    end
    
    
    
end

err = htbClose(fid);                % Close HTB file

%determine number of trials so that we can declare arrays in ReadTEMPOLog
num_trials = good_data.htb_header{EVENT_DB}.sweep

%%%%%%changed by Xiongjie Yu on May 19, 2010
good_trialsindex1=find(temp3(:,:,:) == IN_T1_WIN_CD);
good_trialsindex2=find(temp3(:,:,:) == IN_T2_WIN_CD);
good_trialsindex3=find(temp3(:,:,:) == SUCCESS_CD);
good_trialsnumber1 = ceil(good_trialsindex1/good_data.htb_header{EVENT_DB}.period);
good_trialsnumber2 = ceil(good_trialsindex2/good_data.htb_header{EVENT_DB}.period);
good_trial1=false(num_trials,1);
good_trial2=false(num_trials,1); 
good_trial1(good_trialsnumber1)=true;
good_trial2(good_trialsnumber2)=true;
good_trialsdex12 =xor(good_trial1,good_trial2);
good_trialsdex=good_trialsdex12 ;
good_trials=find(good_trialsdex==true);

%%%%%%%%%%add the choice that monkey made
choice=nan(num_trials,1);
choice(good_trial1)=IN_T1_WIN_CD;
choice(good_trial2)=IN_T2_WIN_CD;
choice(good_trial1&good_trial2)=nan; %%%nan is the indicator that show the trial is not correct.

  trialpara=[curr_trial',choice,curr_headintangle',curr_out'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
return;
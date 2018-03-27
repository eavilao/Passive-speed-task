function  correct_time(PATH,FILE)
  
 FNAME =  [FILE(1:end-4)  '.smr'];
 PNAME  = PATH;   
 FILENAME=[ PNAME ,FNAME];

addpath('Z:\Users\Xiongjie\Yaw discrimilation\afferent code\son');
addpath('Z:\Users\Xiongjie\Yaw discrimilation\afferent code\son\SON32'); 
addpath('Z:\Users\Xiongjie\Yaw discrimilation\afferent code\son\SON32\C_Sources');

channel_number=4;
fid=fopen( FILENAME,'r');
[ts,h]=SONGetChannel(fid,channel_number);

if(isstruct(ts))    
sortcode=ts.markers;
spiketim=ts.timings;

    
clear ts
clear h
  
sortex=unique(sortcode) ;
sortexhadl=channel_number;

SORTSELCT=1 ;
 
sortcodeindex=false(1, length(sortcode));
 
for i=1:length(SORTSELCT)
sortcodeindexcrrent=(sortcode==SORTSELCT(i));
 
sortcodeindex(sortcodeindexcrrent)=true;
 
clear sortcodeindexcrrent
end
SPIKETIME=spiketim(sortcodeindex);
clear sortcodeindex

else
   SPIKETIME= ts;
    
end



    SPIKETIME=sort(SPIKETIME) ; 
 %%%%%%###$$#$#$#$#$&***&&*&*&*&*&*&*&*&*&*&**&&*#$#$#$#$#$#$###$##$##$#$##$$#$##$#$    
 %%%%%%stor all the informaiton
      k = findstr(PNAME , '\');
      ksize=length(k)';
      parerentfold=PNAME (1:k(ksize-1)) ;

     addfolder=[ parerentfold,'batch TXT data','\',FNAME(1:length(FNAME)-4),'\'] 
    if(~isdir(addfolder))
    mkdir(addfolder)
    end
Allspikefile=[ addfolder,'ALLspikes' '.txt'];
ALLspikecrtfold=[PNAME ,'ALLspikes' '.txt'];
% dlmwrite(Allspikefile, SPIKETIME,'precision', 8)  
     dlmwrite(ALLspikecrtfold, SPIKETIME,'precision', 8)  
      
 %###############################     
      
%    onsetsynpulse = findobj(gcbf, 'Tag', 'Onset Option') ;
% onsetoptionhandle=get(onsetsynpulse, 'SelectedObject') ;
% =get(onsetoptionhandle, 'String');
onsetoptionstring = 'Event Code';
if(strcmp(onsetoptionstring, 'Synpulse'))
    [ts,T,t0]=load_named_event_channel_from_smr(FILENAME,'SynPulse');
sizeofsyn=length(ts);
synpusbetw=ts(1:sizeofsyn)-[0;ts(1:sizeofsyn-1)];
trialstartime=ts(synpusbetw>2);
end
 
 if(strcmp(onsetoptionstring, 'Event Code'))
try
     [ts,T,t0]=load_event11(FILENAME,'DigEvents'); 


catch  
   [ts,T,t0]=load_event11(FILENAME,'DigEvent');
end
eventcodes=ts.markers(:,1);
startcodes=find(eventcodes==4);
trialstartime=ts.timings(startcodes); 

 end
  
 
 %###############################   
 
 
  l = length(FNAME);
 if(iscell(FNAME))
        fnamenow=FNAME{1};
        forsortfile=findstr(fnamenow,'_');
       if(isempty(forsortfile))
        logfile = [PNAME  fnamenow(1:l-4) '.log'];   %the TEMPO log file
       else
      logfile = [PNAME  fnamenow(1:forsortfile-1) '.log'] ;
       end 
      
 else
        forsortfile=findstr(FNAME,'_');
        if(isempty(forsortfile))
       logfile = [PNAME  FNAME(1:l-4) '.log'];   %the TEMPO log file
        else
        logfile = [PNAME  FNAME(1:forsortfile-1) '.log'] ;
         end 
 end
 
 
trialpar = ReadLog(logfile);
[keys, data] = textread(logfile, '%s %[^\n]', 'bufsize', 300000);
num_lines_read = size(keys, 1);
kum=1;
protoc=sscanf(data{1},'%d');
PROTOCOL=protoc(1);
 clear keys data 
sweeplimit=size(trialpar, 1) ;
num_trials=length(trialstartime);
trial_selected=false(num_trials,1);
trial_selected(1:sweeplimit)=true;



forbatchfile=[addfolder,'parameter' '.txt'];  
forcurentfoldfile=[PNAME ,'parameter' '.txt']; 
 outcomeoftrial=trialpar(:,4);
 outcomecomparefile=[PNAME ,'out_comp' '.xls'];  


trialparsize=size(trialpar);
trialstartimesize=size(trialstartime); 

if(trialparsize(1)<trialstartimesize(1))
    diffsize=trialstartimesize-trialparsize;
    trialpar=[trialpar;nan(diffsize,4)];
end 
%  trialstartime=[0;trialstartime]
  try,
          ALLINFORMATION=[ trialpar(:,1:3),trialstartime,trialpar(:,4),trial_selected];
  catch,
      NN = min([size(trialpar, 1)  size(trialstartime , 1)    size(trial_selected , 1) ] );
       ALLINFORMATION=[ trialpar(1 : NN,1:3), trialstartime(1 : NN),trialpar(1 : NN,4),  trial_selected(1 : NN)  ];   
  end
%           ALLINFORMATION=[trialpar(1:117,1:3),trialstartime(1:117),trialpar(1:117,4),trial_selected(1:117)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
logpath=PNAME ;
file_name='parameter.txt';                      
formatstr='%12s %12s%12s%12s%12s%12s';
info_data=ALLINFORMATION;
formatdata='%12.8f\t %12.8f\t %12.8f\t  %12.8f\t  %12.8f \t %12.8f ';
if(PROTOCOL==100|PROTOCOL==112) %%%%3-D TUNING
info_str=sprintf('%12s\t%12s\t%12s\t%12s\t%12s\t%12s','Trial number','Elevation','Azimuth','Trial onset','Trial outcome','Trial selected index');
else
    if(PROTOCOL==101)|(PROTOCOL==134)%%%%Heading
    info_str=sprintf('%12s\t%12s\t%12s\t%12s\t%12s\t%12s','Trial number','Choice','Heading angle','Trial onset','Trial outcome','Trial selected index');
    else
     info_str=sprintf('%12s\t%12s\t%12s\t%12s\t%12s\t%12s','Trial number','SIN_ R/T_AMPLITUDE','SIN_FREQUENCY','Trial onset','Protocol','Trial selected index');   
    end
end

%%%%this m code convert the spike train to the format of tempo, the nth
%%%%trial 's start time is the starttime(n);
trialnum=length(starttime);
data=nan(1,5000,trialnum);
edges=((0:1:5000)-995)/1000;
minbin=min(edges);
maxbin=max(edges);
for i=1:trialnum
    
    spikenow=spike-starttime(i);
    spikenowtouse=spikenow((spikenow>=minbin)&(spikenow<=maxbin));
    clear spikenow
    if(~isempty(spikenowtouse))
 currentconvert=histc(spikenowtouse,edges);
 data(1,:,i)=currentconvert(1:5000);
    else
   data(1,:,i)=zeros(5000,1);
    end
clear spikenowtouse
clear currentconvert
end
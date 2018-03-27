function PeakTimes = f3dsimple_spikes_aff(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath, WindowInterval,  m_pval);
% get each trials of spikes


% close(2);
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust the system time by 115
StartEventBin = StartEventBin +115;
StopEventBin  = StopEventBin + 115;
Step=25;

PeakTimes = [0 0];


%get the column of values for azimuth and elevation and stim_type
switch Protocol
    case 100 % DIRECTION_TUNING_3D 
        temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ELEVATION,:,MOOG);
        CurrentProtocol=['Translation'];
    case 112 %ROTATION_TUNING_3D
        temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
        CurrentProtocol=['Rotation'];
    case 104 %DIR3D_VARY_FIXATION 
        temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ELEVATION,:,MOOG);
        CurrentProtocol=['DIR3D VARY FIXATION '];
    case 107
        temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ELEVATION,:,MOOG);
        CurrentProtocol=['Translation 1D'];
end

temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);
% temp_spike_rates = data.spike_rates(SpikeChan, :);    

% add
% 4/17/09~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m_path = [PATH]
spikefile = [PATH  'ALLspikes.txt'];
parafile = [PATH  'parameter.txt'];
allspikes = dlmread(spikefile);

try
   parameter = dlmread(parafile,'\t',1,0);
catch
   parameter = dlmread(parafile);
end

cand_trials = parameter(:,6);
good_trials = find(cand_trials > 0);

trialstarttime = parameter(:,4);
alldatas = convert(allspikes, trialstarttime);

alldata = [];
for ii = 1 : length(good_trials)
    alldata(1, :, ii) = [alldatas(1+(good_trials(ii)-1)*5000: good_trials(ii)*5000   ) ]';
end

pre_spike_data = alldata(1,:,:);

tmpdt = alldata(1,:,:);
temp_spike_data =  reshape(tmpdt, 1, size(tmpdt,2)*size(tmpdt,3));

temp_azimuth = parameter(good_trials,3)';
temp_elevation = parameter(good_trials,2)';


data.spike_data = alldata;
%% CV star
             nostim_trials = squeeze(alldatas(:,:, temp_elevation < -99 & temp_azimuth < -99));
             isi = [];
             for i = 1:size(nostim_trials,2),
             cur_trial = nostim_trials(:,1);
                 while ~isempty(cur_trial),
                     idx = find(cur_trial > 0, 1, 'first');
                     cur_trial = cur_trial(idx+1:end);
                     isi = [isi idx];
                 end
                 m_cv_star(i) = cvtransfer(std(isi)/mean(isi), mean(isi));
                 
             end
             cv_star = cvtransfer(std(isi)/mean(isi), mean(isi));
             m_cv_star = median(m_cv_star);
%%





%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
null_trialsindex  = find(null_trials == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices

% bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response

bad_trials = find(cand_trials == 0);

if  length(bad_trials) > 1
    ss =0;
end


clear select_trials;
select_trials= ( (trials >= 1) & (trials <= length(good_trials)) ); 



for ii = 1 : floor((StopEventBin(1) - StartEventBin(1))/Step)+1;
   tempindex = StartEventBin(1) + (ii-1)*Step;
   for jj = 1 : length(null_trialsindex)
       m_NullResp(jj, ii) =  sum(alldata(SpikeChan,  tempindex - WindowInterval*0.5 : tempindex + WindowInterval*0.5, null_trialsindex(jj)));
   
   end
    
end



% find spontaneous trials which azimuth,elevation,stim_type=-9999
azimuth = temp_azimuth(~null_trials & select_trials);unique_azimuth = munique(azimuth');
elevation = temp_elevation(~null_trials & select_trials);unique_elevation = munique(elevation');

%     stim_type = temp_stim_type(~null_trials & select_trials);unique_stim_type = munique(stim_type');

stim_type = ones(1, length(azimuth));



spike_data = pre_spike_data(1,:,  ~null_trials & select_trials  );
% spike_rates= temp_spike_rates(~null_trials & select_trials);



condition_num = stim_type;unique_condition_num = munique(condition_num');
h_title{1}='Vestibular';h_title{2}='Visual';h_title{3}='Combined';
stim_duration = length(temp_spike_data)/length(temp_azimuth);


 
 unique_condition_num = 1;
 unique_stim_type = 1;
 
if Protocol == 100 | Protocol == 112
    trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_stim_type)+1;
    repetitions = floor( length(good_trials) / trials_per_rep); 
else if Protocol == 107
    trials_per_rep = length(unique_azimuth)* length(unique_stim_type) + 1;
    repetitions = floor( length(good_trials) / trials_per_rep); 
    end
end
unique_azimuth0=unique_azimuth;

%%  only for vestibular cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




 
%%

unique_azimuth0=[0:45:315]';
binsize=50;
for k=1: length(unique_condition_num)
    for j=1:length(unique_elevation)  
        clear spikeTracker;spikeTracker=zeros(repetitions,5000);
        for i=1: length(unique_azimuth0)
            clear select; select = logical( (azimuth==unique_azimuth0(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            act_found = find( select==1 );
%             clear spikeTracker;spikeTracker=zeros(repetitions,5000);
            if length(act_found)>repetitions
                NumRepetition=repetitions;            
            else
                NumRepetition=length(act_found);
            end
            for repeat=1:NumRepetition
                spikeTracker(repeat,:)= spike_data(1,stim_duration*(act_found(repeat)-1)+1:stim_duration*(act_found(repeat)));
            end
%             SpikeArray{k,j,i}(:,:)=spikeTracker;
            SpikeArray{k,j,i}=spikeTracker;
            SpikeHistArray{k,j,i}=mean(spikeTracker);   
            clear time_step; time_step=1;
            for n=1:(stim_duration/binsize)
                temp_PSTH(n)=mean(SpikeHistArray{k,j,i}(1,time_step:n*binsize))*1000;
                time_step=time_step+binsize;                
            end
            MaxCount{k}(j,i)=max(temp_PSTH);
        end
    end
    clear tempMaxRate;tempMaxRate=MaxCount{k};
    MaxRate(k)=max(max(tempMaxRate));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimating the peak or center of PSTH

% WindowInterval=400;%running window: ms (timebin for plot PSTH)
StartRespBin=StartEventBin(1,1)-0.5*WindowInterval;
DCWindow=900;

clear StepMatrix;
for k=1:length(unique_condition_num)
    tempDCResponse=[];
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth0)
            clear tempdata;tempdata=SpikeArray{k,j,i}(:,StartEventBin(1,1)-DCWindow+1:StartEventBin(1,1));
            tempDCResponse=[tempDCResponse;mean(tempdata')];            
            clear StartBinReset;StartBinReset=StartRespBin;
            clear TempPeakRate TempTroughRate;TempPeakRate=0.01; TempTroughRate=999;
            clear Index; Index=1;
            
            while StartBinReset<StopEventBin(1,1)-0.5*WindowInterval
                TempValue=mean(SpikeHistArray{k,j,i}(StartBinReset+1:StartBinReset+WindowInterval))*1000;
                if TempValue>TempPeakRate
                    TempPeakRate=TempValue;
                    PeakWindowStart=StartBinReset+1;
                    PeakWindowEnd=StartBinReset+WindowInterval;
                elseif TempValue<TempTroughRate
                    TempTroughRate=TempValue;
                    TroughWindowStart=StartBinReset+1;
                    TroughWindowEnd=StartBinReset+WindowInterval;
                end
                clear CurrentValue;CurrentValue=SpikeArray{k,j,i}(:,StartBinReset+1:StartBinReset+WindowInterval)*1000;
                clear CurrentValueMean; CurrentValueMean=mean(CurrentValue');
                StepMatrix{k,Index}(j,i,:)=CurrentValueMean;
                Index=Index+1;                
                StartBinReset=StartBinReset+Step;
            end                
            PeakRate(k,j,i)=TempPeakRate;PeakResp{k,j,i}=SpikeArray{k,j,i}(:,PeakWindowStart:PeakWindowEnd);%PeakResp{k,j,i}=SpikeHistArray{k,j,i}(PeakWindowStart:PeakWindowEnd);
            TroughRate(k,j,i)=TempTroughRate;TroughResp{k,j,i}=SpikeArray{k,j,i}(:,TroughWindowStart:TroughWindowEnd);%TroughResp{k,j,i}=SpikeHistArray{k,j,i}(TroughWindowStart:TroughWindowEnd);
            
            %try to find the PeakTimeIndex
            clear tempPeakTimeIndex;tempPeakTimeIndex=(0.5*(PeakWindowStart+PeakWindowEnd)-1000)*0.001;
            PeakTimeIndex(k,j,i)=tempPeakTimeIndex;  %PeakTimeIndex(k,j,i)=0.5*(PeakWindowStart+PeakWindowEnd); 
            PeakTimeStartIndex(k,j,i)=PeakWindowStart;PeakTimeEndIndex(k,j,i)=PeakWindowEnd;                            
            clear tempTroughTimeIndex;tempTroughTimeIndex=(0.5*(TroughWindowStart+TroughWindowEnd)-1000)*0.001;
            TroughTimeIndex(k,j,i)=tempTroughTimeIndex;   
            TroughTimeStartIndex(k,j,i)=TroughWindowStart;TroughTimeEndIndex(k,j,i)=TroughWindowEnd;
        end
    end
%     DCResponse{k}=tempDCResponse([1 length(unique_azimuth0)+1:length(unique_azimuth0)*4 length(unique_azimuth0)*4+1],:);
    if Protocol == 100
        DCResponse{k}=tempDCResponse([1 length(unique_azimuth0)+1:length(unique_azimuth0)*4 length(unique_azimuth0)*4+1],:);
        else if Protocol == 107 | Protocol == 112
          DCResponse{k}=tempDCResponse(:,:);
        end
    end
end


       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Significant test: Decide whether the response is significant or not 
for k=1:length(unique_condition_num)
    clear tempdata;tempdata=DCResponse{k};tempdata=reshape(tempdata,size(tempdata,1)*size(tempdata,2),1);
    DC=mean(tempdata);SD=std(tempdata);
    clear temp_p_peak temp_p_trough;
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth0)
            clear tempDC;tempDC=SpikeArray{k,j,i}(:,StartEventBin(1,1)-DCWindow+1:StartEventBin(1,1));       
            clear DCMean DCSD; DCMean=mean(tempDC');DCSD=std(tempDC');
%             clear tempDC2;tempDC2=SpikeArray{k,j,i}(:,StartEventBin(1,1)-200+1:StartEventBin(1,1)+200);       
            clear tempDC2;tempDC2=SpikeArray{k,j,i}(:,StartEventBin(1,1)-100+1:StartEventBin(1,1)+300);       
            
%             clear tempDC2;tempDC2=SpikeArray{k,j,i}(:,StartEventBin(1,1)+1:StartEventBin(1,1)+400);       
            clear DCMean2 DCSD2; DCMean2=mean(tempDC2');DCSD2=std(tempDC2');
            
            clear tempDC_Pre;tempDC_Pre=SpikeArray{k,j,i}(:,StartEventBin(1,1)-400+1:StartEventBin(1,1));       
            clear DCMean_Pre DCSD_Pre; DCMean_Pre=mean(tempDC_Pre');DCSD_Pre=std(tempDC_Pre');
            clear tempDC_Pre_1;tempDC_Pre_1=SpikeArray{k,j,i}(:,StartEventBin(1,1)-100+1:StartEventBin(1,1)+300);       
            clear DCMean_Pre_1 DCSD_Pre_1; DCMean_Pre_1=mean(tempDC_Pre_1');DCSD_Pre_1=std(tempDC_Pre_1');            
            clear tempDC_Pre_2;tempDC_Pre_2=SpikeArray{k,j,i}(:,StartEventBin(1,1)-100+1:StartEventBin(1,1)+100);       
            clear DCMean_Pre_2 DCSD_Pre_2; DCMean_Pre_2=mean(tempDC_Pre_2');DCSD_Pre_2=std(tempDC_Pre_2');
            
            clear tempDC_Aft;tempDC_Aft=SpikeArray{k,j,i}(:,StopEventBin(1,1)+1:StopEventBin(1,1)+400);                   
            clear DCMean_Aft DCSD_Aft; DCMean_Aft=mean(tempDC_Aft');DCSD_Aft=std(tempDC_Aft');            
            clear tempDC_Aft_1;tempDC_Aft_1=SpikeArray{k,j,i}(:,StopEventBin(1,1)+1:StopEventBin(1,1)+200);       
            clear DCMean_Aft_1 DCSD_Aft_1; DCMean_Aft_1=mean(tempDC_Aft_1');DCSD_Aft_1=std(tempDC_Aft_1');
            clear tempDC_Aft_2;tempDC_Aft_2=SpikeArray{k,j,i}(:,StopEventBin(1,1)+1:StopEventBin(1,1)+200);       
            clear DCMean_Aft_2 DCSD_Aft_2; DCMean_Aft_2=mean(tempDC_Aft_2');DCSD_Aft_2=std(tempDC_Aft_2');
            
            clear h_DC temp_p_DC;
            %[h_DC,temp_p_DC,ci] = ttest2(DCMean_Pre,DCMean_Aft,0.05);clear h_DC ci;            
            [temp_p_DC,h_DC] = ranksum(DCMean_Pre,DCMean_Aft,0.05);clear h_DC ci;            

            
            clear CurrentPeakResp CurrentRespMean CurrentRespSD; 
            CurrentPeakResp=PeakResp{k,j,i};CurrentRespMean=mean(CurrentPeakResp');%CurrentRespSD=std(CurrentPeakResp');            
            CurrentPeakResp1=SpikeArray{k,j,i}(:,PeakTimeStartIndex(k,j,i)-Step:PeakTimeEndIndex(k,j,i)-Step);
            CurrentRespMean1=mean(CurrentPeakResp1'); 
            CurrentPeakResp2=SpikeArray{k,j,i}(:,PeakTimeStartIndex(k,j,i)-2*Step:PeakTimeEndIndex(k,j,i)-2*Step);
            CurrentRespMean2=mean(CurrentPeakResp2');
            CurrentPeakResp3=SpikeArray{k,j,i}(:,PeakTimeStartIndex(k,j,i)+Step:PeakTimeEndIndex(k,j,i)+Step);
            CurrentRespMean3=mean(CurrentPeakResp3');
            CurrentPeakResp4=SpikeArray{k,j,i}(:,PeakTimeStartIndex(k,j,i)+2*Step:PeakTimeEndIndex(k,j,i)+2*Step);
            CurrentRespMean4=mean(CurrentPeakResp4');
            [temp_p_peak1(1),h_peak] = ranksum(CurrentRespMean,DCMean2,m_pval); clear h_peak;
            [temp_p_peak1(2),h_peak] = ranksum(CurrentRespMean1,DCMean2,m_pval); clear h_peak;
            [temp_p_peak1(3),h_peak] = ranksum(CurrentRespMean2,DCMean2,m_pval); clear h_peak;
            [temp_p_peak1(4),h_peak] = ranksum(CurrentRespMean3,DCMean2,m_pval); clear h_peak;
            [temp_p_peak1(5),h_peak] = ranksum(CurrentRespMean4,DCMean2,m_pval); clear h_peak;                                            
%             temp_p_peak(j,i)=max(max(temp_p_peak1), max(temp_p_peak2));
            temp_p_peak(j,i)=max(max(temp_p_peak1));            
            temp_PeakValue(j,i)=mean(CurrentRespMean)*1000;
            temp_PeakValueAll{j,i}=CurrentRespMean*1000;
            temp_PreStiValue(j,i)=mean(DCMean2)*1000;
            
            clear CurrentTroughResp CurrentRespMean CurrentRespSD; 
            CurrentTroughResp=TroughResp{k,j,i};CurrentRespMean=mean(CurrentTroughResp');CurrentRespSD=std(CurrentTroughResp');     
            CurrentTroughResp1=SpikeArray{k,j,i}(:,TroughTimeStartIndex(k,j,i)-Step:TroughTimeEndIndex(k,j,i)-Step);
            CurrentRespMean1=mean(CurrentTroughResp1');            
            CurrentTroughResp2=SpikeArray{k,j,i}(:,TroughTimeStartIndex(k,j,i)-2*Step:TroughTimeEndIndex(k,j,i)-2*Step);
            CurrentRespMean2=mean(CurrentTroughResp2');
            CurrentTroughResp3=SpikeArray{k,j,i}(:,TroughTimeStartIndex(k,j,i)+Step:TroughTimeEndIndex(k,j,i)+Step);
            CurrentRespMean3=mean(CurrentTroughResp3');
            CurrentTroughResp4=SpikeArray{k,j,i}(:,TroughTimeStartIndex(k,j,i)+2*Step:TroughTimeEndIndex(k,j,i)+2*Step);
            CurrentRespMean4=mean(CurrentTroughResp4');              
            
            [temp_p_trough1(1),h_trough] = ranksum(CurrentRespMean,DCMean2,m_pval);clear h_trough
            [temp_p_trough1(2),h_trough] = ranksum(CurrentRespMean,DCMean2,m_pval);clear h_trough
            [temp_p_trough1(3),h_trough] = ranksum(CurrentRespMean,DCMean2,m_pval);clear h_trough
            [temp_p_trough1(4),h_trough] = ranksum(CurrentRespMean,DCMean2,m_pval);clear h_trough
            [temp_p_trough1(5),h_trough] = ranksum(CurrentRespMean,DCMean2,m_pval);clear h_trough
            
            temp_p_trough(j,i)=max(max(temp_p_trough1));            
            temp_TroughValue(j,i)=mean(CurrentRespMean)*1000;
            temp_TroughValueAll{j,i}=CurrentRespMean*1000;
            
            temp_DCMean_Pre(j,i)=mean(DCMean_Pre)*1000;
            temp_DCMean_PreAll{j,i}=DCMean_Pre*1000;
            temp_DCMean_Pre_1(j,i)=mean(DCMean_Pre_1)*1000;
            temp_DCMean_PreAll_1{j,i}=DCMean_Pre_1*1000;            
            temp_DCMean_Pre_2(j,i)=mean(DCMean_Pre_2)*1000;
            temp_DCMean_PreAll_2{j,i}=DCMean_Pre_2*1000;
            
            temp_DCMean_Aft(j,i)=mean(DCMean_Aft)*1000; 
            temp_DCMean_AftAll{j,i}=DCMean_Aft*1000;             
            temp_DCMean_Aft_1(j,i)=mean(DCMean_Aft_1)*1000; 
            temp_DCMean_AftAll_1{j,i}=DCMean_Aft_1*1000; 
            
            temp_DCMean_p(j,i)=temp_p_DC;
        end
    end
    p_peak{k}=temp_p_peak;    
    p_trough{k}=temp_p_trough;
    PeakValue{k}=temp_PeakValue;
    PeakValueAll{k}=temp_PeakValueAll;
    TroughValue{k}=temp_TroughValue;
    TroughValueAll{k}=temp_TroughValueAll;
    PreStiValue{k}=temp_PreStiValue;
    DC_Pre{k}=temp_DCMean_Pre;
    DC_PreAll{k}=temp_DCMean_PreAll;
    DC_Pre_100300{k}=temp_DCMean_Pre_1;
    DC_PreAll_100300{k}=temp_DCMean_PreAll_1;
    DC_Pre_200{k}=temp_DCMean_Pre_2;
    DC_PreAll_200{k}=temp_DCMean_PreAll_2;   
    
    DC_Aft_400{k}=temp_DCMean_Aft;
    DC_AftAll_400{k}=temp_DCMean_AftAll;
    DC_Aft_200{k}=temp_DCMean_Aft_1;
    DC_AftAll_200{k}=temp_DCMean_AftAll_1;
    
    DC_p{k}=temp_DCMean_p;
end

%*******************************plot the figure*******************************%
% define figure
xoffset=0;yoffset=0;
x_time=[-1000:binsize:4000]*0.001;
x_time_smooth=([0:Step:5000-WindowInterval]-1000+0.5*WindowInterval)*0.001;
x_start = [(StartEventBin(1,1)-1000)*0.001, (StartEventBin(1,1)-1000)*0.001];
x_stop =  [(StopEventBin(1,1)-1000)*0.001,  (StopEventBin(1,1)-1000)*0.001];
y_marker=[0,max(MaxRate)];
Stimulus{1}='Vestibular';Stimulus{2}='Visual';Stimulus{3}='Combined';

for k=1:length(unique_condition_num)
    figure(k+1);set(k+1,'Position', [5,5 1250,950], 'Name', '3D Direction Tuning');%set(2,'Position', [5,5 1000,700], 'Name', '3D Direction Tuning');    
    orient landscape;
    if SpikeChan == 1
         text(-0.05,1.08,[FILE(1:end-4) 'Chan1'  ' // ' CurrentProtocol ' //' Stimulus{unique_condition_num(k)}]); axis off;    
    else
        text(-0.05,1.08,[FILE(1:end-4) 'Chan2' ' // ' CurrentProtocol ' //' Stimulus{unique_condition_num(k)}]); axis off;    
    end
   
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth0)
            %Compute PSTH
            clear time_step;time_step=1;
            for n=1:(stim_duration/binsize)
                temp_PSTH(n)=mean(SpikeHistArray{k,j,i}(1,time_step:n*binsize))*1000;
                time_step=time_step+binsize;
            end 
            figure(k+1);axes('position',[0.05*i+0.01+xoffset (0.98-0.05*j)+yoffset 0.045 0.045]);            
            bar(x_time(1:length(x_time)-1), temp_PSTH);xLim([min(x_time) max(x_time)]); yLim([0 MaxRate(k)]);                         
            hold on;plot( x_start, y_marker, 'r-','LineWidth',1.0);
            hold on;plot( x_stop,  y_marker, 'r-','LineWidth',1.0);
            if j==5
                if i==1
                    %                     xlabel(['Azimuth (deg): ' num2str(unique_azimuth0(i))]);
                    xlabel(['Azimuth: ' num2str(unique_azimuth0(i))]);
                else
                    xlabel(num2str(unique_azimuth0(i)));
                end
            else
                set(gca, 'xtick', [] ); 
            end 
            
            if i==1
                ylabel(num2str(unique_elevation(j)));
            else
                set(gca, 'ytick',[]);
            end
            
            clear StartBinReset;StartBinReset=0;
            clear nn;nn=1;
            clear TempValue
            while StartBinReset<5000-WindowInterval
                Temp_PSTH2(nn)=mean(SpikeHistArray{k,j,i}(StartBinReset+1:StartBinReset+WindowInterval))*1000;
                nn=nn+1;
                StartBinReset=StartBinReset+Step;
            end
            hold on;plot(x_time_smooth(1:length(x_time_smooth)-1),Temp_PSTH2,'b', 'LineWidth',2.0);            
            TempSti_PSTH2{j,i}=Temp_PSTH2;
            clear temp_p_peak; temp_p_peak=p_peak{k}(j,i);
            if temp_p_peak < m_pval             
                hold on;plot(PeakTimeIndex(k,j,i), PeakRate(k,j,i),'ro','LineWidth',2.0)
            end

            clear temp_p_trough;temp_p_trough=p_trough{k}(j,i);
            if temp_p_trough< m_pval               
                hold on;plot(TroughTimeIndex(k,j,i),TroughRate(k,j,i),'go','LineWidth',2.0);
            end
        end
    end
    PSTH_Smooth{k}=TempSti_PSTH2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PeakRateAmp PeakRateTime TroughRateAmp TroughRateTime;
if Protocol ~= 107
    
AziCir=0:45:360;
EleCir=-90:45:90;
[AziGrid,EleGrid]=meshgrid(AziCir,EleCir);      
IndexDir=[1 9:16 17:24 25:32 33];
Azi_Interp=0:5:360;
Ele_Interp=-90:5:90;
[Azi_InterpGrid,Ele_InterpGrid]=meshgrid(Azi_Interp,Ele_Interp);

for k=1: length(unique_condition_num)
    %One 3D tuning plot based on the peak value
    clear tempPeakRate; tempPeakRate=squeeze(PeakRate(k,:,:));
    tempPeakRate1=tempPeakRate;
    tempPeakRate1(:,9)=tempPeakRate(:,1);
    clear temp_p_peak; temp_p_peak=p_peak{k};
    clear Index; Index=find(temp_p_peak<  m_pval);  % ?????????
    PeakRateAmp{k}=tempPeakRate;    
    PeakRate_Interp=interp2(AziGrid,EleGrid,tempPeakRate1,Azi_InterpGrid,Ele_InterpGrid);
    %%%%%%%%%%%%%%%%%%%
    %2D to 1D transform (with spline interpolation)
    [maxRow_Interp,maxCol_Interp]=find(PeakRate_Interp==max(max(PeakRate_Interp)));    
    PeakAzi_Interp=Azi_InterpGrid(maxRow_Interp(1,1),maxCol_Interp(1,1));
    PeakEle_Interp=Ele_InterpGrid(maxRow_Interp(1,1),maxCol_Interp(1,1));    
    clear DiffResp DiffDir;
    NumCount=0;
    for m=1:size(PeakRate_Interp,1)
        for n=1:size(PeakRate_Interp,2)-1
            NumCount=NumCount+1;            
            DiffResp(NumCount)= PeakRate_Interp(m,n)/max(max(PeakRate_Interp));
            DiffDir(NumCount)=(180/pi) * acos( sin(PeakEle_Interp*pi/180).* sin(Ele_InterpGrid(m,n)*pi/180)  +  cos(Ele_InterpGrid(m,n)*pi/180).* sin(Azi_InterpGrid(m,n)*pi/180).* cos(PeakEle_Interp*pi/180).* sin(PeakAzi_Interp*pi/180) + cos(Ele_InterpGrid(m,n)*pi/180).* cos(Azi_InterpGrid(m,n)*pi/180).* cos(PeakEle_Interp*pi/180).* cos(PeakAzi_Interp*pi/180));            
        end
    end

    %Construct the histrogram data
    Azimuth_Bin=-7.5:15:180+7.5;
    clear Azimuth_Hist Resp_Hist;
    for m=2:length(Azimuth_Bin)
        clear Index;Index=find(DiffDir>Azimuth_Bin(m-1) & DiffDir<=Azimuth_Bin(m));
        Azimuth_Hist(m-1)=round(0.5*(Azimuth_Bin(m-1)+Azimuth_Bin(m)));
        Resp_Hist(m-1)=round(20*mean(DiffResp(Index)));       
    end
%     clear ModeTestData; ModeTestData=[];
%     for m=1:length(Azimuth_Hist)
%         ModeTestData=[ModeTestData Azimuth_Hist(m)*ones(1,Resp_Hist(m))];
%     end
%     p_Uniform{k}(1)=UniformTest_cah(ModeTestData',1);
%     [mode{k}(1,1), pval{k}(1,:)]=modalityTestForYong(ModeTestData',[1:3],1000); 
    
      %%%%%%%%%%%%%%%%%%%     
      %compute pAnova   
      clear temp_PeakValueAll0;temp_PeakValueAll0=reshape(PeakValueAll{k},size(PeakValueAll{k},1)*size(PeakValueAll{k},2),1);
      clear temp_PeakValueAll;temp_PeakValueAll=temp_PeakValueAll0(IndexDir);
      clear tempdata;
      for m=1:size(temp_PeakValueAll)    
          tempdata(m,:)=temp_PeakValueAll{m};
      end
      tempdata=tempdata';
      clear p_anova;[p_anova, table, stats]=anova1(tempdata,[],'off');
      P_anova_3D{k}(1,1) = p_anova;     
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      %Another plot based on the certain timing of the biggest response;
%       figure(k+1);axes('position',[0.38 0.41 0.23 0.21]);
      clear tempPeakRate;tempPeakRate=squeeze(PeakRate(k,:,:));
      %     clear temp_p_peak; temp_p_peak=p_peak{k};
%       clear Index; Index=find(temp_p_peak>0.01);
       clear Index; Index=find(temp_p_peak> m_pval);
      tempPeakRate(Index)=NaN;    
      [Index1,Index2] = find(tempPeakRate==max(max(tempPeakRate)));
      
      if length(Index1)>=1
        Index1=Index1(1,1);
        Index2=Index2(1,1);
        tempTimeIndexStart=PeakTimeStartIndex(k,Index1,Index2);
        tempTimeIndexEnd=PeakTimeEndIndex(k,Index1,Index2);
        %PeakTimeStartIndex(k,j,i)=PeakWindowStart;PeakTimeEndIndex(k,j,i)=PeakWindowEnd;
        for j=1:length(unique_elevation)
            for i=1:length(unique_azimuth0)
                clear CurrentPeakResp;CurrentPeakResp=mean(SpikeHistArray{k,j,i}(1,tempTimeIndexStart:tempTimeIndexEnd))*1000;                       
                CurrentPeakRate(j,i)=CurrentPeakResp;
                clear CurrentPeakRespAll;CurrentPeakRespAll=SpikeArray{k,j,i}(:,tempTimeIndexStart:tempTimeIndexEnd)*1000;                       
                RatePeakTimeAll{j,i}=mean(CurrentPeakRespAll');
            end
        end
        CurrentPeakRate1=CurrentPeakRate;CurrentPeakRate1(:,9)=CurrentPeakRate(:,1);                
        PeakRateTime{k}(:,:)=CurrentPeakRate;
        PeakRateTimeAll{k}=RatePeakTimeAll;
        clear PeakRate_Interp;PeakRate_Interp=interp2(AziGrid,EleGrid,CurrentPeakRate1,Azi_InterpGrid,Ele_InterpGrid);
    else 
        PeakRateTime{k}(:,:)=[];
        PeakRateTimeAll{k}=[];
    end    
    
    clear tempTroughRate; tempTroughRate=squeeze(TroughRate(k,:,:));    
    TroughRateAmp{k}=tempTroughRate;   
    clear tempTroughRate;tempTroughRate=squeeze(TroughRate(k,:,:));
    clear temp_p_trough; temp_p_trough=p_trough{k};
    clear Index; Index=find(temp_p_trough> m_pval);
    tempTroughRate(Index)=NaN;    
    clear Index1 Index2;[Index1,Index2] = find(tempTroughRate==min(min(tempTroughRate)));
    
    if length(Index1)>=1
        Index1=Index1(1,1);
        Index2=Index2(1,1);
        tempTimeIndexStart=TroughTimeStartIndex(k,Index1,Index2);
        tempTimeIndexEnd=TroughTimeEndIndex(k,Index1,Index2);   
        for j=1:length(unique_elevation)
            for i=1:length(unique_azimuth0)
                clear CurrentTroughResp;CurrentTroughResp=mean(SpikeHistArray{k,j,i}(1,tempTimeIndexStart:tempTimeIndexEnd))*1000;                       
                CurrentTroughRate(j,i)=CurrentTroughResp;
                clear CurrentPeakRespAll;CurrentTroughRespAll=SpikeArray{k,j,i}(:,tempTimeIndexStart:tempTimeIndexEnd)*1000;                       
                RateTroughTimeAll{j,i}=mean(CurrentTroughRespAll');
            end
        end
        TroughRateTime{k}(:,:)=CurrentTroughRate;
        TroughRateTimeAll{k}=RateTroughTimeAll;
    else
        TroughRateTime{k}(:,:)=[];
        TroughRateTimeAll{k}=[];
    end

    %%%%%%%%%%%%%%%%%%%
    %2D to 1D transform (with spline interpolation)
    [maxRow_Interp,maxCol_Interp]=find(PeakRate_Interp==max(max(PeakRate_Interp)));    
    PeakAzi_Interp=Azi_InterpGrid(maxRow_Interp(1,1),maxCol_Interp(1,1));
    PeakEle_Interp=Ele_InterpGrid(maxRow_Interp(1,1),maxCol_Interp(1,1));    
    clear DiffResp DiffDir;
    NumCount=0;
    for m=1:size(PeakRate_Interp,1)
        for n=1:size(PeakRate_Interp,2)-1
            NumCount=NumCount+1;            
            DiffResp(NumCount)= PeakRate_Interp(m,n)/max(max(PeakRate_Interp));
            DiffDir(NumCount)=(180/pi) * acos( sin(PeakEle_Interp*pi/180).* sin(Ele_InterpGrid(m,n)*pi/180)  +  cos(Ele_InterpGrid(m,n)*pi/180).* sin(Azi_InterpGrid(m,n)*pi/180).* cos(PeakEle_Interp*pi/180).* sin(PeakAzi_Interp*pi/180) + cos(Ele_InterpGrid(m,n)*pi/180).* cos(Azi_InterpGrid(m,n)*pi/180).* cos(PeakEle_Interp*pi/180).* cos(PeakAzi_Interp*pi/180));            
        end
    end

    %Construct the histrogram data
    Azimuth_Bin=-7.5:15:180+7.5;
    clear Azimuth_Hist Resp_Hist;
    for m=2:length(Azimuth_Bin)
        clear Index;Index=find(DiffDir>Azimuth_Bin(m-1) & DiffDir<=Azimuth_Bin(m));
        Azimuth_Hist(m-1)=round(0.5*(Azimuth_Bin(m-1)+Azimuth_Bin(m)));
        Resp_Hist(m-1)=round(20*mean(DiffResp(Index)));       
    end
%     clear ModeTestData; ModeTestData=[];
%     for m=1:length(Azimuth_Hist)
%         ModeTestData=[ModeTestData Azimuth_Hist(m)*ones(1,Resp_Hist(m))];
%     end
%     p_Uniform{k}(2,1)=UniformTest_cah(ModeTestData',1);
%     [mode{k}(2,1), pval{k}(2,:)]=modalityTestForYong(ModeTestData',[1:3],1000);    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute pAnova 
    if isempty(PeakRateTimeAll{k})        
        P_anova_3D{k}(1,2) = NaN;         
    else
        clear temp_PeakValueAll0;temp_PeakValueAll0=reshape(PeakRateTimeAll{k},size(PeakRateTimeAll{k},1)*size(PeakRateTimeAll{k},2),1);
        clear temp_PeakValueAll;temp_PeakValueAll=temp_PeakValueAll0(IndexDir);
        clear tempdata;
        for m=1:size(temp_PeakValueAll)    
            tempdata(m,:)=temp_PeakValueAll{m};
        end
        tempdata=tempdata';
        clear p_anova;[p_anova, table, stats]=anova1(tempdata,[],'off');
        P_anova_3D{k}(1,2) = p_anova; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot the another 3D tuning based on Peak-to-Trough
    clear tempPeakTroughRate; tempPeakTroughRate=squeeze(PeakRate(k,:,:))-squeeze(TroughRate(k,:,:));%    tempPeakRate(:,9)=tempPeakRate(:,1);
    PeakTroughRateAmp{k}=tempPeakTroughRate;    
    tempPeakTroughRate1=tempPeakTroughRate;tempPeakTroughRate1(:,9)=tempPeakTroughRate(:,1);
    clear PeakRate_Interp;PeakRate_Interp=interp2(AziGrid,EleGrid,tempPeakTroughRate1,Azi_InterpGrid,Ele_InterpGrid);
  
    %%%%%%%%%%%%%%%%%%%
    %2D to 1D transform (with spline interpolation)
    [maxRow_Interp,maxCol_Interp]=find(PeakRate_Interp==max(max(PeakRate_Interp)));    
    PeakAzi_Interp=Azi_InterpGrid(maxRow_Interp(1,1),maxCol_Interp(1,1));
    PeakEle_Interp=Ele_InterpGrid(maxRow_Interp(1,1),maxCol_Interp(1,1));    
    clear DiffResp DiffDir;
    NumCount=0;
    for m=1:size(PeakRate_Interp,1)
        for n=1:size(PeakRate_Interp,2)-1
            NumCount=NumCount+1;            
            DiffResp(NumCount)= PeakRate_Interp(m,n)/max(max(PeakRate_Interp));
            DiffDir(NumCount)=(180/pi) * acos( sin(PeakEle_Interp*pi/180).* sin(Ele_InterpGrid(m,n)*pi/180)  +  cos(Ele_InterpGrid(m,n)*pi/180).* sin(Azi_InterpGrid(m,n)*pi/180).* cos(PeakEle_Interp*pi/180).* sin(PeakAzi_Interp*pi/180) + cos(Ele_InterpGrid(m,n)*pi/180).* cos(Azi_InterpGrid(m,n)*pi/180).* cos(PeakEle_Interp*pi/180).* cos(PeakAzi_Interp*pi/180));            
        end
    end    
    %Construct the histrogram data
    Azimuth_Bin=-7.5:15:180+7.5;
    clear Azimuth_Hist Resp_Hist;
    for m=2:length(Azimuth_Bin)
        clear Index;Index=find(DiffDir>Azimuth_Bin(m-1) & DiffDir<=Azimuth_Bin(m));
        Azimuth_Hist(m-1)=round(0.5*(Azimuth_Bin(m-1)+Azimuth_Bin(m)));
        Resp_Hist(m-1)=round(20*mean(DiffResp(Index)));       
    end
%     clear ModeTestData; ModeTestData=[];
%     for m=1:length(Azimuth_Hist)
%         ModeTestData=[ModeTestData Azimuth_Hist(m)*ones(1,Resp_Hist(m))];
%     end
%     p_Uniform{k}(3,1)=UniformTest_cah(ModeTestData',1);
%     [mode{k}(3,1), pval{k}(3,:)]=modalityTestForYong(ModeTestData',[1:3],1000);
    %%%%%%%%%%%%%%%%%%%
    %compute pAnova   
    clear temp_PeakValueAll0;temp_PeakValueAll0=reshape(PeakValueAll{k},size(PeakValueAll{k},1)*size(PeakValueAll{k},2),1);
    clear temp_PeakValueAll;temp_PeakValueAll=temp_PeakValueAll0(IndexDir);    
    clear temp_TroughValueAll0;temp_TroughValueAll0=reshape(TroughValueAll{k},size(TroughValueAll{k},1)*size(TroughValueAll{k},2),1);
    clear temp_TroughValueAll;temp_TroughValueAll=temp_TroughValueAll0(IndexDir);
    
    clear tempdata;
    for m=1:size(temp_PeakValueAll)    
        tempdata_Peak(m,:)=temp_PeakValueAll{m};
        tempdata_Trough(m,:)=temp_TroughValueAll{m};
    end  
    tempdata=tempdata_Peak'-tempdata_Trough';
    clear p_anova;[p_anova, table, stats]=anova1(tempdata,[],'off');
    P_anova_3D{k}(1,3) = p_anova;        
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(unique_condition_num)   
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth0)
            clear StartBinReset;StartBinReset=StartRespBin;
            clear Index; Index=1;
            while StartBinReset<StopEventBin(1,1)-0.5*WindowInterval
                clear CurrentValue;CurrentValue=SpikeArray{k,j,i}(:,StartBinReset+1:StartBinReset+WindowInterval)*1000;
                clear CurrentValueMean; CurrentValueMean=mean(CurrentValue');
                StepMatrix{Index,unique_condition_num(k)}(j,i)=mean(CurrentValueMean);                
                clear CurrentDC;CurrentDC=DC_PreAll_100300{k}{j,i};
                clear temp_p;[temp_p,h_DC] = ranksum(CurrentValueMean,CurrentDC,0.05);clear h_DC ci;                
                Step_p{Index,unique_condition_num(k)}(j,i)=temp_p;                
                Index=Index+1;                
                StartBinReset=StartBinReset+Step;
            end                
        end
    end
end

TuningStep=MOOG_TuningStep_cah_simple_aff(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath,  WindowInterval,  Step);
TuningStep.Step_p=Step_p;
addpath('C:\work\poisson')
 spikes = getspikes1(TuningStep);
      
      for jj = 1 :  length(spikes)
          meanvar(jj, :) = [mean(spikes{jj}) var(spikes{jj})];
          pars(jj, :) = fminsearch('fitnegbinom',[10 1],[],spikes{jj});
      end 
  figure(2);
  axes('position',[0.02 (0.98-0.05*17)+0.23+yoffset 0.25 0.25]);
  

  plot(meanvar(:,1), meanvar(:,2), '.')
  hold on;
 plot( [min(min(meanvar))-2  max(max(meanvar))+2],  [min(min(meanvar))-2 max(max(meanvar))+2] ); 

 axes('position',[0.02 (0.98-0.05*17)-0.05 0.25 0.25]);
xx = [1+ pars(:,1).* pars(:,2)  meanvar(:,2)./meanvar(:,1) ];
plot(xx(:,1),  xx(:,2),  '.')
hold on;
 plot( [min(min(xx))-2  max(max(xx))+2],  [min(min(xx))-2  max(max(xx))+2] ); 
%  if  Protocol ~= 107
%    for j=1:length(unique_elevation)
%         for i=1:length(unique_azimuth0)
%             axes('position',[0.30+0.06*i 0.65-0.06*j 0.045 0.045]);   
%               
%              if  i == 1 && j == 1
%                    index_vm = 1; 
%                    plot(spikes(index_vm, :));
%              end
%              if  i == 1 && j == 5
%                    index_vm = 26; 
%                    plot(spikes(index_vm, :));
%              end
%              if j ~= 1 &&  j  ~= 5
%                    index_vm = 1 + (j-2)*8 + i;
%                    plot(spikes(index_vm, :))
%              end
%              set(gca, 'XTick', []);
%              set(gca, 'YTick', []);
% %              title([num2str( var(spikes(index_vm, :))/mean(spikes(index_vm, :)))   ] );
%         end
%    end
%    axes('position',[0.30 0.08 0.2 0.2]);  
%    plot(1:size(spikes,2), spikes(1:9, :));
%    axes('position',[0.52 0.08 0.2 0.2]);  
%    plot(1:size(spikes,2), spikes(1:18, :));
%    axes('position',[0.74 0.08 0.2 0.2]);  
%    plot(1:size(spikes,2), spikes(19:26, :));
%    
%  else
%         for i=1:length(unique_azimuth0)        
%             index_vm = i;
%              axes('position',[0.30+0.06*i 0.65-0.06*3 0.045 0.045]);   
%             plot(spikes(index_vm, :));
%             title([num2str( var(spikes(index_vm, :))/mean(spikes(index_vm, :)) )     ] );
%         end
%  end
                    
            
 
%%  add 4/17/09
% if SpikeChan == 1
%     SaveFileName=[OutputPath FILE(1:end-4) '_TuningStep_C1'];
%     save(SaveFileName,'TuningStep'); clear SaveFileName;
% end
% 
% if SpikeChan == 3
%     SaveFileName=[OutputPath FILE(1:end-4) '_TuningStep_C2'];
%     save(SaveFileName,'TuningStep'); clear SaveFileName;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output some values for further analysis
for k=1:length(unique_condition_num)    
    clear temptrial; 
    
    temptrial.name=[FILE(1:end-4) '_' Stimulus{unique_condition_num(k)}];
    temptrial.p_peak=p_peak{k};
    temptrial.p_trough=p_trough{k};
    temptrial.peaktimeindex=squeeze(PeakTimeIndex(k,:,:));
    temptrial.troughtimeindex=squeeze(TroughTimeIndex(k,:,:));
    temptrial.DC_p=DC_p{k};
    temptrial.PSTH_Smooth=PSTH_Smooth{k};
    temptrial.PeakValue=PeakValue{k};
    temptrial.PeakValueAll=PeakValueAll{k};
    temptrial.TroughValue=TroughValue{k};
    temptrial.TroughValueAll=TroughValueAll{k};
    temptrial.PreStiValue=PreStiValue{k};
    temptrial.DC_Pre=DC_Pre{k};
    temptrial.DC_PreAll=DC_PreAll{k};
    temptrial.DC_Pre_100_300=DC_Pre_100300{k};
    temptrial.DC_PreAll_100_300=DC_PreAll_100300{k};
    temptrial.DC_Pre_200=DC_Pre_200{k};
    temptrial.DC_PreAll_200=DC_PreAll_200{k};    
    temptrial.DC_Aft_400=DC_Aft_400{k};    
    temptrial.DC_AftAll_400=DC_AftAll_400{k};        
    temptrial.DC_Aft_200=DC_Aft_200{k};    
    temptrial.DC_AftAll_200=DC_AftAll_200{k};   

    if Protocol ~= 107
      temptrial.pAnova3D=P_anova_3D{k};
      temptrial.PeakRate_Time=PeakRateTime{k};
      temptrial.PeakRate_Amp=PeakRateAmp{k};
      temptrial.TroughRate_Amp=TroughRateAmp{k};
      temptrial.TroughRate_Time=TroughRateTime{k};
    end

%     temptrial.pUniform=p_Uniform{k};
%     temptrial.mode=mode{k};
%     temptrial.pval=pval{k}; 
   
        
%     SaveFileName=[OutputPath FILE(1:end-4) '_' Stimulus{unique_condition_num(k)}];
%     save(SaveFileName, 'temptrial');  %eval([temptrial.name '=temptrial']);
    
    FigureIndex=k+1;
    CurrentSti=Stimulus{unique_condition_num(k)};
%     CurrentSti=Stimulus{k};
    if Protocol == 100 || Protocol == 112
         m_peaktime  = ModeChecking_simple_3d(FILE(1:end-4), temptrial, TuningStep,CurrentProtocol, CurrentSti,FigureIndex, SpikeChan,  OutputPath, Step, m_pval);
    else
        if  Protocol == 107
             m_peaktime  = ModeChecking_simple_1d(FILE(1:end-4), temptrial, TuningStep,CurrentProtocol, CurrentSti,FigureIndex, SpikeChan,  OutputPath, Step, m_pval);
        end
    end
         TuningStep.peaktime=m_peaktime;
    

    
    figure(FigureIndex); 
    set(gcf, 'PaperOrientation', 'portrait');

    PeakTimeIndex =  m_peaktime;
    resp_VMtemp = TuningStep.StepResp_VM{PeakTimeIndex}'
    resp_VMtemp =  [resp_VMtemp;  mean(m_NullResp(:, PeakTimeIndex))  var(m_NullResp(:, PeakTimeIndex))];
    
  
    if  Protocol == 100
           SaveFileName=[OutputPath FILE(1:end-4) '_C'  num2str(SpikeChan)    '_3D'];
    else
        if  Protocol == 107 || Protocol == 112
%             SaveFileName=[OutputPath FILE(1:end-6) '_C'  num2str(SpikeChan)  '_' num2str(WindowInterval)  '_1D'];
            SaveFileName=[OutputPath FILE(1:end-4) '_C'  num2str(SpikeChan)   '_1D'];
        end
    end
    
   eval(['save ' SaveFileName ' PeakTimeIndex m_cv_star  cv_star TuningStep resp_VMtemp FILE m_NullResp']);


    
   
    
    
%     SaveFileName=[OutputPath FILE(1:end-4) '_C'  num2str(SpikeChan) '_'  num2str(WindowInterval)];
%     eval(['save ' SaveFileName ' TuningStep FILE m_slope m_intercept m_stat rr pp']);
 
    saveas(gcf,SaveFileName,'png');
    
    close(FigureIndex);     
    
end




return;
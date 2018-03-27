function Heading_PSTH1(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  FigureIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SysDelay=0;
StartEventBin=StartEventBin+SysDelay;
StopEventBin=StopEventBin+SysDelay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP


temp_rotamplitude = data.moog_params(ROT_AMPLITUDE,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_spike_data = data.spike_data(SpikeChan,:);   % spike rasters
 
for i=1:length(temp_heading)
    %[velmax maxAcc1s]=calpeakVel(1, 4, abs(temp_rotamplitude(i)));
    velmax = 2*abs(temp_rotamplitude(i))/(0.7+0.5);  % Trapzoid profile. Jing 09/06/2013
 
    temp_heading(i) = velmax *sin(pi*temp_heading(i)/180);
end
% % %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_heading);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type(select_trials );
heading = temp_heading(select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
% spike_rates = temp_spike_rates( select_trials);

% % total_trials = temp_total_trials( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');



% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_heading);
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 99;
end
spike_data = temp_spike_data( temp_spike_data~=99 );
spike_data( find(spike_data>90) ) = 1; % something is absolutely wrong 

%%%%%%%%%%%%%%%%%%%%%%%%neurometric dataset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the spikeArray
binsize=50;
Step=25;%Step=100;
WindowInterval=400; %timebin for plot PSTH
x_time=([-1000:binsize:4000]-115)*0.001;% x_time=[-1000:binsize:4000]*0.001;
SmoothBin_L=[StartEventBin(1,1):-Step:0.5*WindowInterval];SmoothBin_L=fliplr(SmoothBin_L);SmoothBin_L=SmoothBin_L(1:end-1);
SmoothBin_R=[StartEventBin(1,1):Step:5000-0.5*WindowInterval];
SmoothBin=[SmoothBin_L SmoothBin_R];
x_time_smooth=(SmoothBin-StartEventBin(1,1))*0.001;%x_time_smooth=([0:Step:5000-WindowInterval]-1000+0.5*WindowInterval)*0.001;

StartPoint=0;
clear tempdata; tempdata=SmoothBin-(StartEventBin(1,1)+StartPoint);
clear IndexStart; IndexStart=find(abs(tempdata)==min(abs(tempdata)));
clear tempdata; tempdata=SmoothBin-StopEventBin(1,1);
clear IndexEnd; IndexEnd=find(abs(tempdata)==min(abs(tempdata)));
 
trials_per_rep = length(unique_heading)*length(unique_stim_type);
repetition = floor( length(select_trials)/trials_per_rep); % take minimum repetition

resp_mat_anova = [];
for k=1:1 
    clear temp_p_peak temp_p_trough;
    for i=1:length(unique_heading)
        clear select; select = logical( (heading==unique_heading(i)) & (stim_type==unique_stim_type(k)) );  
        act_found = find( select==1 );
        for repeat=1:length(act_found)
            spikeTracker(repeat,:)=double(spike_data(1,stim_duration*(act_found(repeat)-1)+1:stim_duration*(act_found(repeat))));
        end
        SpikeArray{k,i}=spikeTracker;
        SpikeHistArray{k,i}=mean(spikeTracker);  
        
        %Get the spontaneous response      
        clear temp_Spont; temp_Spont=mean(SpikeHistArray{k,i}(StartEventBin(1,1)-0.5*WindowInterval: StartEventBin(1,1)+0.5*WindowInterval))*1000;
        SpontResp(k,i)=temp_Spont;
        
        %Get the individual PSTH        
        clear time_step;time_step=1;
        for n=1:(stim_duration/binsize)
            temp_PSTH(n)=mean(SpikeHistArray{k,i}(1,time_step:n*binsize))*1000;
            time_step=time_step+binsize;
        end        
        PSTH_Ori{k,i}=temp_PSTH;
        
        for nn=1:length(SmoothBin)            
            temp_PSTH2(nn)=mean(SpikeHistArray{k,i}(SmoothBin(nn)-0.5*WindowInterval: SmoothBin(nn)+0.5*WindowInterval))*1000;
            clear tempdata1; tempdata1=spikeTracker(:,SmoothBin(nn)-0.5*WindowInterval: SmoothBin(nn)+0.5*WindowInterval);
            clear tempdata2; tempdata2=sum(tempdata1')*1000/WindowInterval;tempdata2=tempdata2';
%             resp_mat_anova{k,nn}(1:repetition,i)=tempdata2(1:repetition); clear tempdata1 tempdata2;
        end      
        
        clear tempdata; tempdata=temp_PSTH2(IndexStart:IndexEnd);
        clear currentPeakIndex; currentPeakIndex=find(tempdata==max(tempdata));
        clear PeakIndex; PeakIndex=currentPeakIndex(1,1)+IndexStart-1;
        
        tempDC=spikeTracker(:,SmoothBin(IndexStart)-0.5*WindowInterval: SmoothBin(IndexStart)+0.5*WindowInterval);
        clear DCMean; DCMean=mean(tempDC');
        
        CurrentPeakResp=spikeTracker(:,SmoothBin(PeakIndex)-0.5*WindowInterval: SmoothBin(PeakIndex)+0.5*WindowInterval);
        CurrentRespMean=mean(CurrentPeakResp');
        CurrentPeakResp1=spikeTracker(:,SmoothBin(PeakIndex-1)-0.5*WindowInterval: SmoothBin(PeakIndex-1)+0.5*WindowInterval);
        CurrentRespMean1=mean(CurrentPeakResp1');
        CurrentPeakResp2=spikeTracker(:,SmoothBin(PeakIndex-2)-0.5*WindowInterval: SmoothBin(PeakIndex-2)+0.5*WindowInterval);
        CurrentRespMean2=mean(CurrentPeakResp2');        
        CurrentPeakResp3=spikeTracker(:,SmoothBin(PeakIndex+1)-0.5*WindowInterval: SmoothBin(PeakIndex+1)+0.5*WindowInterval);
        CurrentRespMean3=mean(CurrentPeakResp3');
        CurrentPeakResp4=spikeTracker(:,SmoothBin(PeakIndex+2)-0.5*WindowInterval: SmoothBin(PeakIndex+2)+0.5*WindowInterval);
        CurrentRespMean4=mean(CurrentPeakResp4');        
        [temp_p_peak1(1),h_peak] = ranksum(CurrentRespMean,DCMean,0.01); clear h_peak;
        [temp_p_peak1(2),h_peak] = ranksum(CurrentRespMean1,DCMean,0.01); clear h_peak;
        [temp_p_peak1(3),h_peak] = ranksum(CurrentRespMean2,DCMean,0.01); clear h_peak;
        [temp_p_peak1(4),h_peak] = ranksum(CurrentRespMean3,DCMean,0.01); clear h_peak;
        [temp_p_peak1(5),h_peak] = ranksum(CurrentRespMean4,DCMean,0.01); clear h_peak;                                            
        p_peak(k,i)=max(temp_p_peak1);            
        PeakRate(k,i)=mean(CurrentRespMean)*1000;
        
        clear currentTroughIndex; currentTroughIndex=find(tempdata==min(tempdata));
        clear TroughIndex; TroughIndex=currentTroughIndex(1,1)+IndexStart-1;
        CurrentTroughResp=spikeTracker(:,SmoothBin(TroughIndex)-0.5*WindowInterval: SmoothBin(TroughIndex)+0.5*WindowInterval);
        CurrentRespMean=mean(CurrentTroughResp');                
        CurrentTroughResp1=spikeTracker(:,SmoothBin(TroughIndex-1)-0.5*WindowInterval: SmoothBin(TroughIndex-1)+0.5*WindowInterval);
        CurrentRespMean1=mean(CurrentTroughResp1');        
        CurrentTroughResp2=spikeTracker(:,SmoothBin(TroughIndex-2)-0.5*WindowInterval: SmoothBin(TroughIndex-2)+0.5*WindowInterval);
        CurrentRespMean2=mean(CurrentTroughResp2');                 
        CurrentTroughResp3=spikeTracker(:,SmoothBin(TroughIndex+1)-0.5*WindowInterval: SmoothBin(TroughIndex+1)+0.5*WindowInterval);
        CurrentRespMean3=mean(CurrentTroughResp3');        
        CurrentTroughResp4=spikeTracker(:,SmoothBin(TroughIndex+2)-0.5*WindowInterval: SmoothBin(TroughIndex+2)+0.5*WindowInterval);
        CurrentRespMean4=mean(CurrentTroughResp4');
        [temp_p_trough1(1),h_trough] = ranksum(CurrentRespMean,DCMean,0.01);clear h_trough
        [temp_p_trough1(2),h_trough] = ranksum(CurrentRespMean,DCMean,0.01);clear h_trough
        [temp_p_trough1(3),h_trough] = ranksum(CurrentRespMean,DCMean,0.01);clear h_trough
        [temp_p_trough1(4),h_trough] = ranksum(CurrentRespMean,DCMean,0.01);clear h_trough
        [temp_p_trough1(5),h_trough] = ranksum(CurrentRespMean,DCMean,0.01);clear h_trough        
        p_trough(k,i)=max(max(temp_p_trough1));            
        TroughRate(k,i)=mean(CurrentRespMean)*1000;          
        Time_Peak(k,i)=x_time_smooth(PeakIndex);
        Time_Trough(k,i)=x_time_smooth(TroughIndex);             
        PSTH_Smooth{k,i}=temp_PSTH2;        
        MaxCount(k,i)=max(temp_PSTH);         
    end          
    clear tempMaxRate;tempMaxRate=MaxCount(k,:);
    MaxRate(k)=max(max(tempMaxRate));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the PSTH
x_start = [(StartEventBin(1,1)-1000)*0.001, (StartEventBin(1,1)-1000)*0.001];
x_stop =  [(StopEventBin(1,1)-1000)*0.001,  (StopEventBin(1,1)-1000)*0.001];

% FigureIndex=2;
% figure(FigureIndex);set(FigureIndex,'Position', [5,50 1200,750], 'Name', 'HeadingDiscrimination');    
% orient landscape;text(-0.05,1.08,FILE); axis off;

% h_title{1}='Vest';h_title{2}='Visu';h_title{3}='Comb';

for k=1:1  
    for i=1:length(unique_heading)
        figure(FigureIndex); 
        axes('position',[0.38 (0.6-0.06*i) 0.11 0.04]);%figure(FigureIndex); axes('position',[0.14*(k-1)+0.3 (0.78-0.076*i) 0.12 0.06]);  
        clear temp_PSTH; temp_PSTH=PSTH_Ori{k,i};        
%         bar(x_time(1:length(x_time)-1), temp_PSTH);
        area(x_time(1:length(x_time)-1), temp_PSTH)
        xLim([min(x_time) max(x_time)]); yLim([MaxRate(k)-MaxRate(k)/3 MaxRate(k)+MaxRate(k)/8]);           
        clear temp_PSTH2;temp_PSTH2=PSTH_Smooth{k,i};
%         hold on;plot(x_time_smooth,temp_PSTH2,'b', 'LineWidth',2.0);  
%         y_marker=[0, MaxRate(k)];         
%         hold on;plot( x_start, y_marker, 'k-','LineWidth',2.0);
%         hold on;plot( x_stop,  y_marker, 'k-','LineWidth',2.0);
        
        %Add temporal modulation information
%         clear temp_p_peak; temp_p_peak=p_peak(k,i);
%         if temp_p_peak <0.01             
%             hold on;plot(Time_Peak(k,i), PeakRate(k,i),'ro','LineWidth',2.0)
%         end
%         clear temp_p_trough;temp_p_trough=p_trough(k,i);
%         if temp_p_trough<0.01                
%             hold on;plot(Time_Trough(k,i),TroughRate(k,i),'go','LineWidth',2.0);
%         end
        
%         %Add peaktime infomration
%         clear tempPeakTime; 
%         tempPeakTime=[tempmaker(1), tempmaker(1)];
%         hold on;
%         plot(tempPeakTime,y_marker,'r--','LineWidth',1.0);
%         clear tempPeakTime; 
%         tempPeakTime=[tempmaker(2), tempmaker(2)];
%         hold on;
%         plot(tempPeakTime,y_marker,'b--','LineWidth',1.0);
        
%         if i==1
%             title(h_title{k});
%         else            
%         end
        xlim([-0.25 1])
        if i==length(unique_heading)   
           set(gca, 'xtick', [0 0.7] );  
            xlabel('Time (s)');
        else
            set(gca, 'xtick', [] );            
        end        
        if k==1
            ylabel(num2str(unique_heading(i),  '%.1f'));
        else
            set(gca, 'ytick',[]);
        end           
    end    
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output some values for further analysis
clear HeadingPSTH;
HeadingPSTH.name=FILE(1:end-4); 
HeadingPSTH.Coherence=unique_motion_coherence;  
HeadingPSTH.Repets=repetition; 
HeadingPSTH.SpontResp=mean(SpontResp);
HeadingPSTH.PSTHOri=PSTH_Ori;
HeadingPSTH.PSTHSmooth=PSTH_Smooth;

HeadingPSTH.p_peak=p_peak;
HeadingPSTH.TimePeak=Time_Peak;
HeadingPSTH.PeakValue=PeakRate;
HeadingPSTH.p_trough=p_trough;
HeadingPSTH.TimeTrough=Time_Trough;
HeadingPSTH.TroughValue=TroughRate;

% OutputPath=['C:\Aihua\z_TempOutputs\'];
% SaveFileName=[OutputPath FILE(1:end-4) '_HeadingPSTH'];
% save(SaveFileName, 'HeadingPSTH');

% %save the figure
% figure(FigureIndex); 
% set(gcf, 'PaperOrientation', 'portrait');
% saveas(gcf,[OutputPath FILE(1:end-4) '_Heading.png'],'png');
% close(FigureIndex)

return;



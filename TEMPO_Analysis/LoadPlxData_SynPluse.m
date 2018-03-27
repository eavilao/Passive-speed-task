function [good_data outputChan]=LoadPlxData_SynPluse(good_data,FName,goodTrialInd)

load(FName);

if(isempty(data))
    return;
end

%****************************************
% Find all the good trials index.
ts = data.marker.timestamps;
len=length(ts);
j=1;
k=1;
flag4 = 0;
flag5 = 0;
for i= 1:len-1
    dt = ts(i+1)-ts(i);
    if dt < 0.5 %0.017
        if flag4 == 0 
            mkt4(k)= ts(i);
            flag4 = 1;
            flag5 = 0;
            j=j+1;
        end;
    else
        if flag5 == 0
            mkt5(k)= ts(i);
            flag5 = 1;
            flag4 = 0;
            k=k+1;
        end;        
    end;
end;
%mkt5(k)= ts(len);
%dt45=mkt5-mkt4
%find(dt45>0.95)   for YF
%find(dt45>0.65);   %for Fine/Coarse Yaw

disp(['evt2=' num2str(k)]);
disp(['goodTrialInd=', num2str(goodTrialInd')]);

goodTrialNum = length(goodTrialInd)
for i=1:goodTrialNum
    startCT(i)=mkt4(goodTrialInd(i))-0.0205; % offset 20.5ms is the median time from marker04 to syncpulse
    startT(i)=startCT(i)-1;
    endT(i)=startCT(i)+4;
    %eventCode{i,:}=[4 5];
end;
%****************************************

for i=1:16
    for j=1:4
        if data.spikes.counts(i,j+1) ~= 0
            spsData{i,j}.spsUnit=[i,j];
            indL(i,j)=0;
            for k=1:goodTrialNum
                %spsData{i,j}.spikeInfo(k).eventCode=eventCode{k,:};
                
                spsData{i,j}.spikeInfo(k).startCodeTime = startCT(k);
                spsData{i,j}.spikeInfo(k).startTime = startT(k);
                spsData{i,j}.spikeInfo(k).endTime =  endT(k);
        
                spts=data.spikes.timestamps{i,j+1};
                sptsInd=find(spts>=startT(k) & spts<=endT(k));
                spsData{i,j}.spikeInfo(k).SpikeTimes=spts(sptsInd);
                indL(i,j)=indL(i,j)+length(sptsInd);
            end
            spsData{i,j}.spsNum =indL(i,j);
        else
           spsData{i,j}=[];
        end;        
    end;
end;

outputChan=5;
for i=1:16
    for j=1:4
        if ~isempty(spsData{i,j})
            outputChan=outputChan+1;
            good_data.spike_unit(outputChan,:)=spsData{i,j}.spsUnit;
            for k=1:goodTrialNum
                startT=round(1000*spsData{i,j}.spikeInfo(k).startTime);
                endT=round(1000*spsData{i,j}.spikeInfo(k).endTime);
                
                edges=startT:1:endT;
                sps=spsData{i,j}.spikeInfo(k).SpikeTimes;
                if isempty(sps)
                    n=zeros(1,5000);
                else
                    n=histc(1000*spsData{i,j}.spikeInfo(k).SpikeTimes,edges);
                    n=n(1:length(n)-1);
                end;
                %good_data.spike_data(outputChan,:,k)=n;
                good_data.spike_data(outputChan,:,k)=n(1:5000);
            end;
        end;
    end;
end;






        
        
        
        
        
        


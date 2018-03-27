function [good_data outputChan]=LoadPlxData(good_data,FName)

load(FName);

if(isempty(data))
    return;
end

% Find all the good trials index.
j=0;
i=1;
goodTrialInd = [];
while (i <= data.marker.counts)
    % If we find a stimuls start code, then we look to see if we find
    % a success code before we get to another start code.
    if(data.marker.values(i) == 5)
        j = i + 1;
        if(data.marker.values(i+1) == 5)
            j = i + 2;
        end;
        
        while (j <= data.marker.counts & data.marker.values(j) ~= 5)
            if data.marker.values(j) == 12
                goodTrialInd = [goodTrialInd, i];
                break;
            else
                j = j + 1;
            end
        end
        i = j;
    else
        i=i+1;
    end  
end

goodTrialInd=goodTrialInd-1;

goodTrialNum = length(goodTrialInd)
for i=1:goodTrialNum
    startCT(i)=data.marker.timestamps(goodTrialInd(i));
    startT(i)=startCT(i)-1;
    endT(i)=startCT(i)+4;
    mkt=data.marker.timestamps;
    mkInd=find(mkt>=startT(i) & mkt<=endT(i));
    eventCode{i,:}=data.marker.values(mkInd);
end;


for i=1:16
    for j=1:4
        if data.spikes.counts(i,j+1) ~= 0
            spsData{i,j}.spsUnit=[i,j];
            indL(i,j)=0;
            for k=1:goodTrialNum
                spsData{i,j}.spikeInfo(k).eventCode=eventCode{k,:};
                
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
                good_data.spike_data(outputChan,:,k)=n;
            end;
        end;
    end;
end;






        
        
        
        
        
        


function data=convert(spike,starttime)

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
function data=convert_1(spike,starttime)
%  it is designed for rotatin protocols,  From 0.1 to 0.6 s
%%%%this m code convert the spike train to the format of tempo, the nth
%%%%trial 's start time is the starttime(n);
trialnum=length(starttime);
clear data
edges=((0:1:5000)-995)/1000;
minbin=min(edges);
maxbin=max(edges);
for i=1:trialnum
    
    spikenow=spike-starttime(i);
    spikenowtouse=spikenow((spikenow>=0.115+0.1)&(spikenow<=0.115+0.6));
    clear spikenow
    if(~isempty(spikenowtouse))  
       data{1,i} = spikenowtouse;
    else
       data{1,i}=-1;
    end
clear spikenowtouse
clear currentconvert
end
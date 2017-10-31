function [t,spktrain] = wf2spktrain(spkwf,tspk,fs); 

% spkwf = units{1}.wfspk; % waveforms
% tspk = units{1}.tspk; % timing of spikes
% spkwf2 = units{2}.wfspk; % waveforms
% tspk2 = units{2}.tspk; % timing of spikes


% fs= 40000; % Sampling rate



nsamp = (max(tspk)+1)*fs; % seconds to samples
spklength = size(spkwf,2); % length of wavef samples
spktrain = nan(1,nsamp); % Create a vector of NaNs of the lenght of nsamp
for i=1:size(spkwf,1)
    spktrain(round(tspk(i)*fs):round(tspk(i)*fs)+spklength-1) = spkwf(i,:);
end
t = linspace(0,max(tspk)+1,nsamp);



% Plot
%plot(t2,spktrain2, 'Color', 'k')
%hold on
%plot(t1,spktrain1, 'Color', 'c')


% Plot
%plot(t2(1:22:end),spktrain2(1:22:end), 'Color', 'c'); hold on
%plot(t1(1:22:end),spktrain1(1:22:end), 'Color', 'k'); 
% box off
% set(gca, 'TickDir', 'out'); 

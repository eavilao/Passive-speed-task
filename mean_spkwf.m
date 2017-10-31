function spkwf2 = mean_spkwf(spkwf)

% invert waveform if needed
if abs(max(mean(spkwf))) > abs(min(mean(spkwf))), spkwf = -spkwf; end

% remove outliers and align troughs
for i=1:size(spkwf,1)
    minsamp(i) = find(spkwf(i,:)==min(spkwf(i,:)),1);
end
minsamp_bounds = prctile(minsamp,[10 90]);
spkwf(minsamp<minsamp_bounds(1) | minsamp>minsamp_bounds(2),:)=NaN;
spkwf = [NaN(size(spkwf,1),10) spkwf NaN(size(spkwf,1),10)];
for i=1:size(spkwf,1)
    spkwf(i,:) = circshift(spkwf(i,:),[0 9-minsamp(i)]);
end
spkwf = spkwf(:,11:end-10);

% copmute mean and sig
spkwf2.mu = nanmean(spkwf);
spkwf2.sig = nanstd(spkwf);

% compute spike half-width (in microseconds)
minsamp = find(spkwf2.mu == min(spkwf2.mu));
maxsamp = minsamp + find(spkwf2.mu(minsamp:end) == max(spkwf2.mu(minsamp:end)));
spkwf2.width = (maxsamp - minsamp)*(1/40000)*1e6; % 40k samp/s
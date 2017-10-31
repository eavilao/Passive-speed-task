function [spktrain_t,spktrain_r] = spktimes2spktrain(spktimes,tbeg,tend,binwidth)

spktrain_t = linspace(tbeg,tend,(tend-tbeg)/binwidth);
spktrain_r = zeros(1,(tend-tbeg)/binwidth);
spktimes = spktimes(spktimes>=0 & spktimes<1.5);
spkbins = ceil(spktimes/binwidth);
for i=1:length(spktrain_r)
    spktrain_r(i) = sum(ceil(spktimes/binwidth)==i);
end
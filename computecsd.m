function csd = computecsd(lfps,t,dx)    
% dx in mm
nchans = size(lfps,1);
% baseline subtracted lfps
for i=1:nchans
    lfps(i,:) = lfps(i,:) - mean(lfps(i,t>0 & t<0.4));
end
for i=2:nchans-1
    csd(i-1,:) = (lfps(i-1,:) - 2*lfps(i,:) + lfps(i+1,:))/(dx^2);
end
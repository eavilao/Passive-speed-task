function xr = spiketrains_xcorr(tspk1,tspk2)
tmax = max([tspk1 ; tspk2]);
ts1 = zeros(1,round(tmax+1e3)*1e3); ts2 = ts1;
ts1(1+round(tspk1*1e3)) = 1;
ts2(1+round(tspk2*1e3)) = 1;
maxlag=200; xr = xcorr(ts1,ts2,maxlag); % maxlag=200ms
% xr.t = (-maxlag:1:maxlag)/1e3;
% xr(maxlag+1)=0;
% plot(-maxlag:1:maxlag,xr);
% xlabel('lag (ms)'); ylabel('spike counts');
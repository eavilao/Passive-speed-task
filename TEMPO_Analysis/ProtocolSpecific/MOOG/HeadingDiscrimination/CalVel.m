function [return_vel return_acc]=CalVel(dur, sig, mag, windows)


vel = GenGaussian(dur, sig);
dist = cumtrapz(vel)/600;
for i = 1 : length(mag)
    dist = dist/max(dist)*mag(i);
    vel = diff(dist)*600;    
    accel= diff(vel)*600;   
end
velmax = abs(max(vel)) ;       %%
maxAcc1s = abs(max(accel))*100;%% cm/s/s

midpt =  floor(length(vel)/2);
window = floor(length(vel)*windows/dur);
return_vel = mean( vel(midpt -  floor(window/2) : midpt + floor(window/2) )  );
return_acc = mean( accel(midpt -  floor(window/2) : midpt + floor(window/2) )  );

function [data] = GenGaussian(dur, sig)
t0 = dur/2;
t = 0:1/600:dur;

% Generate the Gaussian.
data = exp(-sqrt(2)* ((t-t0) / (dur/sig)).^2);
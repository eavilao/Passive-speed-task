function [r2,t2] = smooth_colormap_1DAzi(r,t,cond)

% smooth psth
twin = 0.1;
dt = diff(t); dt = dt(1);
nwin = round(twin/dt);
nt = length(t);
for i=1:nt-nwin+1
    r2(:,i) = mean(r(:,i:i+nwin-1),2);
    t2(i) = mean(t(i:i+nwin-1));
end

nunits = size(r2,1);
% normalise psth
for i=1:nunits
%     r2(i,:) = r2(i,:) - mean(r2(i,t2>0 & t2<0.4));
r2(i,:) = r2(i,:) - mean(r2(i,:));
    r2(i,:) = r2(i,:)/max(r2(i,:));
end

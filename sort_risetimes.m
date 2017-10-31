function [r2,t2,t_on2,t_off2] = sort_risetimes(r,t,t_on,t_off,cond)

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
    r2(i,:) = r2(i,:) - mean(r2(i,t2>0 & t2<0.4));
    r2(i,:) = r2(i,:)/max(r2(i,:));
end

% sustained neurons
% switch cond
%     case 'ves'
%         indx_sus = [1 3 19 22 23 24 55 66 68];
%         t_off(indx_sus) = [1.135 1.195 1.165 1.135 1.165 1.195 1.115 1.205 1.235];
%     case 'vis'
%         indx_sus = [3 10 16 19 27 28 34];
%         t_on(indx_sus) = [0.5613 0.5312 0.5211 0.5513 0.4708 0.5513 0.5714];
%         indx = randperm(nunits);
%         t_on = t_on(indx); r2 = r2(indx,:);
%         t_off(~isnan(t_off))=nan;
%     case 'com'
%         indx_sus = [2 17 20 36 54 58 68];
%         t_off(indx_sus) = [1.185 1.205 1.074 1.145 1.145 1.155 nan];
% end

% sort neurons
off_units = ~isnan(t_off);
on_units = isnan(t_off);
r2 = [r2(off_units,:) ; r2(on_units,:)]; 
t_on2 = [t_on(off_units) t_on(on_units)];
t_off2 = [t_off(off_units) t_off(on_units)];
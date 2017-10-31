function [x2,y2] = groupmean(x,y,groupsize)
if nargin<3, groupsize=1; end
if ~isempty(x)
    x2 = unique(x);
    x2 = reshape(x2,[groupsize length(x2)/groupsize])';
    for i=1:size(x2,1)
        y2.mu(i) = nanmean(y(ismember(x,x2(i,:))));
        y2.std(i) = nanstd(y(ismember(x,x2(i,:))))/sqrt(sum(ismember(x,x2(i,:))));
    end
    x2 = mean(x2,2);
else
    x2=[]; y2=[];
end
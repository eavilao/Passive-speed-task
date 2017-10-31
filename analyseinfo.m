function fisherinfo = analyseinfo(stim,nspk,modality)

%%
modalities = stim.modality;
switch modality
    case 'ves_nofix'
        indx = modalities==0;
    case 'ves'
        indx = modalities==1;
    case 'vis'
        indx = modalities==2;
    case 'com'
        indx = modalities==3;
end
stim = stim.speed(indx);
stims = unique(stim);
nspk = nspk(:,indx);
nunits = size(nspk,1);

%% all speeds
% compute signal: dfds
indx1 = stim==stims(1)|stim==stims(2);
indx2 = stim==stims(end)|stim==stims(end-1);
s1 = mean(stim(indx1)); s2 = mean(stim(indx2));
f1 = mean(nspk(:,indx1)'); f2 = mean(nspk(:,indx2)');
dfds = (f2-f1)/(s2-s1);
% remove signal from nspk to compute noise
nspk(:,indx1) = nspk(:,indx1) - repmat(f1',[1 sum(indx1)]);
nspk(:,indx2) = nspk(:,indx2) - repmat(f2',[1 sum(indx2)]);
for k=1:nunits
    clear J;
    these_units = nchoosek(1:nunits,k);
    if size(these_units,1)>2000
        these_units = these_units(randperm(size(these_units,1)),:);
        these_units = these_units(1:2000,:);
    end
    for j=1:size(these_units,1)
        these_slopes = dfds(these_units(j,:));
        these_nspk = nspk(these_units(j,:),:);
        [c,~,l]=pca(these_nspk');
        J(j,1) = ((these_slopes*c(:,1))^2)/l(1); % fisher info in the first component
        if k>=2
            J(j,2) = ((these_slopes*c(:,2))^2)/l(2); % fisher info in the second component
        else
            J(j,2) = 0;
        end
    end
    fisherinfo.mean(k,:) = mean(J,1);
    fisherinfo.raw(k).J = J;
end

%% group low and high speeds
% if any(sign(stims)==-1)
%     stims2(1,:) = stims(stims<0); stims2(2,:) = stims(stims>0); stims = stims2;
%     for i=1:size(stims,1)
%         indx1 = stim==stims(i,1)|stim==stims(i,2);
%         indx2 = stim==stims(i,end)|stim==stims(i,end-1);
%         s1 = mean(stim(indx1)); s2 = mean(stim(indx2));
%         f1 = mean(nspk(:,indx1)'); f2 = mean(nspk(:,indx2)');
%         dfds = (f2-f1)/(s2-s1);
%         for k=1:nunits
%             clear J;
%             these_units = nchoosek(1:nunits,k);
%             if size(these_units,1)>2000
%                 these_units = these_units(randperm(size(these_units,1)),:);
%                 these_units = these_units(1:2000,:);
%             end
%             for j=1:size(these_units,1)
%                 J(j) = 0;
%                 these_slopes = dfds(these_units(j,:));
%                 these_nspk = nspk(these_units(j,:),:);
%                 [c,~,l]=pca(these_nspk');
%                 J(j) = ((these_slopes*c(:,1))^2)/l(1); % fisher info in the first component
%             end
%             fisherinfo.mean(k,i) = mean(J);
%             fisherinfo.raw(k,i).J = J;
%         end
%     end
% else
%     indx1 = stim==stims(1)|stim==stims(2);
%     indx2 = stim==stims(end)|stim==stims(end-1);
%     s1 = mean(stim(indx1)); s2 = mean(stim(indx2));
%     f1 = mean(nspk(:,indx1)'); f2 = mean(nspk(:,indx2)');
%     dfds = (f2-f1)/(s2-s1);
%     for k=1:nunits
%         clear J;
%         these_units = nchoosek(1:nunits,k);
%         if size(these_units,1)>2000
%             these_units = these_units(randperm(size(these_units,1)),:);
%             these_units = these_units(1:2000,:);
%         end
%         for j=1:size(these_units,1)
%             these_slopes = dfds(these_units(j,:));
%             these_nspk = nspk(these_units(j,:),:);
%             [c,~,l]=pca(these_nspk');
%             J(j,1) = ((these_slopes*c(:,1))^2)/l(1); % fisher info in the first component
%             if k>=2
%                 J(j,2) = ((these_slopes*c(:,2))^2)/l(2); % fisher info in the second component
%             else
%                 J(j,2) = 0;
%             end
%         end
%         fisherinfo.mean(k,:) = mean(J,1);
%         fisherinfo.raw(k).J = J;
%     end
% end
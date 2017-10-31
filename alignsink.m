function alignsink

load csd;
%
sessions(6).ves_tsink = 229;
sessions(6).ves_zsink = 551;
sessions(34).ves_tsink = 236;
sessions(34).ves_zsink = 578;

nsessions = length(sessions);
count = 0;
for i=[6 7 8 13:21 34 35]
    fprintf(['session ' num2str(i) '\n']);
    % load session
    count = count + 1;
    csd = sessions(i).ves_csd;
    zsink = sessions(i).ves_zsink;
    tsink = sessions(i).ves_tsink;
    ves_mua = sessions(i).ves_mua;
    vis_mua = sessions(i).vis_mua;
    
    % surround by nans
    [nt,nz] = size(csd);
    csd = [nan(1000,nz) ; csd ; nan(1500-nt,nz)];
    [nt,nz] = size(csd);
    csd = [nan(nt,nz) csd nan(nt,nz)];
    [nt,nz] = size(ves_mua);
    ves_mua = [nan(1000,nz) ; ves_mua ; nan(1500-nt,nz)];
    [nt,nz] = size(ves_mua);
    ves_mua = [nan(nt,nz) ves_mua nan(nt,nz)];
    [nt,nz] = size(vis_mua);
    vis_mua = [nan(1000,nz) ; vis_mua ; nan(1500-nt,nz)];
    [nt,nz] = size(vis_mua);
    vis_mua = [nan(nt,nz) vis_mua nan(nt,nz)];
    % shift to align sink
    csd2(count,:,:) = circshift(csd,[-tsink -zsink]);
    ves_mua2(count,:,:) = circshift(ves_mua,[0 -zsink]);
    vis_mua2(count,:,:) = circshift(vis_mua,[0 -zsink]);
end
x=1;

figure; imagesc(squeeze(nanmean(csd2))');
figure; imagesc(squeeze(nanmean(ves_mua2))',[0 1]);
figure; imagesc(squeeze(nanmean(vis_mua2([1:6 9 13 14],:,:)))',[0 1]);
set(gca,'XTick',[1100 1180 1320 1420],'XTickLabel',[0 0.5 1.2 1.6]);
set(gca,'XTick',[1100 1180 1320 1420]-250,'XTickLabel',[0 0.5 1.2 1.6]);
set(gca,'YTick',[950-308 950-154 950 950+154 950+308 950+462],'YTickLabel',[-0.4 -0.2 0 0.2 0.4 0.6]);
hline(950,'--k');
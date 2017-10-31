function lfps = addlfps(lfp,tstim,ntrls,prs)

for i=1:ntrls
    lfps(i).v = lfp.v(lfp.t>(tstim(i).on + prs.tspk(1)) & lfp.t<(tstim(i).off + prs.tspk(2)));
    lfps(i).t = lfp.t(lfp.t>(tstim(i).on + prs.tspk(1)) & lfp.t<(tstim(i).off + prs.tspk(2))) - tstim(i).on;
    lfps(i).v = downsample(lfps(i).v,prs.Nd);
    lfps(i).t = downsample(lfps(i).t,prs.Nd);
    nt(i) = length(lfps(i).t);
end

nt = min(nt);
for i=1:ntrls
    lfps(i).v = lfps(i).v(1:nt);
    lfps(i).t = lfps(i).t(1:nt)';
end
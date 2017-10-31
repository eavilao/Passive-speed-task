function spks = addspks(unit,tstim,ntrls,prs)

for i=1:ntrls
    spks(i).tspk = unit.tspk(unit.tspk>(tstim(i).on + prs.tspk(1)) & unit.tspk<(tstim(i).off + prs.tspk(2)))-tstim(i).on;
    spks(i).nspk = numel(unit.tspk(unit.tspk>(tstim(i).on + prs.nspk(1)) & unit.tspk<(tstim(i).off + prs.nspk(2))));
end
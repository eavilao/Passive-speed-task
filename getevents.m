function [ntrls,tfix,tstim] = getevents(fname)

[n, ts, sv] = plx_event_ts(fname, 257); count=0;
for i=1:length(sv)
    if sv(i)==12
        count=count+1;
        tfix(count).on=ts(i-3);
        tstim(count).on=ts(i-2);
        tstim(count).off=ts(i-1);
    end
end
ntrls = count;
function [ntrls,tfix,tstim,sacc] = getevents_HD(fname)

[n, ts, sv] = plx_event_ts(fname, 257); count=0;
for i=1:length(sv)
    if sv(i)==12
        count=count+1;
        tfix(count).on=ts(i-6);
        tstim(count).on=ts(i-5);
        tstim(count).off=ts(i-4);
        tstim(count).targetOn = ts(i-3);
        tstim(count).saccadeOnset = ts(i-2);
        tstim(count).saccOnTarg = ts(i-1);
        if sv(i-1)==8
            sacc(count).direction = 1;
        elseif sv(i-1)==9
            sacc(count).direction = -1;
        end
    end
end
ntrls = count;
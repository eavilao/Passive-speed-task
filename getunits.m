function units = getunits(fname,getall,getch)

nunits=0;
tscounts = plx_info(fname,0);                                               % changed from fullread=1 to 0 on 051815
[nrow,ncol] = size(tscounts);
if getall==1, k=1; else k=2; end
if nargin<3, getch=[]; end

if isempty(getch)
    for i=1:ncol
        for j=k:nrow
            if tscounts(j,i)>0
                nunits = nunits+1;
                units{nunits}.nspk = tscounts(j,i);
                [~,~,units(nunits).tspk,units(nunits).spkwf] = plx_waves(fname,i-1,j-1);
                units(nunits).chnl = i-1;
                units(nunits).unit = j-1;
            end
        end
    end
else
    for j=k:nrow
        if tscounts(j,getch+1)>0
            nunits = nunits+1;
            units(nunits).nspk = tscounts(j,getch+1);
            [~,~,units(nunits).tspk,units(nunits).spkwf] = plx_waves(fname,getch,j-1);
            units(nunits).chnl = getch;
            units(nunits).unit = j-1;
        end
    end
end
function stim = getstim_HD(fname,prs,ntrls)

fid = fopen(fname, 'r');
eof=0; newline = 'nonewline'; count=0;
count=0;

while ~eof
    while ~strcmp(newline(1:7),prs.tempo(1).label)
        newline = fgetl(fid);
    end
    if strcmp(newline(9),'0') || strcmp(newline(9),'5') % 0=correct 5=incorrect
        count = count+1;
        get_speed = 0; get_modality = 0;
        stim(count).choice = str2num(newline(prs.tempo(1).col)); 
        
        while  (~get_speed && ~get_modality)% Use this if you want to run angularspeed %~get_speed %Use this to run linearspeed 
            newline = fgetl(fid);
            get_modality = strcmp(newline(1:9),prs.tempo(2).label); % STIM_TYPE Vis/Ves/Combined
            get_speed = strcmp(newline(1:7),prs.tempo(3).label); % HEADING
        end
        
        
        if get_speed
            stim(count).speed = str2num(newline(prs.tempo(3).col));
            while ~get_modality
                newline = fgetl(fid);
                get_modality = strcmp(newline(1:9),prs.tempo(2).label);
            end
            stim(count).modality = str2num(newline(prs.tempo(2).col));
        elseif get_modality
            stim(count).modality = str2num(newline(prs.tempo(2).col));
        end
         
        while ~get_speed
            newline = fgetl(fid);
            get_speed = strcmp(newline(1:7),prs.tempo(3).label);
            stim(count).speed = str2num(newline(prs.tempo(3).col));
        end
    else
        newline = fgetl(fid);
    end
    if count==ntrls, eof=1; end
end


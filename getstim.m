function stim = getstim(fname,prs,ntrls)

fid = fopen(fname, 'r');
eof=0; newline = 'nonewline'; count=0;
count=0;

while ~eof
    while ~strcmp(newline(1:7),prs.tempo(1).label)
        newline = fgetl(fid);
    end
    if strcmp(newline(prs.tempo(1).col),prs.tempo(1).val)
        count = count+1;
        get_speed = 0; get_modality = 0;
        
        while  (~get_speed && ~get_modality)% Use this if you want to run angularspeed %~get_speed %Use this to run linearspeed 
            newline = fgetl(fid);
            %get_speed = strcmp(newline(1:9),prs.tempo(2).label); % For speed protocol
            get_speed = strcmp(newline(1:7),prs.tempo(2).label); %For 1DAzi
            get_modality = strcmp(newline(1:9),prs.tempo(3).label);
        end
        
        
        if get_speed
            stim(count).speed = str2num(newline(prs.tempo(2).col));
            while ~get_modality
                newline = fgetl(fid);
                get_modality = strcmp(newline(1:9),prs.tempo(3).label);
            end
            stim(count).modality = str2num(newline(prs.tempo(3).col));
        elseif get_modality
            stim(count).modality = str2num(newline(prs.tempo(3).col));
        end
        
        
        while ~get_speed
            newline = fgetl(fid);
            get_speed = strcmp(newline(1:7),prs.tempo(2).label); %For 1DAzi
            %get_speed = strcmp(newline(1:9),prs.tempo(2).label); % For speed protocol
        end
        %stim(count).speed = str2num(newline(prs.tempo(2).col)); % Case in '1DAzi' speed == azimuth        EA 04/05/2017
    else
        newline = fgetl(fid);
    end
    if count==ntrls, eof=1; end
end
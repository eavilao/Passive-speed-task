[FileName,pathname] = uigetfile('*.log','Select Log Files','Z:\Data\MOOG','MultiSelect','on');
if iscell(FileName)
    FileName = sort(FileName);
    nBlocks = length(FileName);
else
    nBlocks = 1;
end

AllCorrect = 45;
AllTotal = 45;
AllFB = 0;
if ~isempty(strfind(pathname,'Moran'))
    n = round(10*rand(1));
    if n==0||n==6,
        disp('Moran ate my homework.')
    elseif n==1||n==7,
        disp('Moran ate his weight in paper today.')
    elseif n==2,
        disp('Moran the moron')
    elseif n==3||n==8,
        disp('Moran is my best friend.')
    elseif n==4||n==9,
        disp('Moran completed a billion trials today.')
    elseif n==5||10 
        disp('Moran is a big disappointment.')
    end
end
if ~isempty(strfind(pathname,'Jasper'))
%     n = round(10*rand(1));
%     if n==0||n==6,
%         disp('Jasper punched me in the nose today.')
%     elseif n==1||n==7,
%         disp('Jasper hurt my feelings.')
%     elseif n==2,
%         disp('Jasper is a great listener.')
%     elseif n==3||n==8,
%         disp('Jasper can count to 100!')
%     elseif n==4||n==9,
%         disp('Jasper is my best friend.')
%     elseif n==5||10 
%         disp('Jasper is the most beautiful monkey.')
%     end
    disp('You are ruining my monkey! - The Ghost of Natascha')
end
if ~isempty(strfind(pathname,'Quigley'))
    n = round(10*rand(1));
    if n==0,
        disp('Quigley is dizzy from spinning in his chair!')
    elseif n==1,
        disp('Quigley is upside down!')
    elseif n==2,
        disp('"Oooooooooh" - Quigley')
      elseif n==3,
        disp('Quigley ate a ball of tape.')
    elseif n==4||n==9,
        disp('I love Quigley')
    elseif n==5||10 
        disp('Quigley is coming to get you!')
    elseif n==6
%         disp('Quigley drew the mural in the monkey room')
        disp('I found Quigleys finger on the floor today')
    elseif n==7
        disp('"Where did my finger go?" - Quigley')
    elseif n==8
        disp('Quigley is Wiggly')
    end
end


for block = 1:nBlocks,
    if nBlocks>1,
        fullName = strcat(pathname,FileName{block});
    else
        fullName = strcat(pathname,FileName);
    end
    fid = fopen(fullName, 'r');

    %skips irregularly formatted rows in beginning of log file
    fseek(fid,54,'bof');

    %import log file into cell array
    Log = textscan(fid,'%s %s %f %f');

    %Find trial, heading, stimulus type, and outcome fields
    TrialIdx = strmatch('TRIAL#',Log{1,1});
    HeadIdx = strmatch('HEADING',Log{1,1});
    OutcomeIdx = strmatch('OUTCOME',Log{1,1});
    StimIdx = strmatch('STIM_TYPE',Log{1,1});
    TTidx = strmatch('TARG_LUM_MULT',Log{1,1});

    nTrials = length(TrialIdx);
    BlockCorrect = 0;
    BlockTotal = 0;
    %Fill in values for trial, heading, and outcome
    Data = zeros(nTrials,3);
    for i = 1:nTrials,
        Data(i,1) = str2double(Log{1,2}{TrialIdx(i)});
        Data(i,2) = str2double(Log{1,2}{HeadIdx(i)});
        Data(i,3) = str2double(Log{1,2}{OutcomeIdx(i)}(1));
        Data(i,4) = Log{1,3}(TTidx(i));
        Data(i,5) = Log{1,4}(TTidx(i));
        Data(i,6) = str2double(Log{1,2}{StimIdx(i)});
    end
    
    %Determine whether vestibular, visual, or combined stimulus
    stimtypes = unique(Data(:,6));
    if ismember(1,stimtypes)&&~ismember(2,stimtypes)&&~ismember(3,stimtypes)
        %counters
        nHead = length(unique(Data(:,2)));
        Headings = unique(Data(:,2));
        Correct = zeros(size(Headings));
        Total = zeros(size(Headings));
        Percent = zeros(size(Headings));
        GoodTrials = sort([find((Data(:,3)==0)); find((Data(:,3)==5))]);

        for j = GoodTrials',
            for k = 1:nHead,
                if Data(j,2)==Headings(k)
                    if Data(j,4)&&Data(j,5)
                        BlockTotal = BlockTotal+1;
                        Total(k) = Total(k)+1;
                        if ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                            Correct(k) = Correct(k)+1;
                        end
                    elseif ~Data(j,3)
                        BlockCorrect = BlockCorrect+1;
                        BlockTotal = BlockTotal+1;
                    else
                        BlockTotal = BlockTotal+1;
                    end
                end
            end
        end

    %     FixBreaks = nTrials-length(GoodTrials);
    %     BlockCorrect = sum(Correct);
    %     BlockTotal = sum(Total);
        Percent = [Correct;BlockCorrect]./[Total;BlockTotal];

        AllCorrect = AllCorrect + BlockCorrect;
        AllTotal = AllTotal + BlockTotal;
    %     AllFB = AllFB + FixBreaks;

        strings = cell(nHead,1);
        for ii = 1:nHead
            if ii == 1
                strings{ii} = sprintf('Block %d:\n%c 2TT; EAP \nVestibular Only:\n%g: %d/%d - %d%c',...
                block,'%',Headings(1),Correct(1),Total(1),round(100*Percent(1)),'%');
            elseif ii == nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nTotal: %d/%d - %d%c\nFB: \nThreshold:\nBias:\n',...
                    Headings(nHead),Correct(nHead),Total(nHead),round(100*Percent(nHead)),'%',...
                BlockCorrect,BlockTotal,round(100*Percent(nHead+1)),'%');
            else
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii),Correct(ii),Total(ii),round(100*Percent(ii)),'%'); 
            end
        end

        if block == nBlocks
            strings{nHead+1} = sprintf('Total Trials: %d/%d - %d%c',AllCorrect,AllTotal,round(100*AllCorrect/AllTotal),'%');
            strings{nHead+2} = sprintf('Juice: \nmL/trial: ');
%             strings{nHead+3} = sprintf('FB: %d',AllFB);
        end

        for jj = 1:length(strings)
            disp(strings{jj})
        end
    elseif ~ismember(1,stimtypes)&&ismember(2,stimtypes)&&~ismember(3,stimtypes)
        %counters
        nHead = length(unique(Data(:,2)));
        Headings = unique(Data(:,2));
        Correct = zeros(size(Headings));
        Total = zeros(size(Headings));
        Percent = zeros(size(Headings));
        GoodTrials = sort([find((Data(:,3)==0)); find((Data(:,3)==5))]);

        for j = GoodTrials',
            for k = 1:nHead,
                if Data(j,2)==Headings(k)
                    if Data(j,4)&&Data(j,5)
                        BlockTotal = BlockTotal+1;
                        Total(k) = Total(k)+1;
                        if ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                            Correct(k) = Correct(k)+1;
                        end
                    elseif ~Data(j,3)
                        BlockCorrect = BlockCorrect+1;
                        BlockTotal = BlockTotal+1;
                    else
                        BlockTotal = BlockTotal+1;
                    end
                end
            end
        end

    %     FixBreaks = nTrials-length(GoodTrials);
    %     BlockCorrect = sum(Correct);
    %     BlockTotal = sum(Total);
        Percent = [Correct;BlockCorrect]./[Total;BlockTotal];

        AllCorrect = AllCorrect + BlockCorrect;
        AllTotal = AllTotal + BlockTotal;
    %     AllFB = AllFB + FixBreaks;

        strings = cell(nHead,1);
        for ii = 1:nHead
            if ii == 1
                strings{ii} = sprintf('Block %d:\n%c 2TT; EAP \nVisual Only\n%g: %d/%d - %d%c',...
                block,'%',Headings(1),Correct(1),Total(1),round(100*Percent(1)),'%');
            elseif ii == nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nTotal: %d/%d - %d%c\nFB: \nThreshold:\nBias:\n',...
                    Headings(nHead),Correct(nHead),Total(nHead),round(100*Percent(nHead)),'%',...
                BlockCorrect,BlockTotal,round(100*Percent(nHead+1)),'%');
            else
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii),Correct(ii),Total(ii),round(100*Percent(ii)),'%'); 
            end
        end

        if block == nBlocks
            strings{nHead+1} = sprintf('Total Trials: %d/%d - %d%c',AllCorrect,AllTotal,round(100*AllCorrect/AllTotal),'%');
            strings{nHead+2} = sprintf('Juice: \nmL/trial: ');
%             strings{nHead+3} = sprintf('FB: %d',AllFB);
        end

        for jj = 1:length(strings)
            disp(strings{jj})
        end
    elseif ismember(1,stimtypes)&&ismember(3,stimtypes)&&~ismember(2,stimtypes)
        
        %counters
        Headings = unique(Data(:,2));
        nHead = length(Headings);
        VestCorrect = zeros(size(Headings));
        VestTotal = zeros(size(Headings));
        VestPercent = zeros(size(Headings));
        CombCorrect = zeros(size(Headings));
        CombTotal = zeros(size(Headings));
        CombPercent = zeros(size(Headings));
        GoodTrials = sort([find((Data(:,3)==0)); find((Data(:,3)==5))]);

        for j = GoodTrials',
            for k = 1:nHead,
                if Data(j,2)==Headings(k)
                     BlockTotal = BlockTotal+1;
                     if Data(j,6)==1
                         if Data(j,4)&&Data(j,5)
                            VestTotal(k) = VestTotal(k)+1;
                            if ~Data(j,3)
                                BlockCorrect = BlockCorrect+1;
                                VestCorrect(k) = VestCorrect(k)+1;
                            end
                         elseif ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                         end
                     elseif Data(j,6)==3
                         if Data(j,4)&&Data(j,5)
                            CombTotal(k) = CombTotal(k)+1;
                            if ~Data(j,3)
                                BlockCorrect = BlockCorrect+1;
                                CombCorrect(k) = CombCorrect(k)+1;
                            end
                         elseif ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                         end
                     end
                end
            end
        end

    %     FixBreaks = nTrials-length(GoodTrials);
    %     BlockCorrect = sum(Correct);
    %     BlockTotal = sum(Total);
        VestPercent = VestCorrect./VestTotal;
        CombPercent = CombCorrect./CombTotal;
        BlockPercent = BlockCorrect./BlockTotal;
        AllCorrect = AllCorrect + BlockCorrect;
        AllTotal = AllTotal + BlockTotal;
    %     AllFB = AllFB + FixBreaks;

    
        strings = cell(2*nHead,1);
        for ii = 1:2*nHead
            if ii == 1
                strings{ii} = sprintf('Block %d:\n%c 2TT; EAP \n\nVestibular Only: \n%g: %d/%d - %d%c',...
                block,'%',Headings(1),VestCorrect(1),VestTotal(1),round(100*VestPercent(1)),'%');
            elseif ii == nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nThreshold:\nBias:\n\nCombined Visual and Vestibular:',...
                    Headings(nHead),VestCorrect(nHead),VestTotal(nHead),round(100*VestPercent(nHead)),'%');
            elseif ii<nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii),VestCorrect(ii),VestTotal(ii),round(100*VestPercent(ii)),'%'); 
            elseif ii==2*nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nThreshold:\nBias:\nTotal: %d/%d - %d%c\nFB:\n',...
                    Headings(ii-nHead),CombCorrect(ii-nHead),CombTotal(ii-nHead),round(100*CombPercent(ii-nHead)),'%',BlockCorrect,BlockTotal,round(100*BlockPercent),'%');
            elseif ii>nHead&&ii<2*nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii-nHead),CombCorrect(ii-nHead),CombTotal(ii-nHead),round(100*CombPercent(ii-nHead)),'%');
            end
        end

        if block == nBlocks
            strings{2*nHead+1} = sprintf('Total Trials: %d/%d - %d%c',AllCorrect,AllTotal,round(100*AllCorrect/AllTotal),'%');
            strings{2*nHead+2} = sprintf('Juice: \nmL/trial: ');
            strings{2*nHead+3} = sprintf('FB: ');
        end

        for jj = 1:length(strings)
            disp(strings{jj})
        end
    elseif ismember(1,stimtypes)&&ismember(3,stimtypes)&&ismember(2,stimtypes)
        
        %counters
        Headings = unique(Data(:,2));
        nHead = length(Headings);
        VestCorrect = zeros(size(Headings));
        VestTotal = zeros(size(Headings));
        VestPercent = zeros(size(Headings));
        CombCorrect = zeros(size(Headings));
        CombTotal = zeros(size(Headings));
        VisCorrect = zeros(size(Headings));
        VisTotal = zeros(size(Headings));
        CombPercent = zeros(size(Headings));
        GoodTrials = sort([find((Data(:,3)==0)); find((Data(:,3)==5))]);

        for j = GoodTrials',
            for k = 1:nHead,
                if Data(j,2)==Headings(k)
                     BlockTotal = BlockTotal+1;
                     if Data(j,6)==1
                         if Data(j,4)&&Data(j,5)
                            VestTotal(k) = VestTotal(k)+1;
                            if ~Data(j,3)
                                BlockCorrect = BlockCorrect+1;
                                VestCorrect(k) = VestCorrect(k)+1;
                            end
                         elseif ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                         end
                     elseif Data(j,6)==3
                         if Data(j,4)&&Data(j,5)
                            CombTotal(k) = CombTotal(k)+1;
                            if ~Data(j,3)
                                BlockCorrect = BlockCorrect+1;
                                CombCorrect(k) = CombCorrect(k)+1;
                            end
                         elseif ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                         end
                     elseif Data(j,6)==2,
                         if Data(j,4)&&Data(j,5)
                            VisTotal(k) = VisTotal(k)+1;
                            if ~Data(j,3)
                                BlockCorrect = BlockCorrect+1;
                                VisCorrect(k) = VisCorrect(k)+1;
                            end
                         elseif ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                         end
                     end
                end
            end
        end

    %     FixBreaks = nTrials-length(GoodTrials);
    %     BlockCorrect = sum(Correct);
    %     BlockTotal = sum(Total);
        VestPercent = VestCorrect./VestTotal;
        CombPercent = CombCorrect./CombTotal;
        VisPercent = VisCorrect./VisTotal;
        BlockPercent = BlockCorrect./BlockTotal;
        AllCorrect = AllCorrect + BlockCorrect;
        AllTotal = AllTotal + BlockTotal;
    %     AllFB = AllFB + FixBreaks;

    
        strings = cell(3*nHead,1);
        for ii = 1:3*nHead
            if ii == 1
                strings{ii} = sprintf('Block %d:\n%c 2TT; EAP \n\nVestibular Only: \n%g: %d/%d - %d%c',...
                block,'%',Headings(1),VestCorrect(1),VestTotal(1),round(100*VestPercent(1)),'%');
            elseif ii == nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nThreshold:\nBias:\n\nVisual Only:',...
                    Headings(nHead),VestCorrect(nHead),VestTotal(nHead),round(100*VestPercent(nHead)),'%');
            elseif ii<nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii),VestCorrect(ii),VestTotal(ii),round(100*VestPercent(ii)),'%'); 
            elseif ii>nHead&&ii<2*nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii-nHead),VisCorrect(ii-nHead),VisTotal(ii-nHead),round(100*VisPercent(ii-nHead)),'%');
            elseif ii==3*nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nThreshold:\nBias:\nTotal: %d/%d - %d%c\nFB:\n',...
                    Headings(ii-2*nHead),CombCorrect(ii-2*nHead),CombTotal(ii-2*nHead),round(100*CombPercent(ii-2*nHead)),'%',BlockCorrect,BlockTotal,round(100*BlockPercent),'%');
            elseif ii>2*nHead&&ii<3*nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii-2*nHead),CombCorrect(ii-2*nHead),CombTotal(ii-2*nHead),round(100*CombPercent(ii-2*nHead)),'%');
            elseif ii==2*nHead
                strings{ii}= sprintf('%g: %d/%d - %d%c\nThreshold:\nBias:\n\nCombined Visual and Vestibular:',...
                    Headings(ii-nHead),VisCorrect(ii-nHead),VisTotal(ii-nHead),round(100*VisPercent(ii-nHead)),'%');
            end
        end

        if block == nBlocks
            strings{3*nHead+1} = sprintf('Total Trials: %d/%d - %d%c',AllCorrect,AllTotal,round(100*AllCorrect/AllTotal),'%');
            strings{3*nHead+2} = sprintf('Juice: \nmL/trial: ');
            strings{3*nHead+3} = sprintf('FB: ');
        end

        for jj = 1:length(strings)
            disp(strings{jj})
        end
    elseif ismember(3,stimtypes)&&~ismember(1,stimtypes)&&~ismember(2,stimtypes)
        
        %counters
        Headings = unique(Data(:,2));
        nHead = length(Headings);
        CombCorrect = zeros(size(Headings));
        CombTotal = zeros(size(Headings));
        CombPercent = zeros(size(Headings));
        GoodTrials = sort([find((Data(:,3)==0)); find((Data(:,3)==5))]);

        for j = GoodTrials',
            for k = 1:nHead,
                if Data(j,2)==Headings(k)
                     BlockTotal = BlockTotal+1;
                     if Data(j,4)&&Data(j,5)
                        CombTotal(k) = CombTotal(k)+1;
                        if ~Data(j,3)
                            BlockCorrect = BlockCorrect+1;
                            CombCorrect(k) = CombCorrect(k)+1;
                        end
                     elseif ~Data(j,3)
                        BlockCorrect = BlockCorrect+1;
                     end
                end
            end
        end

    %     FixBreaks = nTrials-length(GoodTrials);
    %     BlockCorrect = sum(Correct);
    %     BlockTotal = sum(Total);
        CombPercent = CombCorrect./CombTotal;
        BlockPercent = BlockCorrect./BlockTotal;
        AllCorrect = AllCorrect + BlockCorrect;
        AllTotal = AllTotal + BlockTotal;
    %     AllFB = AllFB + FixBreaks;

    
        strings = cell(nHead,1);
        for ii = 1:nHead
            if ii == 1
                strings{ii} = sprintf('Block %d:\n%c 2TT; EAP \n\nCombined Visual and Vestibular: \n%g: %d/%d - %d%c',...
                block,'%',Headings(1),CombCorrect(1),CombTotal(1),round(100*CombPercent(1)),'%');
            elseif ii<nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c',...
                    Headings(ii),CombCorrect(ii),CombTotal(ii),round(100*CombPercent(ii)),'%'); 
            elseif ii==nHead
                strings{ii} = sprintf('%g: %d/%d - %d%c\nThreshold:\nBias:\nTotal: %d/%d - %d%c\nFB:\n',...
                    Headings(nHead),CombCorrect(nHead),CombTotal(nHead),round(100*CombPercent(nHead)),'%',BlockCorrect,BlockTotal,round(100*BlockPercent),'%');
            end
        end

        if block == nBlocks
            strings{nHead+1} = sprintf('Total Trials: %d/%d - %d%c',AllCorrect,AllTotal,round(100*AllCorrect/AllTotal),'%');
            strings{nHead+2} = sprintf('Juice: \nmL/trial: ');
            strings{nHead+3} = sprintf('FB: ');
        end

        for jj = 1:length(strings)
            disp(strings{jj})
        end
    end
  
end





fclose(fid);
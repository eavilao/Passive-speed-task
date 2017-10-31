function analyse_linearspeed5

analyseAll = input('Analyse all files? Yes (1) or No (0): ');
unit_type = input('singleunit or multiunit? 1-multiunit, 2-singleunit: ');
if analyseAll
    monkeyInfoFile_linearSpeed
else
    monkeyInfo(1).monkey=input('monkey #: ');
    monkeyInfo(1).session = input('session #: ','s');
    monkeyInfo(1).run = input('run #: ');
    monkeyInfo(1).row = input('row #: ');
    monkeyInfo(1).column = input('column #: ');
    monkeyInfo(1).depth = input('Depth below gray matter (in microns) #: ');
    monkeyInfo(1).channel = input('channel #? (all channels=0) : ');
    monkeyInfo(1).units = input('All units? (sorted units=0; all units=1): = ');
end

recList = squeeze(struct2cell(monkeyInfo));
recordingsToAnalyse = find(cellfun(@(x) ~isempty(x), recList(1, :)));

for numOfFile = recordingsToAnalyse;
    monk = monkeyInfo(numOfFile).monkey;
    sess = monkeyInfo(numOfFile).session;
    run = monkeyInfo(numOfFile).run;
    loc(1) = monkeyInfo(numOfFile).row;
    loc(2) = monkeyInfo(numOfFile).column;
    loc(3) = monkeyInfo(numOfFile).depth;
    getch = monkeyInfo(numOfFile).channel;
    all = monkeyInfo(numOfFile).units;
    
    %Read out the plexon and tempo files
    file1 = ['m' num2str(monk) 'c' num2str(sess) 'r' num2str(run) '-01.plx'];
    file2 = ['m' num2str(monk) 'c' num2str(sess) 'r' num2str(run) '.log'];
    
    if monk == 45
        path1 = 'C:\Users\erico\Documents\Jasper\Jasper_spike_sorted_files\'; % For Jasper
        path2 = 'Z:\Data\MOOG\Jasper\Raw\'; % For Jasper
    elseif monk == 44
        path1 = 'C:\Users\erico\Documents\Quigley\Quigley_spike_sorted\'; % For Quigley
        path2 = 'Z:\Data\MOOG\Quigley\Eric_Project\'; % For Quigley
        
    else
        path1 = 'C:\Users\erico\Documents\Oscar\Oscar_LinSpeed\'; % For Oscar
        path2 = 'C:\Users\erico\Documents\Oscar\Tempo\'; % For Oscar
        
    end
    
    
    disp([' Loading file ' num2str(file1) '  --   Processing file  ' num2str(numOfFile) ' of ' num2str(max(recordingsToAnalyse)) '  --  ' ])
    
    clear ch

    nchan = 16;
    [ntrls,tfix,tstim] = getevents([path1 file1]);
    if getch==0, getch=1:nchan; end
    for k=getch
        fprintf(['channel # ' num2str(k) '\n']);
        units = getunits([path1 file1],all,k); nunits(k)=length(units);
        for j=1:nunits(k)
            ch(k).unit(j).wfspk = units{j}.wfspk;
            for i=1:ntrls
                ch(k).unit(j).trl(i).tspk = units{j}.tspk(units{j}.tspk>(tfix(i).on-0.2) & units{j}.tspk<(tstim(i).off+0.3))-tstim(i).on;
                ch(k).unit(j).trl(i).nspk = numel(units{j}.tspk(units{j}.tspk>(tstim(i).on+0.5) & units{j}.tspk<(tstim(i).off-0.5)));
            end
            linSpeedUnits(k,:) = [k j];
            
            % Plot waveforms
            % Mean + std dev
            fname2 = ['m' num2str(monk) 's' num2str(sess) 'r' num2str(run) 'ch' num2str(k) 'u' num2str(j)];
            shadedErrorBar([],mean(ch(k).unit(j).wfspk),std(ch(k).unit(j).wfspk), 'k');
            box off
            title(fname2)
            set(gcf, 'Color', 'w');
            fname3 = ['Wf_m' num2str(monk) 's' num2str(sess) 'r' num2str(run)]; 
            print(fname3,'-append', '-dpsc2')
            
        end
    end
    
    linSpeedUnits
    clear stim stimtype stimref stim_cmp
    
    
    fid = fopen([path2 file2], 'r'); eof=0; newline = 'nonewline'; count=0;
    while ~eof
        while ~strcmp(newline(1:7),'OUTCOME')
            newline = fgetl(fid);
        end
        if strcmp(newline(9),'0')
            count = count+1;
            while ~strcmp(newline(1:9),'AMPLITUDE')
                newline = fgetl(fid);
            end
            stim(count)=str2num(newline(17:21));
            while ~strcmp(newline(1:9),'STIM_TYPE')
                newline = fgetl(fid);
            end
            stimtype(count)=str2num(newline(11:15));
        else
            newline = fgetl(fid);
        end
        if count==ntrls, eof=1; end
    end
    
    stim2=repmat(stim,numel(unique(stim)),1);
    stimref = repmat(unique(stim)',1,ntrls);
    stim_cmp = (stim2==stimref);   %ERROR
    
    stimtype2 = repmat(stimtype,numel(unique(stimtype)),1);
    stimtypes = repmat(unique(stimtype)',1,ntrls);
    stimtype_cmp = (stimtype2==stimtypes);
    
    % correlograms
    % for i=1:nunits(k)
    %     for j=1:nunits(k)
    %         Rij=crosscorrelogram(units{i}.tspk,units{j}.tspk,[-0.2 0.2]);
    %         R(i,j,:)=hist(Rij,100);
    %     end
    % end
    
    % tuning curves and psth
    for k=getch
        for j=1:nunits(k)
            for l=1:3
                for i=1:5
                    trl2=ch(k).unit(j).trl(stimtype_cmp((i==1)*i + (i~=1)*(l+1),:) & stim_cmp(i,:));
                    ch(k).unit(j).fr(l,i).mu=mean([trl2.nspk]/.5);
                    ch(k).unit(j).fr(l,i).sig=std([trl2.nspk]/.5);
                    if ch(k).unit(j).fr(l,i).mu<2,
                        ch(k).unit(j).psth(l,i).r = [];
                        ch(k).unit(j).psth(l,i).t = [];
                        ch(k).unit(j).psth(l,i).e = [];
                    else
                        %[ch(k).unit(j).psth(l,i).r,ch(k).unit(j).psth(l,i).t,ch(k).unit(j).psth(l,i).e]=psth(trl2,0.025,'n');
                        [ch(k).unit(j).psth(l,i).r,ch(k).unit(j).psth(l,i).t,ch(k).unit(j).psth(l,i).e]=psth(trl2,0.05,'n');
                    end
                end
            end
        end
    end
    
    % statistics
    for k=getch
        for j=1:nunits(k)
            for l=1:3
                SSE = 0;
                for i=2:5
                    trl2=ch(k).unit(j).trl(stimtype_cmp(l+1,:) & stim_cmp(i,:));
                    heading(i-1).nspk=cell2mat({trl2.nspk});
                    SSE = SSE + sum((cell2mat({trl2.nspk})-mean(cell2mat({trl2.nspk}))).^2);
                end
                ch(k).unit(j).DDI(l) = (max([ch(k).unit(j).fr(l,2:5).mu])-min([ch(k).unit(j).fr(l,2:5).mu]))/((max([ch(k).unit(j).fr(l,2:5).mu])-min([ch(k).unit(j).fr(l,2:5).mu]))+2*sqrt(SSE/16));
                groupname=[];
                for m=1:length(heading)
                    numtrials(m) = numel(heading(m).nspk);
                    groupname=[groupname repmat({num2str(m)},[1 numtrials(m)])];
                end
                nspk= cell2mat({heading.nspk});
                ch(k).unit(j).pval(l)=anova1(nspk,groupname,'off');   % Check for speed tuning
            end
        end
    end
    
    ntimes=200; ts2=linspace(-0.3,1.8,ntimes);
    % multivariate statistics
    multi_nspk = []; multi_tspk=[];
    for k=getch
        if unit_type==1, theseunits = 1; else theseunits=2:nunits(k); end
        for j=theseunits
            for l=1:3
                SSE = 0;
                for i=1:5
                    trl2=ch(k).unit(j).trl(stimtype_cmp((i==1)*i + (i~=1)*(l+1),:)& stim_cmp(i,:));
                    heading(i).nspk=cell2mat({trl2.nspk});
                    heading(i).r=ch(k).unit(j).psth(l,i).r;
                    heading(i).e=ch(k).unit(j).psth(l,i).e;
                    ts=ch(k).unit(j).psth(l,i).t;
                    if ~isempty(ts)
                        r2(i,:)=interp1(ts,heading(i).r,ts2,'spline',0);
                        e2(i,:)=interp1(ts,heading(i).e,ts2,'spline',0);
                    else
                        r2(i,:) = zeros(1,ntimes);
                        e2(i,:) = zeros(1,ntimes);
                    end
                end
                groupname=[];
                for m=1:length(heading)
                    numtrials(m) = numel(heading(m).nspk);
                    groupname=[groupname repmat({num2str(m)},[1 numtrials(m)])];
                end
                multi_nspk = [multi_nspk cell2mat({heading.nspk})'];
                multi_tspk = cat(3,multi_tspk,r2);
                clear r2 e2;
            end
        end
    end
    
    multi = []; multi_nspk2 = []; multi_tspk2 = []; Xr = []; Yr = [];
    nstim = numel(unique(groupname)); nreps = size(multi_nspk,1)/nstim;
    if size(multi_nspk,2)>3
        for l=1:3
            multi_nspk2(l,:,:) = multi_nspk(:,mod(1:size(multi_nspk,2),3)==mod(l,3));
            multi_tspk2(l,:,:,:) = multi_tspk(:,:,mod(1:size(multi_nspk,2),3)==mod(l,3));
            X = squeeze(multi_nspk2(l,:,:));
            Y = squeeze(multi_tspk2(l,:,:,:)); Y(isnan(Y)) = 0;
            multi.nspk(l,:,:) = X;
            multi.tspk(l,:,:,:) = Y;
            Xm = (reshape(mean(reshape(X,[nstim,size(X,2)*nstim])),[nstim,size(X,2)]));
            [coeff,score,latent] = pca(Xm);
            multi.coeff(l,:,:) = coeff;
            Xr(l,:,:) = X*coeff;
            for m=1:ntimes
                Yr(l,:,m,:) = squeeze(Y(:,m,:))*coeff;
            end
            multi.mu(l,:,1) = mean(reshape(squeeze(Xr(l,:,1)),[nstim nreps]));    % first dimension
            multi.sig(l,:,1) = std(reshape(squeeze(Xr(l,:,1)),[nstim nreps]))/sqrt(nreps);
            multi.mu(l,:,2) = mean(reshape(squeeze(Xr(l,:,2)),[nstim nreps]));    % second dimension
            multi.sig(l,:,2) = std(reshape(squeeze(Xr(l,:,2)),[nstim nreps]))/sqrt(nreps);
            multi.r(l,:,:,:) = Yr(l,:,:,1:2); multi.t = ts2;
            [multi.dim(l),multi.pval(l,:)]=manova1(squeeze(Xr(l,:,:)),groupname);
        end
    elseif size(multi_nspk,2)==3
        for l=1:3
            multi_nspk2(l,:) = squeeze(multi_nspk(:,mod(1:size(multi_nspk,2),3)==mod(l,3)));
            multi_tspk2(l,:,:) = squeeze(multi_tspk(:,:,mod(1:size(multi_nspk,2),3)==mod(l,3)));
            X = squeeze(multi_nspk2(l,:));
            Y = squeeze(multi_tspk2(l,:,:)); Y(isnan(Y)) = 0;
            multi.nspk(l,:,1) = X;
            multi.tspk(l,:,:,1) = Y;
        end
    end
    
    %
    % save
    if length(getch)==1
        fOut = ['m' num2str(monk) 's' num2str(sess) 'r' num2str(run) 'ch' num2str(getch) '.mat'];
    else
        fOut = ['m' num2str(monk) 's' num2str(sess) 'r' num2str(run) '.mat'];
    end
    save(fOut,'stim','stimtype','ch','loc','multi', 'linSpeedUnits');
    
end
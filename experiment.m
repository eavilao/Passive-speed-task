classdef experiment < handle
    %%
    properties
        name                                                                % protocol
        multiunits = multiunit.empty();                                     % multiunit
        singleunits = singleunit.empty();                                   % singleunit
        lfps = lfp.empty();                                                 % lfp
        sessions = session.empty();
        populations = population.empty();
    end
    
    %%
    methods
        %% class constructor
        function this = experiment(exp_name)
            this.name = exp_name;
        end
        %% function to add sessions
        function addunits(this,monk_id,session_id,prs)
            % get params
            if nargin<2
                monk_id=input('monkey id (Oscar - 41, Quigley - 44, Jasper - 45, All - 0): ');
                session_id=input('session id (type 0 to add all sessions): ');
                prs = default_prs(this.name,monk_id);
            elseif nargin<3
                session_id=input('session id (type 0 to add all sessions): ');
                prs = default_prs(this.name,monk_id);
            elseif nargin<4
                prs = default_prs(this.name,monk_id);
            end
            % check params
            checkparams = 0;
            while ~checkparams
                checkparams = 1;
                if monk_id==0 && session_id~=0
                    checkparams = 0;
                    fprintf(['cannot add session ' num2str(session_id) ' from all monkeys \n']);
                    monk_id=input('monkey id (Oscar - 41, Quigley - 44, Jasper - 45, All - 0): ');
                    session_id=input('session id (type 0 to add all sessions): ');
                end
            end
            % warning
            if monk_id==0
                if ~isempty(this.multiunits)
                    input('replace everything? yes (1) or no (0)? ');
                    this.multiunits = multiunit.empty();
                    this.singleunits = singleunit.empty();
                    this.lfps = lfp.empty();
                    this.sessions = session.empty();
                    this.populations = population.empty();
                end
            end
            % load and analyse
            run(prs.InfoFile);
            if monk_id==0
                monk_id = unique([monkeyInfo.monk_id]);
            end
            for i=1:length(monk_id)
                monk_indx = [monkeyInfo.monk_id]==monk_id(i);
                prs = default_prs(this.name,monk_id(i));
                if session_id==0
                    session_indx = monk_indx;
                else
                    session_indx = [monkeyInfo.session_id]==session_id;
                end
                indx = find(monk_indx & session_indx);
                for j=1:length(indx)
                    this_monk_id = monkeyInfo(indx(j)).monk_id;
                    this_session_id = monkeyInfo(indx(j)).session_id;
                    this_run = monkeyInfo(indx(j)).run;
                    this_coord = monkeyInfo(indx(j)).coord;
                    this_channels = monkeyInfo(indx(j)).channels; maxchannels = max(this_channels);
                    this_units = monkeyInfo(indx(j)).units;
                    f_plx = ['m' num2str(this_monk_id) 'c' num2str(this_session_id) 'r' num2str(this_run) '-01.plx'];
                    f_log = ['m' num2str(this_monk_id) 'c' num2str(this_session_id) 'r' num2str(this_run) '.log'];
                    % read plexon file
                    fprintf(['reading ' f_plx '...   ']);
                    if strcmp(prs.InfoFile(end-1:end),'HD')
                        [ntrls,tfix,tstim] = getevents_HD([prs.filepath f_plx]);
                    else
                        [ntrls,tfix,tstim] = getevents([prs.filepath f_plx]);
                    end
                    % read log file
                    fprintf(['reading ' f_log '...\n']);
                    if strcmp(prs.InfoFile(end-1:end),'HD')
                        stim = getstim_HD([prs.filepath f_log],prs,ntrls);
                    else
                        stim = getstim([prs.filepath f_log],prs,ntrls);
                    end
                    count_multiunits = 0; count_singleunits = 0; count_lfps = 0; % counters for this session
                    try
                        for k=this_channels
                            count_multiunits = count_multiunits + 1;
                            fprintf(['... reading SPK channel ' num2str(k) '...\n']);
                            units = getunits([prs.filepath f_plx],this_units,k);
                            nunits = length(units);
                            % add multiunit
                            fprintf(['... adding multiunit ' num2str(length(this.multiunits)+1) '\n']);
                            this.multiunits(end+1) = multiunit(this_monk_id,this_session_id,k,maxchannels,this_coord,units(1).spkwf); % create a new multiunit and fills basic properties
                            
                            this.multiunits(end).add_stim(stim,tfix,tstim,prs); % fill in the attribute called stim
                            
                            this.multiunits(end).add_spks(units(1),tstim,ntrls,prs); % fill in spike times
                            this.multiunits(end).analyse_spks(this.name,prs);
                            for l=2:nunits
                                count_singleunits = count_singleunits + 1;
                                % add singleunit
                                fprintf(['... adding singleunit ' num2str(length(this.singleunits)+1),'\n']);
                                this.singleunits(end+1) = singleunit(this_monk_id,this_session_id,k,maxchannels,this_coord,units(l).spkwf);
                                this.singleunits(end).add_stim(stim,tfix,tstim,prs);
                                this.singleunits(end).add_spks(units(l),tstim,ntrls,prs);
                                this.singleunits(end).analyse_spks(this.name,prs);
                            end
                            
                            % add lfp
                            count_lfps = count_lfps + 1;
                            fprintf(['\n  ... reading AD channel ' num2str(k) '... ']);
                            this_lfp = getlfp([prs.filepath f_plx],k);
                            fprintf(['... adding lfp ' num2str(length(this.lfps)+1) '\n']);
                            this.lfps(end+1) = lfp(this_monk_id,this_session_id,k,maxchannels,this_coord);
                            this.lfps(end).add_stim(stim,tfix,tstim,prs);
                            this.lfps(end).add_lfps(this_lfp,tstim,ntrls,prs);
                            this.lfps(end).analyse_lfps(this.name,prs);
                            
                        end
                        fprintf(['***** adding session ' num2str(length(this.sessions)+1) ' *****\n']);
                        this.sessions(end+1) = session(this_monk_id,this_session_id,this_coord);
                        % multiunits
                        these_multiunits = this.multiunits(end - count_multiunits + 1:end);
                        if ~isempty(these_multiunits)
                            this.sessions(end).analyse_multiunits(this.name,these_multiunits);
                        end
                        % singleunits
                        these_singleunits = this.singleunits(end - count_singleunits + 1:end);
                        if ~isempty(these_singleunits)
                            this.sessions(end).analyse_singleunits(this.name,these_singleunits);
                            this.sessions(end).analyse_clustering(this.name,these_singleunits,these_multiunits)
                        end
                        % destroy spikes
                        n_singleunits = length(this.singleunits);
                        for ind=n_singleunits - count_singleunits + 1:n_singleunits
                            this.singleunits(ind).destroy_spks;
                            this.singleunits(ind).destroy_stim;
                        end
                        % lfps
                        these_lfps = this.lfps(end - count_lfps + 1:end);
                        if ~isempty(these_lfps)
                            this.sessions(end).analyse_lfps(this.name,these_lfps);
                        end   
                        n_lfps = length(this.lfps);
                        for ind=n_lfps - count_lfps + 1:n_lfps
                            this.lfps(ind).destroy_wave;
                            this.lfps(ind).destroy_stim;
                        end
                        catch ME
                        fprintf(['Monkey ' num2str(this_monk_id) ' ; Session ' num2str(this_session_id) ' : SESSION FAILED BECAUSE %s\n'],ME.message);
                        continue;
                    end
                end
            end
        end
        %% function to add populations
        function addpopulations(this,groupname,unitname,prs)
            if nargin<4
                prs = default_prs(this.name,0);
            end
            data_units = this.(unitname);
            data_sessions = [this.sessions.(unitname)];
            switch groupname
                case 'all'
                    this.populations(end+1) = population(unitname);
                    this.populations(end).groupname = 'all';
                    thisdata{1} = data_units;
                    thisdata{2} = data_sessions;
                    this.populations(end).analyse_population(this.name,thisdata,unitname,prs);
                case 'monks'
                    monk_ids = unique([data_units.monk_id]);
                    for i=1:length(monk_ids)
                        this.populations(end+1) = population(unitname);
                        this.populations(end).groupname = ['monkey ' num2str(monk_ids(i))];
                        thisdata{1} = data_units([data_units.monk_id]==monk_ids(i));
                        this.populations(end).analyse_population(this.name,thisdata,unitname,prs);
                    end
                case 'spkwidth'
                    for i=1:length(data_units), spkwidth(i) = data_units(i).spkwf.width; end
                    thresh_spkwidth = [prctile(spkwidth,33) prctile(spkwidth,66)];
                    % narrow spiking
                    this.populations(end+1) = population(unitname);
                    this.populations(end).groupname = 'narrow spiking';
                    thisdata{1} = data_units(spkwidth < thresh_spkwidth(1));
                    this.populations(end).analyse_population(this.name,thisdata,unitname,prs);
                    % broad spiking
                    this.populations(end+1) = population(unitname);
                    this.populations(end).groupname = 'broad spiking';
                    thisdata{1} = data_units(spkwidth > thresh_spkwidth(2));
                    this.populations(end).analyse_population(this.name,thisdata,unitname,prs);
            end
        end
        %% function to save existing singleunits
        function savesingleunit(this,unitids)
            if nargin<2
                unitids = 1:length(this.singleunits);
            end
            for i=unitids
                eval([genvarname(this.singleunits(i).tag(1:5)) '= this.singleunits(i)']);
                save(genvarname(this.singleunits(i).tag(1:5)),genvarname(this.singleunits(i).tag(1:5)));
                clear(genvarname(this.singleunits(i).tag(1:5)));
            end
        end
        %% function to plot single/multiunit/lfp data from individual electrodes
        function plotunit(this,plottype,unit_type,unit_num,unit_id)
            if nargin<5, unit_id = {}; end
            switch unit_type
                case 'multiunits'
                    if strcmp(unit_num,'all')
                        for i=1:length(this.multiunits)
                            this.multiunits(i).plot(i,this.name,plottype);
                            %                             fname2 = ['m' num2str(this.multiunits(i).monk_id) 's' num2str(this.multiunits(i).session_id) 'ch' num2str(this.multiunits(i).channel_no)];
                            %                             print(fname2,'-depsc2');
                            waitforbuttonpress; close all;
                        end
                    elseif ~isempty(unit_num)
                        this.multiunits(unit_num).plot(unit_num,this.name,plottype);
                    elseif ~isempty(unit_id)
                        unit_num = find([this.multiunits.monk_id]==unit_id{1} & ...
                            [this.multiunits.session_id]==unit_id{2} & ...
                            [this.multiunits.channel_no]==unit_id{3});
                        for i=1:length(unit_num)
                            this.multiunits(unit_num(i)).plot(unit_num(i),this.name,plottype);
                            waitforbuttonpress; close all;
                        end
                    end
                case 'singleunits'
                    if strcmp(unit_num,'all')
                        for i=1:length(this.singleunits)
                            this.singleunits(i).plot(i,this.name,plottype);
                            waitforbuttonpress; close all;
                        end
                    elseif ~isempty(unit_num)
                        this.singleunits(unit_num).plot(unit_num,this.name,plottype);
                    elseif ~isempty(unit_id)
                        unit_num = find([this.singleunits.monk_id]==unit_id{1} & ...
                            [this.singleunits.session_id]==unit_id{2} & ...
                            [this.singleunits.channel_no]==unit_id{3});
                        for i=1:length(unit_num)
                            this.singleunits(unit_num(i)).plot(unit_num(i),this.name,plottype);
                            waitforbuttonpress; close all;
                        end
                    end
                case 'lfps'
                    if strcmp(unit_num,'all')
                        for i=1:length(this.lfps)
                            this.lfps(i).plot(i,this.name,plottype);
                            waitforbuttonpress; close all;
                        end
                    elseif ~isempty(unit_num)
                        this.lfps(unit_num).plot(unit_num,this.name,plottype);
                    elseif ~isempty(unit_id)
                        unit_num = find([this.lfps.monk_id]==unit_id{1} & ...
                            [this.lfps.session_id]==unit_id{2} & ...
                            [this.lfps.channel_no]==unit_id{3});
                        for i=1:length(unit_num)
                            this.lfps(unit_num(i)).plot(unit_num(i),this.name,plottype);
                            waitforbuttonpress; close all;
                        end
                    end
            end
        end
        %% function to plot all single/multiunit/lfp data recorded simultaneously from individual sessions
        function plotsession(this,plottype,unit_type,session_num,session_id)
            if ~isempty(session_num)
                monk_id = this.sessions(session_num).monk_id;
                session_id = this.sessions(session_num).session_id;
            else
                monk_id = session_id{1};
                session_id = session_id{2};
            end
            switch unit_type
                case 'multiunits'
                    unit_num = find([this.multiunits.monk_id]==monk_id & ...
                        [this.multiunits.session_id]==session_id);
                    for i=1:length(unit_num)
                        this.multiunits(unit_num(i)).plot(unit_num(i),this.name,plottype);
                        waitforbuttonpress; close all;
                    end
                case 'singleunits'
                    unit_num = find([this.singleunits.monk_id]==monk_id & ...
                        [this.singleunits.session_id]==session_id);
                    for i=1:length(unit_num)
                        this.singleunits(unit_num(i)).plot(unit_num(i),this.name,plottype);
                        waitforbuttonpress; close all;
                    end
                case 'lfps'
                    if strcmp(plottype,'psth_average')
                        unit_num = find([this.lfps.monk_id]==monk_id & ...
                            [this.lfps.session_id]==session_id);
                        all_lfps = this.lfps(unit_num);
                        this.lfps(unit_num(1)).plot(all_lfps,this.name,plottype);
                    elseif ~any(strcmp(plottype,'csd'))
                        unit_num = find([this.lfps.monk_id]==monk_id & ...
                            [this.lfps.session_id]==session_id);
                        for i=1:length(unit_num)
                            this.lfps(unit_num(i)).plot(unit_num(i),this.name,plottype);
                            waitforbuttonpress; close all;
                        end
                    else
                        unit_num = [this.multiunits.monk_id]==monk_id & ...
                            [this.multiunits.session_id]==session_id;
                        these_multiunits = this.multiunits(unit_num);
                        this.sessions(session_num).plot(these_multiunits,plottype);
                    end
                case 'all'
                    if strcmp(plottype,'trial')
                        prs = default_prs(this.name,monk_id);
                        run(prs.InfoFile);
                        prs = default_prs(this.name,monk_id);
                        monk_indx = [monkeyInfo.monk_id]==monk_id;
                        session_indx = [monkeyInfo.session_id]==session_id;
                        indx = find(monk_indx & session_indx);
                        this_monk_id = monkeyInfo(indx).monk_id;
                        this_session_id = monkeyInfo(indx).session_id;
                        this_run = monkeyInfo(indx).run;
                        this_coord = monkeyInfo(indx).coord;
                        this_channels = monkeyInfo(indx).channels; maxchannels = max(this_channels);
                        this_units = monkeyInfo(indx).units;
                        f_plx = ['m' num2str(this_monk_id) 'c' num2str(this_session_id) 'r' num2str(this_run) '-01.plx'];
                        f_log = ['m' num2str(this_monk_id) 'c' num2str(this_session_id) 'r' num2str(this_run) '.log'];
                        % read plexon file
                        fprintf(['reading ' f_plx '...   ']);
                        [ntrls,tfix,tstim] = getevents([prs.filepath f_plx]);
                        % read log file
                        fprintf(['reading ' f_log '...\n']);
                        stim = getstim([prs.filepath f_log],prs,ntrls);
                        count_multiunits = 0; count_singleunits = 0; count_lfps = 0; % counters for this session
                        multiunits = multiunit.empty();
                        singleunits = singleunit.empty();
                        lfps = lfp.empty();
                        for k=this_channels
                            count_multiunits = count_multiunits + 1;
                            fprintf(['... reading SPK channel ' num2str(k) '...  ']);
                            units = getunits([prs.filepath f_plx],this_units,k);
                            nunits = length(units);
                            % add multiunit
                            fprintf(['... adding multiunit ' num2str(count_multiunits)]);
                            multiunits(end+1) = multiunit(this_monk_id,this_session_id,k,maxchannels,this_coord,units(1).spkwf);
                            multiunits(end).add_stim(stim,tfix,tstim,prs);
                            multiunits(end).add_spks(units(1),tstim,ntrls,prs);
                            for l=2:nunits
                                count_singleunits = count_singleunits + 1;
                                % add singleunit
                                fprintf(['... adding singleunit ' num2str(count_singleunits)]);
                                singleunits(end+1) = singleunit(this_monk_id,this_session_id,k,maxchannels,this_coord,units(l).spkwf);
                                singleunits(end).add_stim(stim,tfix,tstim,prs);
                                singleunits(end).add_spks(units(l),tstim,ntrls,prs);
                            end
                            % add lfp
                            count_lfps = count_lfps + 1;
                            fprintf(['\n ... reading AD channel ' num2str(k) '... ']);
                            this_lfp = getlfp([prs.filepath f_plx],k);
                            fprintf(['... adding lfp ' num2str(count_lfps) '\n']);
                            lfps(end+1) = lfp(this_monk_id,this_session_id,k,maxchannels,this_coord);
                            lfps(end).add_stim(stim,tfix,tstim,prs);
                            lfps(end).add_lfps(this_lfp,tstim,ntrls,prs);
                        end
                        plottrials(lfps,multiunits,singleunits,prs);
                    end
            end
        end
        %% function to plot single/multiunit/lfp data from individual electrodes (populations)
        function plotpopulation(this,plottype,unit_type,pop_type,pop_num,pop_name)
            unitnames = {this.populations.unitname};
            groupnames = {this.populations.groupname};
            if strcmp(plottype,'csd')
                x=1;
            elseif (strcmp(plottype,'allpsth') || strcmp(plottype(1:4),'comp')) && ~strcmp(unit_type,'lfps')
                pop{1} = this.(unit_type);
                if strcmp(pop_num,'all')
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,'all');
                elseif ~isempty(pop_num)
                    indx = pop_num;
                elseif isempty(pop_num)
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,pop_name);
                end
                pop{2} = this.populations(indx);
            elseif (strcmp(plottype,'allpsth_1DAzi') || strcmp(plottype(1:4),'comp')) && ~strcmp(unit_type,'lfps')
                pop{1} = this.(unit_type);
                if strcmp(pop_num,'all')
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,'all');
                elseif ~isempty(pop_num)
                    indx = pop_num;
                elseif isempty(pop_num)
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,pop_name);
                end
                pop{2} = this.populations(indx);
            elseif strcmp(plottype,'allpsth') && strcmp(unit_type,'lfps')
                pop{1} = [this.sessions];
                if strcmp(pop_num,'all')
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,'all');
                elseif ~isempty(pop_num)
                    indx = pop_num;
                elseif isempty(pop_num)
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,pop_name);
                end
                pop{2} = this.populations(indx);
            elseif strcmp(plottype,'allpsth_1DAzi') && strcmp(unit_type,'lfps')
                pop{1} = [this.sessions];
                if strcmp(pop_num,'all')
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,'all');
                elseif ~isempty(pop_num)
                    indx = pop_num;
                elseif isempty(pop_num)
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,pop_name);
                end
                pop{2} = this.populations(indx);
            elseif strcmp(plottype(1:4),'clus')
                monk_id = [this.singleunits.monk_id];
                session_id = [this.singleunits.session_id];
                channel_no = [this.singleunits.channel_no];
                monk_id2 = [this.multiunits.monk_id];
                session_id2 = [this.multiunits.session_id];
                channel_no2 = [this.multiunits.channel_no];
                for i=1:length(monk_id)
                    these_multiunits(i) = this.multiunits(monk_id2==monk_id(i) & ...
                        session_id2==session_id(i) & channel_no2==channel_no(i));
                end
                pop{1} = these_multiunits;
                pop{2} = this.singleunits;
                pop{3} = this.sessions;
            elseif ~strcmp(plottype(1:4),'corr') || strcmp(unit_type,'singleunits')
                if strcmp(pop_num,'all')
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,'all');
                elseif ~isempty(pop_num)
                    indx = pop_num;
                elseif isempty(pop_num)
                    indx = strcmp(unitnames,unit_type) & strcmp(groupnames,pop_name);
                end
                pop = this.populations(indx);
                %                pop_multiunits=this.multiunits;
                %                pop_singleunits=this.singleunits;
            else
                count = 0;
                for i=1:length(this.sessions)
                    if ~isempty(this.sessions(i).(unit_type))
                        count = count + 1;
                        pop(count) = this.sessions(i).(unit_type).(plottype);
                    end
                end
            end
            plotpop(this.name,unit_type,pop,pop_type,plottype);
        end
    end
end
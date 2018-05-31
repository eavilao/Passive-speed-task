classdef singleunit < handle
    %%
    properties
        monk_id
        session_id
        channel_no
        coord
        spkwf
        stim
        spks
        tspk
        null = null.empty();
        ves = vestibular.empty();
        vis = visual.empty();
        com = combined.empty();
    end
    %%
    methods
        %% class constructor
        function this = singleunit(monk_id,session_id,channel_no,maxchannels,coord,spkwf)            
            this.monk_id = monk_id;
            this.session_id = session_id;
            this.channel_no = channel_no;
            this.get_coord(maxchannels,coord);
            this.spkwf = mean_spkwf(spkwf);
        end
         %% add row, col, depth of recording
        function get_coord(this,maxchannels,coord)
            this.coord = coord;
            this.coord.depth = coord.depth + (this.channel_no-maxchannels)*100;
        end
        %% add stim values
        function add_stim(this,stim,tfix,tstim,prs)
            this.stim = addstim(stim,tfix,tstim,prs);
        end
        %% add spike times
        function add_spks(this,unit,tstim,ntrls,prs)
            this.spks = addspks(unit,tstim,ntrls,prs);
            this.tspk = unit.tspk;
        end
        %% analyse spikes
        function analyse_spks(this,exp_name,prs)
            switch exp_name
                case 'linearspeed'
                    this.null = analysespks(this,'null',prs);
                    this.ves = analysespks(this,'ves',prs);
                    this.vis = analysespks(this,'vis',prs);
                    this.com = analysespks(this,'com',prs);
                case 'angularspeed'
                    this.null = analysespks(this,'null',prs);
                    this.ves = analysespks(this,'ves_nofix',prs);
                    this.vis = analysespks(this,'vis',prs);
                case '1DAzi'
                    this.null = analysespks_1DAzi(this,'null',prs);
                    this.ves = analysespks_1DAzi(this,'ves',prs);
                    this.vis = analysespks_1DAzi(this,'vis',prs);
                    this.com = analysespks_1DAzi(this,'com',prs);
                case 'HD'
                    this.ves = analysespks_HD(this,'ves',prs);
                    this.vis = analysespks_HD(this,'vis',prs);
                    this.com = analysespks_HD(this,'com',prs);
            end
        end
        %% compute psychometric
        function psycho_neuro_metric(this,exp_name,prs)
            switch exp_name.name
                case 'HD'
                    this.psycho_neuro = compute_psychometric(this)
                    save('psycho_neuro', 'psycho_neuro');
            end
        end
        %% destroy spike times
        function destroy_spks(this)
            this.spks = [];
            this.tspk = [];
        end
        %% destroy stim
        function destroy_stim(this)
            this.stim = [];
        end
        %% plot
        function plot(this,unit_num,exp_name,plottype)
            plotunit(this,unit_num,exp_name,plottype);
        end
    end
end
classdef lfp < handle
    %%
    properties
        monk_id
        session_id
        channel_no
        coord
        stim
        wave
        null = null.empty();
        ves = vestibular.empty();
        vis = visual.empty();
        com = combined.empty();
    end
    %%
    methods
        %% class constructor
        function this = lfp(monk_id,session_id,channel_no,maxchannels,coord)
            this.monk_id = monk_id;
            this.session_id = session_id;
            this.channel_no = channel_no;
            this.get_coord(maxchannels,coord);
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
        function add_lfps(this,lfp,tstim,ntrls,prs)
            this.wave = addlfps(lfp,tstim,ntrls,prs);
        end
        %% analyse lfp
        function analyse_lfps(this,exp_name,prs)
            switch exp_name
                case 'linearspeed'
                    this.null = analyselfps(this,'null',prs);
                    this.ves = analyselfps(this,'ves',prs);
                    this.vis = analyselfps(this,'vis',prs);
                    this.com = analyselfps(this,'com',prs);
                case 'angularspeed'
                    this.null = analyselfps(this,'null',prs);
                    this.ves = analyselfps(this,'ves_nofix',prs);
                    this.vis = analyselfps(this,'vis',prs);
                case '1DAzi'
                    this.null = analyselfps(this,'null',prs);
                    this.ves = analyselfps(this,'ves',prs);
                    this.vis = analyselfps(this,'vis',prs);
                    this.com = analyselfps(this,'com',prs);
                case 'HD'
                    this.null = analyselfps(this,'null',prs);
                    this.ves = analyselfps(this,'ves',prs);
                    this.vis = analyselfps(this,'vis',prs);
                    this.com = analyselfps(this,'com',prs);
            end
        end
        %% destroy lfps
        function destroy_wave(this)
            this.wave = [];
        end
        %% destroy stim
        function destroy_stim(this)
            this.stim = [];
        end
        %% plot
        function plot(this,lfp_num,exp_name,plottype)
            plotlfp(this,lfp_num,exp_name,plottype);
        end
    end
end
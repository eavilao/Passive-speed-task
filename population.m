classdef population < handle
    %%
    properties
        groupname        % all/monk X/narrow-spiking etc.
        unitname         % single/multiunit/lfp
        ves
        vis
        com
        all
    end
    
    %%
    methods
        %% class constructor
        function this = population(unitname)
            this.unitname = unitname;
        end
        %% function to add sessions
        function analyse_population(this,expname,data,unitname,prs)
            switch expname
                case 'linearspeed'
                    this.ves = analysepopulation(data,unitname,'ves',prs);
                    this.vis = analysepopulation(data,unitname,'vis',prs);
                    this.com = analysepopulation(data,unitname,'com',prs);
                    if length(data)==2
                        this.all = analysesession(data);
                    end
                case 'angularspeed'
                    this.ves = analysepopulation(data,unitname,'ves',prs);
                    this.vis = analysepopulation(data,unitname,'vis',prs);
                    if length(data)==2
                        this.all = analysesession(data);
                    end
                case '1DAzi'
                    this.ves = analysepopulation_1DAzi(data,unitname,'ves',prs);
                    this.vis = analysepopulation_1DAzi(data,unitname,'vis',prs);
                    this.com = analysepopulation_1DAzi(data,unitname,'com',prs);
                    if length(data)==2
                        this.all = analysesession(data);
                    end
                case 'HD'
                    this.ves = analysepopulation(data,unitname,'ves',prs);
                    this.vis = analysepopulation(data,unitname,'vis',prs);
                    this.com = analysepopulation(data,unitname,'com',prs);
                    if length(data)==2
                        this.all = analysesession(data);
                    end
            end
        end
    end
end
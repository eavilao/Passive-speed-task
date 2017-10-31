classdef session < handle
    %%
    properties
        monk_id
        session_id
        coord
        multiunits
        singleunits
        lfps
        csd
    end
    %%
    methods
        %% class constructor
        function this = session(monk_id,session_id,coord)
            this.monk_id = monk_id;
            this.session_id = session_id;
            this.coord.row = coord.row; this.coord.column = coord.column;
        end
        %% analyse multiunits
        function analyse_multiunits(this,exp_name,multiunits)
            stim = multiunits(1).stim;
            for i=1:length(multiunits) 
                spkwf(i,:) = multiunits(i).spkwf;
                nspk(i,:) = [multiunits(i).spks.nspk];
            end
            for i=1:length(multiunits)
                for j=1:length(multiunits)
                    nspk1 = [multiunits(i).spks.nspk];
                    nspk2 = [multiunits(j).spks.nspk];
                    for k=1:length(multiunits(i).spks)
                        trls(k).tspk1 = [multiunits(i).spks(k).tspk];
                        trls(k).tspk2 = [multiunits(j).spks(k).tspk];
                    end
                    % spike-count correlation (all trials)
                    [this.multiunits.corr_spatial.r_nspk(i,j),this.multiunits.corr_spatial.p_nspk(i,j)] = ...
                        nspkcorr(nspk1(:),nspk2(:));
                    % temporal correlation (all trials)
                    [this.multiunits.corr_spatial.r_tspk(i,j),this.multiunits.corr_spatial.p_tspk(i,j)] = ...
                        tspkcorr(trls,'true');
                    % signal & noise correlation
                    switch exp_name
                        case 'linearspeed'
                            [this.multiunits.corr_spatial.ves.r_sig(i,j),this.multiunits.corr_spatial.ves.p_sig(i,j),...
                                this.multiunits.corr_spatial.ves.r_noise(i,j),this.multiunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.multiunits.corr_spatial.vis.r_sig(i,j),this.multiunits.corr_spatial.vis.p_sig(i,j),...
                                this.multiunits.corr_spatial.vis.r_noise(i,j),this.multiunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.multiunits.corr_spatial.com.r_sig(i,j),this.multiunits.corr_spatial.com.p_sig(i,j),...
                                this.multiunits.corr_spatial.com.r_noise(i,j),this.multiunits.corr_spatial.com.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                        case 'angularspeed'
                            [this.multiunits.corr_spatial.ves.r_sig(i,j),this.multiunits.corr_spatial.ves.p_sig(i,j),...
                                this.multiunits.corr_spatial.ves.r_noise(i,j),this.multiunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves_nofix');
                            [this.multiunits.corr_spatial.vis.r_sig(i,j),this.multiunits.corr_spatial.vis.p_sig(i,j),...
                                this.multiunits.corr_spatial.vis.r_noise(i,j),this.multiunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                        case '1DAzi'
                            [this.multiunits.corr_spatial.ves.r_sig(i,j),this.multiunits.corr_spatial.ves.p_sig(i,j),...
                                this.multiunits.corr_spatial.ves.r_noise(i,j),this.multiunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.multiunits.corr_spatial.vis.r_sig(i,j),this.multiunits.corr_spatial.vis.p_sig(i,j),...
                                this.multiunits.corr_spatial.vis.r_noise(i,j),this.multiunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.multiunits.corr_spatial.com.r_sig(i,j),this.multiunits.corr_spatial.com.p_sig(i,j),...
                                this.multiunits.corr_spatial.com.r_noise(i,j),this.multiunits.corr_spatial.com.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                        case 'HD'
                            [this.multiunits.corr_spatial.ves.r_sig(i,j),this.multiunits.corr_spatial.ves.p_sig(i,j),...
                                this.multiunits.corr_spatial.ves.r_noise(i,j),this.multiunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.multiunits.corr_spatial.vis.r_sig(i,j),this.multiunits.corr_spatial.vis.p_sig(i,j),...
                                this.multiunits.corr_spatial.vis.r_noise(i,j),this.multiunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.multiunits.corr_spatial.com.r_sig(i,j),this.multiunits.corr_spatial.com.p_sig(i,j),...
                                this.multiunits.corr_spatial.com.r_noise(i,j),this.multiunits.corr_spatial.com.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                    end
                    % spike-train cross correlogram to detect duplicates
                    tspk1 = multiunits(i).tspk;
                    tspk2 = multiunits(j).tspk;
                    this.multiunits.corr_spatiotemporal.r_tspk(i,j,:) = spiketrains_xcorr(tspk1,tspk2);
                end
            end
            this.multiunits.corr_spatial.x = 1:length(multiunits);
            this.multiunits.corr_distance = corr_ch2dist(this.multiunits.corr_spatial,0.1);
            this.multiunits.corr_spatiotemporal.x = 1:length(multiunits);
            this.multiunits.corr_spatiotemporal.t = linspace(-0.2,0.2,401);
            switch exp_name
                case 'linearspeed'
                    this.multiunits.ves = analysecorr(stim,nspk,'ves');
                    this.multiunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves');
                    this.multiunits.vis = analysecorr(stim,nspk,'vis');
                    this.multiunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                    this.multiunits.com = analysecorr(stim,nspk,'com');
                    this.multiunits.com.fisherinfo = analyseinfo(stim,nspk,'com');
                case 'angularspeed'
                    this.multiunits.ves = analysecorr(stim,nspk,'ves_nofix');
                    this.multiunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves_nofix');
                    this.multiunits.vis = analysecorr(stim,nspk,'vis');
                    this.multiunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                case '1DAzi'
                    this.multiunits.ves = analysecorr(stim,nspk,'ves');
                    this.multiunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves');
                    this.multiunits.vis = analysecorr(stim,nspk,'vis');
                    this.multiunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                    this.multiunits.com = analysecorr(stim,nspk,'com');
                    this.multiunits.com.fisherinfo = analyseinfo(stim,nspk,'com');
            end
        end
        %% analyse singleunits
        function analyse_singleunits(this,exp_name,singleunits)
            stim = singleunits(1).stim;
            for i=1:length(singleunits)
                spkwf(i,:) = singleunits(i).spkwf;
                nspk(i,:) = [singleunits(i).spks.nspk];
            end
            for i=1:length(singleunits)
                for j=1:length(singleunits)
                    nspk1 = [singleunits(i).spks.nspk];
                    nspk2 = [singleunits(j).spks.nspk];
                    for k=1:length(singleunits(i).spks)
                        trls(k).tspk1 = [singleunits(i).spks(k).tspk];
                        trls(k).tspk2 = [singleunits(j).spks(k).tspk];
                    end
                    % spike-count correlation (all trials)
                    [this.singleunits.corr_spatial.r_nspk(i,j),this.singleunits.corr_spatial.p_nspk(i,j)] = ...
                        nspkcorr(nspk1(:),nspk2(:));
                    % temporal correlation (all trials)
                    [this.singleunits.corr_spatial.r_tspk(i,j),this.singleunits.corr_spatial.p_tspk(i,j)] = ...
                        tspkcorr(trls,'true');
                    % signal & noise correlation
                    switch exp_name
                        case 'linearspeed'
                            [this.singleunits.corr_spatial.ves.r_sig(i,j),this.singleunits.corr_spatial.ves.p_sig(i,j),...
                                this.singleunits.corr_spatial.ves.r_noise(i,j),this.singleunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.singleunits.corr_spatial.vis.r_sig(i,j),this.singleunits.corr_spatial.vis.p_sig(i,j),...
                                this.singleunits.corr_spatial.vis.r_noise(i,j),this.singleunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.corr_spatial.com.r_sig(i,j),this.singleunits.corr_spatial.com.p_sig(i,j),...
                                this.singleunits.corr_spatial.com.r_noise(i,j),this.singleunits.corr_spatial.com.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                        case 'angularspeed'
                            [this.singleunits.corr_spatial.ves.r_sig(i,j),this.singleunits.corr_spatial.ves.p_sig(i,j),...
                                this.singleunits.corr_spatial.ves.r_noise(i,j),this.singleunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves_nofix');
                            [this.singleunits.corr_spatial.vis.r_sig(i,j),this.singleunits.corr_spatial.vis.p_sig(i,j),...
                                this.singleunits.corr_spatial.vis.r_noise(i,j),this.singleunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                        case '1DAzi'
                            [this.singleunits.corr_spatial.ves.r_sig(i,j),this.singleunits.corr_spatial.ves.p_sig(i,j),...
                                this.singleunits.corr_spatial.ves.r_noise(i,j),this.singleunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.singleunits.corr_spatial.vis.r_sig(i,j),this.singleunits.corr_spatial.vis.p_sig(i,j),...
                                this.singleunits.corr_spatial.vis.r_noise(i,j),this.singleunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.corr_spatial.com.r_sig(i,j),this.singleunits.corr_spatial.com.p_sig(i,j),...
                                this.singleunits.corr_spatial.com.r_noise(i,j),this.singleunits.corr_spatial.com.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                        case 'HD'
                            [this.singleunits.corr_spatial.ves.r_sig(i,j),this.singleunits.corr_spatial.ves.p_sig(i,j),...
                                this.singleunits.corr_spatial.ves.r_noise(i,j),this.singleunits.corr_spatial.ves.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.singleunits.corr_spatial.vis.r_sig(i,j),this.singleunits.corr_spatial.vis.p_sig(i,j),...
                                this.singleunits.corr_spatial.vis.r_noise(i,j),this.singleunits.corr_spatial.vis.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.corr_spatial.com.r_sig(i,j),this.singleunits.corr_spatial.com.p_sig(i,j),...
                                this.singleunits.corr_spatial.com.r_noise(i,j),this.singleunits.corr_spatial.com.p_noise(i,j)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                    end
                    % spike-train cross correlogram to detect duplicates
                    tspk1 = singleunits(i).tspk;
                    tspk2 = singleunits(j).tspk;
                    this.singleunits.corr_spatiotemporal.r(i,j,:) = spiketrains_xcorr(tspk1,tspk2);
                end
            end
            this.singleunits.corr_spatial.x =  [singleunits.channel_no];
            this.singleunits.corr_distance = corr_ch2dist(this.singleunits.corr_spatial,0.1);
            this.singleunits.corr_spatiotemporal.x =  [singleunits.channel_no];
            this.singleunits.corr_spatiotemporal.t = linspace(-0.2,0.2,401);
            switch exp_name
                case 'linearspeed'
                    this.singleunits.ves = analysecorr(stim,nspk,'ves');
                    this.singleunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves');
                    this.singleunits.vis = analysecorr(stim,nspk,'vis');
                    this.singleunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                    this.singleunits.com = analysecorr(stim,nspk,'com');
                    this.singleunits.com.fisherinfo = analyseinfo(stim,nspk,'com');
                case 'angularspeed'
                    this.singleunits.ves = analysecorr(stim,nspk,'ves_nofix');
                    this.singleunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves_nofix');
                    this.singleunits.vis = analysecorr(stim,nspk,'vis');
                    this.singleunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                case '1DAzi'
                    this.singleunits.ves = analysecorr(stim,nspk,'ves');
                    this.singleunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves');
                    this.singleunits.vis = analysecorr(stim,nspk,'vis');
                    this.singleunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                    this.singleunits.com = analysecorr(stim,nspk,'com');
                    this.singleunits.com.fisherinfo = analyseinfo(stim,nspk,'com');
                case 'HD'
                    this.singleunits.ves = analysecorr(stim,nspk,'ves');
                    this.singleunits.ves.fisherinfo = analyseinfo(stim,nspk,'ves');
                    this.singleunits.vis = analysecorr(stim,nspk,'vis');
                    this.singleunits.vis.fisherinfo = analyseinfo(stim,nspk,'vis');
                    this.singleunits.com = analysecorr(stim,nspk,'com');
                    this.singleunits.com.fisherinfo = analyseinfo(stim,nspk,'com');
            end
        end
        %% analyse clustering_true
        function analyse_clustering(this,exp_name,singleunits,multiunits)
            stim = singleunits(1).stim;
            for i=1:length(singleunits)
                channel_no = singleunits(i).channel_no; 
                randchannel_no = randperm(16); randchannel_no = randchannel_no(1);
                nspk1 = [singleunits(i).spks.nspk];
                nspk2 = [multiunits(channel_no).spks.nspk];
                nspk_rand = [multiunits(randchannel_no).spks.nspk];
                for k=1:length(singleunits(i).spks)
                    trls(k).tspk1 = [singleunits(i).spks(k).tspk];
                    trls(k).tspk2 = [multiunits(channel_no).spks(k).tspk];
                    trls(k).tspk_rand = [multiunits(randchannel_no).spks(k).tspk];
                end
                % spike-count correlation (all trials)
                [this.singleunits.clustering.r_nspk(i),this.singleunits.clustering.p_nspk(i)] = ...
                    nspkcorr(nspk1(:),nspk2(:));
                [this.singleunits.clustering_shuffled.r_nspk(i),this.singleunits.clustering_shuffled.p_nspk(i)] = ...
                    nspkcorr(nspk1(:),nspk_rand(:));
                % temporal correlation (all trials)
                [this.singleunits.clustering.r_tspk(i),this.singleunits.clustering.p_tspk(i)] = ...
                    tspkcorr(trls,'true');
                [this.singleunits.clustering_shuffled.r_tspk(i),this.singleunits.clustering_shuffled.p_tspk(i)] = ...
                    tspkcorr(trls,'shuffled');
                    % signal & noise correlation
                    switch exp_name
                        case 'linearspeed'
                            [this.singleunits.clustering.ves.r_sig(i),this.singleunits.clustering.ves.p_sig(i),...
                                this.singleunits.clustering.ves.r_noise(i),this.singleunits.clustering.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.singleunits.clustering_shuffled.ves.r_sig(i),this.singleunits.clustering_shuffled.ves.p_sig(i),...
                                this.singleunits.clustering_shuffled.ves.r_noise(i),this.singleunits.clustering_shuffled.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'ves');
                            
                            [this.singleunits.clustering.vis.r_sig(i),this.singleunits.clustering.vis.p_sig(i),...
                                this.singleunits.clustering.vis.r_noise(i),this.singleunits.clustering.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.clustering_shuffled.vis.r_sig(i),this.singleunits.clustering_shuffled.vis.p_sig(i),...
                                this.singleunits.clustering_shuffled.vis.r_noise(i),this.singleunits.clustering_shuffled.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'vis');
                            
                            [this.singleunits.clustering.com.r_sig(i),this.singleunits.clustering.com.p_sig(i),...
                                this.singleunits.clustering.com.r_noise(i),this.singleunits.clustering.com.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                            [this.singleunits.clustering_shuffled.com.r_sig(i),this.singleunits.clustering_shuffled.com.p_sig(i),...
                                this.singleunits.clustering_shuffled.com.r_noise(i),this.singleunits.clustering_shuffled.com.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'com');
                        case 'angularspeed'
                            [this.singleunits.clustering.ves.r_sig(i),this.singleunits.clustering.ves.p_sig(i),...
                                this.singleunits.clustering.ves.r_noise(i),this.singleunits.clustering.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves_nofix');
                            [this.singleunits.clustering_shuffled.ves.r_sig(i),this.singleunits.clustering_shuffled.ves.p_sig(i),...
                                this.singleunits.clustering_shuffled.ves.r_noise(i),this.singleunits.clustering_shuffled.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'ves_nofix');
                            
                            [this.singleunits.clustering.vis.r_sig(i),this.singleunits.clustering.vis.p_sig(i),...
                                this.singleunits.clustering.vis.r_noise(i),this.singleunits.clustering.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.clustering_shuffled.vis.r_sig(i),this.singleunits.clustering_shuffled.vis.p_sig(i),...
                                this.singleunits.clustering_shuffled.vis.r_noise(i),this.singleunits.clustering_shuffled.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'vis');
                             case '1DAzi'
                            [this.singleunits.clustering.ves.r_sig(i),this.singleunits.clustering.ves.p_sig(i),...
                                this.singleunits.clustering.ves.r_noise(i),this.singleunits.clustering.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.singleunits.clustering_shuffled.ves.r_sig(i),this.singleunits.clustering_shuffled.ves.p_sig(i),...
                                this.singleunits.clustering_shuffled.ves.r_noise(i),this.singleunits.clustering_shuffled.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'ves');
                            
                            [this.singleunits.clustering.vis.r_sig(i),this.singleunits.clustering.vis.p_sig(i),...
                                this.singleunits.clustering.vis.r_noise(i),this.singleunits.clustering.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.clustering_shuffled.vis.r_sig(i),this.singleunits.clustering_shuffled.vis.p_sig(i),...
                                this.singleunits.clustering_shuffled.vis.r_noise(i),this.singleunits.clustering_shuffled.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'vis');
                            
                            [this.singleunits.clustering.com.r_sig(i),this.singleunits.clustering.com.p_sig(i),...
                                this.singleunits.clustering.com.r_noise(i),this.singleunits.clustering.com.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                            [this.singleunits.clustering_shuffled.com.r_sig(i),this.singleunits.clustering_shuffled.com.p_sig(i),...
                                this.singleunits.clustering_shuffled.com.r_noise(i),this.singleunits.clustering_shuffled.com.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'com');
                            
                        case 'HD'
                            [this.singleunits.clustering.ves.r_sig(i),this.singleunits.clustering.ves.p_sig(i),...
                                this.singleunits.clustering.ves.r_noise(i),this.singleunits.clustering.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'ves');
                            [this.singleunits.clustering_shuffled.ves.r_sig(i),this.singleunits.clustering_shuffled.ves.p_sig(i),...
                                this.singleunits.clustering_shuffled.ves.r_noise(i),this.singleunits.clustering_shuffled.ves.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'ves');
                            
                            [this.singleunits.clustering.vis.r_sig(i),this.singleunits.clustering.vis.p_sig(i),...
                                this.singleunits.clustering.vis.r_noise(i),this.singleunits.clustering.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'vis');
                            [this.singleunits.clustering_shuffled.vis.r_sig(i),this.singleunits.clustering_shuffled.vis.p_sig(i),...
                                this.singleunits.clustering_shuffled.vis.r_noise(i),this.singleunits.clustering_shuffled.vis.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'vis');
                            
                            [this.singleunits.clustering.com.r_sig(i),this.singleunits.clustering.com.p_sig(i),...
                                this.singleunits.clustering.com.r_noise(i),this.singleunits.clustering.com.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk2(:),'com');
                            [this.singleunits.clustering_shuffled.com.r_sig(i),this.singleunits.clustering_shuffled.com.p_sig(i),...
                                this.singleunits.clustering_shuffled.com.r_noise(i),this.singleunits.clustering_shuffled.com.p_noise(i)] = ...
                                signoisecorr(stim,nspk1(:),nspk_rand(:),'com');
                            
                            
                    end
            end
        end
        %% analyse lfps
        function analyse_lfps(this,exp_name,lfps)
            this.lfps.time = lfps(1).ves.time;
            dt = diff(this.lfps.time); dt = dt(1);
            for i=1:length(lfps)
                wave(i,:) = reshape([lfps(i).wave.v],[numel([lfps(i).wave.v]),1]);
                switch exp_name
                    case 'linearspeed'
                        wave_null(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==-1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==-1).v]),1]);
                        wave_ves(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==1).v]),1]);
                        wave_vis(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==2).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==2).v]),1]);
                    case 'angularspeed'
                        wave_null(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==-1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==-1).v]),1]);
                        wave_ves(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==0).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==0).v]),1]);
                        wave_vis(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==2).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==2).v]),1]);
                    case '1DAzi'
                        wave_null(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==-1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==-1).v]),1]);
                        wave_ves(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==1).v]),1]);
                        wave_vis(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==2).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==2).v]),1]);
                    case 'HD'
                        wave_null(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==-1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==-1).v]),1]);
                        wave_ves(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==1).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==1).v]),1]);
                        wave_vis(i,:) = reshape([lfps(i).wave(lfps(i).stim.modality==2).v],...
                            [numel([lfps(i).wave(lfps(i).stim.modality==2).v]),1]);
                end
            end
            this.lfps.corr_spatial.r = corr(wave');
            this.lfps.corr_spatial.null.r = corr(wave_null');
            this.lfps.corr_spatial.ves.r = corr(wave_ves');
            this.lfps.corr_spatial.vis.r = corr(wave_vis');
            this.lfps.corr_spatial.x = 1:length(lfps);
            for i=1:length(lfps)
                for j=1:length(lfps)
                    wave1 = wave(i,:);
                    wave2 = wave(j,:);
                    this.lfps.corr_spatiotemporal.r(i,j,:) = xcorr(wave1,wave2,round(1/dt));
                end
            end
            this.lfps.corr_spatiotemporal.x = 1:length(lfps);
            this.lfps.corr_distance = corr_ch2dist(this.lfps.corr_spatial,0.1);
            this.lfps.corr_spatiotemporal.t = linspace(-1,1,2*round(1/dt)+1);
            clear wave_ves wave_vis;
            switch exp_name
                case 'linearspeed'
                    for i=1:length(lfps)
                        wave_ves(i,:) = mean(lfps(i).ves.wave_pst);
                        wave_vis(i,:) = mean(lfps(i).vis.wave_pst);
                        wave_com(i,:) = mean(lfps(i).com.wave_pst);
                    end
                    this.csd.ves = computecsd(wave_ves,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.ves = mean(wave_ves);
                    this.csd.vis = computecsd(wave_vis,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.vis = mean(wave_vis);
                    this.csd.com = computecsd(wave_com,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.com = mean(wave_com);
                case 'angularspeed'
                    for i=1:length(lfps)
                        wave_ves(i,:) = mean(lfps(i).ves.wave_pst);
                        wave_vis(i,:) = mean(lfps(i).vis.wave_pst);
                    end
                    this.csd.ves = computecsd(wave_ves,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.ves = mean(wave_ves);
                    this.csd.vis = computecsd(wave_vis,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.vis = mean(wave_vis);
                     case '1DAzi'
                    for i=1:length(lfps)
                        wave_ves(i,:) = mean(lfps(i).ves.wave_pst);
                        wave_vis(i,:) = mean(lfps(i).vis.wave_pst);
                        wave_com(i,:) = mean(lfps(i).com.wave_pst);
                    end
                    this.csd.ves = computecsd(wave_ves,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.ves = mean(wave_ves);
                    this.csd.vis = computecsd(wave_vis,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.vis = mean(wave_vis);
                    this.csd.com = computecsd(wave_com,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.com = mean(wave_com);
                     case 'HD'
                    for i=1:length(lfps)
                        wave_ves(i,:) = mean(lfps(i).ves.wave_pst);
                        wave_vis(i,:) = mean(lfps(i).vis.wave_pst);
                        wave_com(i,:) = mean(lfps(i).com.wave_pst);
                    end
                    this.csd.ves = computecsd(wave_ves,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.ves = mean(wave_ves);
                    this.csd.vis = computecsd(wave_vis,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.vis = mean(wave_vis);
                    this.csd.com = computecsd(wave_com,this.lfps.time,0.1); %dx = 0.1 mm
                    this.lfps.com = mean(wave_com);
            end
        end
        %% plot
        function plot(this,these_multiunits,plottype)
            plotsession(this,these_multiunits,plottype);
        end
    end
end
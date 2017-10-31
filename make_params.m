function params = make_params

%% params for spike-field coherency estimation
params.pad=0;                       % pad data to next highest power of 2
params.tapers=[4 7];                % use time-bandwidth(TW) product of 4, T being (tend-tbeg); 7 tapers
params.trialave=1;                  % get a single SFC estimate across all trials 
params.fpass=1:200;                 % estimate for all frequencies between 1 and 200 Hz
params.err=0;                       % no error bars for coherence measures
params.tbeg=0.4;                    % tbeg=0.4, tend=1; therefore TW=4 => W=4/0.6 = 6.67 Hz of spectral smoothing
params.tend=1;
params.tfix=0.3;                    % duration of fixation
params.gsmooth=0.00625;               % width of gaussian for estimating psth
params.movingwin = [0.3 0.002];
params.tf_tbeg=-0.25;
params.tf_tend=1.5;

%% optional fields for pre-processing
params.interplfp=1;                 % interpolates lfp around spikes to remove spurious spike-LFP correlations
params.interpwin = 4e-3;            % interpolation window of 4ms
params.rmvtrls_saturatedlfps = 1;   % remove trials with saturated lfp
params.threshduration = 4e-3;       % lfp constant for more than 4ms => saturation

%% optional fields to append to existing singleunits
params.computetfspectra = 0;
params.checksfcsignificance = 0;
params.checksfcmodsignificance = 0;
params.checklfpmodsignificance = 0;
params.permtestsfc = 0;
params.computesfcmodulation = 0;
params.computepsfc = 1;
params.computespikephase = 0;
params.extractspikephase = 0;
params.lfp_filt = [33 43];
params.fitlfp = 0;
params.runsimulation = 0;
params.runsimulation_filtdata = 0;
params.recomputespectra = 0;
params.computelfpstats = 0;
params.loadspikewidth = 0;
params.computegammashifting = 0;
params.phasematch = 0;
params.fitAGmodel = 0;
params.phasematched_spectra = 0;

%% params for spike-LFP model - simulation_I
params.exponent = 1;
params.gain = 1;
params.numthresh = 20;
params.numlambdain = 20;
params.lambdamax = 2500;

%% params for spike-LFP model - simulation_II
params.exponent2 = 1;
params.gain2 = 1;
params.npts = 41;
params.ntrl = 200;
params.trllen = 0.6;

%% params for AGmodel
params.Amax = 100;
params.Gmax = 100;

%% params for SFC significance tests
params.fmin_gamma = 28;
params.fmax_gamma = 50;
params.bwmin_gamma = 16;
params.bwmin_gamma_mod = 7;
params.pcrit = 0.005;
params.pcrit_mod = 0.05;
params.min_trls = 30;

%% about data
cond{1}=[1 3]; cond{2}=[2 4];
cond{3}=[5 7]; cond{4}=[6 8];
params.cond = cond;
params.numsessions_d98 = 70;
params.numsessions_f03 = 142;
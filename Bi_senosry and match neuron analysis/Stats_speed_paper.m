% Script to run some statistical analysis for the 7a speed paper
% All are times from start of static motion period, so need to subtract 0.4
% load experiments.mat


%% Linear Speed
%% Extract latency values 
latency_ves = experiments(1).populations(1).ves.all.latency.t_on-0.4;
latency_vis = experiments(1).populations(1).vis.all.latency.t_on-0.4;
latency_com = experiments(1).populations(1).com.all.latency.t_on-0.4;

%% Extract rise time 
rise_ves = experiments(1).populations(1).ves.all.rise.t_on-0.4; 
rise_vis = experiments(1).populations(1).vis.all.rise.t_on-0.4; 
rise_com = experiments(1).populations(1).com.all.rise.t_on-0.4; 
%% Extract peak responses
peak_ves = experiments(1).populations(1).ves.all.peak.t_on-0.4;
peak_vis = experiments(1).populations(1).vis.all.peak.t_on-0.4;
peak_com = experiments(1).populations(1).com.all.peak.t_on-0.4; 


%% Compute ttest
[h_ves_vis,p_ves_vis]=ttest(rise_ves, rise_vis); 

%% compute correlations for rise times
% [r_ves_com,p_ves_com]=corr(rise_ves(~isnan(rise_ves) & ~isnan(rise_com)),rise_com(~isnan(rise_ves) & ~isnan(rise_com))); 
% [r_vis_com,p_ves_com]=corr(rise_vis(~isnan(rise_vis) & ~isnan(rise_com)),rise_com(~isnan(rise_vis) & ~isnan(rise_com))); 

[r_ves_com,p_ves_com]=nancorr(rise_ves(:),rise_com(:)); 
[r_vis_com,p_vis_com]=nancorr(rise_com(:),rise_vis(:)); 

%% SSI
%% Extract discindx for all neurons
ves_discindx = experiments(1).populations(1).ves.all.discindx; 
vis_discindx = experiments(1).populations(1).vis.all.discindx; 
com_discindx = experiments(1).populations(1).com.all.discindx;
% Stat
[hVesVisDiscinx,pVesVisDiscinx] = ttest(ves_discindx,vis_discindx);
[hComVisDiscinx,pComVisDiscinx] = ttest(com_discindx,vis_discindx);
[hComVesDiscinx,pComVesDiscinx] = ttest(com_discindx,ves_discindx);



%% Discindx for responsive neurons
% Vestibular discindx
resp_discindx = [experiments(1).populations(1).ves.resp.discindx experiments(1).populations(1).vis.resp.discindx experiments(1).populations(1).com.resp.discindx]; 
unresp_discinx = [experiments(1).populations(1).ves.unresp.discindx experiments(1).populations(1).vis.unresp.discindx experiments(1).populations(1).com.unresp.discindx]; 
%Stat
[hDiscindx,pDiscindx] = ranksum(resp_discindx,unresp_discinx) 



%% Angular speed
%% Extract latency times 
ang_latency_ves = experiments(2).populations(1).ves.all.latency.t_on-0.4;
ang_latency_vis = experiments(2).populations(1).vis.all.latency.t_on-0.4;
%% Extract rise time 
ang_rise_ves = experiments(2).populations(1).ves.all.rise.t_on-0.4; 
ang_rise_vis = experiments(2).populations(1).vis.all.rise.t_on-0.4; 
%% Extract peak responses
ang_peak_ves = experiments(2).populations(1).ves.all.peak.t_on-0.4;
ang_peak_vis = experiments(2).populations(1).vis.all.peak.t_on-0.4;

%% Compute ttest
[h_ang_ves_vis,p_ang_ves_vis]=ttest(ang_latency_ves, ang_latency_vis); 



%% Calculate tuninig for direction of rotation
%% Extract info to calc corr for fastest speed to CW and CCW on speed tuned neurons

tuned_CW_indx = experiments(2).populations(1).ves.tuned.stim_pref == 45; 
tuned_CCW_indx = experiments(2).populations(1).ves.tuned.stim_pref == -45;

tuned_all = experiments(2).populations(1).ves.tuned.indx; 

for i = 1:length(experiments(2).singleunits)
rate_CCW_ves(i) = experiments(2).singleunits(i).ves.rate_avg(1).mu; 
rate_CW_ves(i) = experiments(2).singleunits(i).ves.rate_avg(10).mu;
rate_CCW_vis(i) = experiments(2).singleunits(i).vis.rate_avg(1).mu; 
rate_CW_vis(i) = experiments(2).singleunits(i).vis.rate_avg(10).mu;
end 

[r_speed_ves,p_speed_ves]= corr(rate_CCW_ves(:),rate_CW_ves(:));
[r_speed_vis,p_speed_vis]= corr(rate_CCW_vis(:),rate_CW_vis(:));


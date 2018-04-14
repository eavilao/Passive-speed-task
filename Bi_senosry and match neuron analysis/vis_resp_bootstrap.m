% Compute variability in visual resp by bootstrapping (100 times)
exp_name= 'linearspeed'
load('vis_spks.mat'); run(default_prs, 'linearspeed', 44); 
ntrls = length(vis_spks(1,:,1)); 

for cellNum = 1:length(vis_spks(1,1,:))
    for boot_count = 1:100;
        for indx_count = randsample(length(vis_spks(1,:,1)),length(vis_spks(1,:,1)),'true')
            
            tspk = {spk.tspk}; t_prestim = prs.tbeg_acc - prs.tstim_on;

            rate_pst = hist(vis_spks{1,indx_count,cellNum},stim_time)/(ntrls*prs.dt); % extract 20 random trials and store
            rate_pst = smooth_pst(results.rate_pst(i,:),dt,prs.tsmooth);

        end
    end
end


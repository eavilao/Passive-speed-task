% % Extract info from peak times using bootstrapping
% nunits = 1:length(experiments(1).singleunits);
% 
% for i=1:length(nunits)
%     exc_ves(i) = experiments(1).singleunits(i).ves.stats.flags.exc;
%     exc_vis(i) = experiments(1).singleunits(i).vis.stats.flags.exc;
% end
% % % k=1:3
% % units_ves = experiments(1).singleunits(exc_ves);
% % units_vis = experiments(1).singleunits(exc_vis);
% % % end
% 
% %per condition
% cond = {'ves','vis', 'com'};
% for i=1:length(cond)
%     for k = 1:length(units)
%         peak_time_mu(k,i)= units(k).(cond{i}).peak_time.mu;
%         peak_time_sem(k,i)= units(k).(cond{i}).peak_time.sem;
%     end
%     figure;hist(peak_time_mu(:,i));
%     title(cond{i})
% end
% 
% mean(peak_time_mu(:,2)) % mean of peak times for vis
% mean(peak_time_sem(:,2))  % mean of sem
% 
% % compute Fano factor (std^2/mean)
% fano_factor = std(peak_time_mu(:,2))^2/mean(peak_time_mu(:,2))
% % coefficient of variation (std/mean)
% 
% cv=mean(peak_time_sem(:,2))/mean(peak_time_mu(:,2))

%% Extract info from peak times using bootstrapping

% extract indices for exc neurons in ves and vis
ves_indx = experiments(1).populations.ves.exc.indx; 
vis_indx = experiments(1).populations.vis.exc.indx; 
com_indx = experiments(1).populations.com.exc.indx;

% extract mu and sem
ves = [experiments(1).singleunits.ves]; peak_ves = [ves.peak_time]; 
vis = [experiments(1).singleunits.vis]; peak_vis = [vis.peak_time];
com = [experiments(1).singleunits.com]; peak_com = [com.peak_time];

ves_mu = [peak_ves.mu]; vis_mu = [peak_vis.mu]; com_mu = [peak_com.mu];
ves_sem = [peak_ves.sem]; vis_sem = [peak_vis.sem]; com_sem = [peak_com.sem];

% compute mean for mu and sem
ves_mu_pop = mean(ves_mu(ves_indx)); 
vis_mu_pop = mean(vis_mu(vis_indx));
com_mu_pop = mean(com_mu(com_indx));

ves_sem_pop = mean(ves_sem(ves_indx)); 
vis_sem_pop = mean(vis_sem(vis_indx)); 
com_sem_pop = mean(com_sem(com_indx)); 

% compute CV
cv_ves = ves_sem(ves_indx)./ves_mu(ves_indx); 
cv_vis = vis_sem(vis_indx)./vis_mu(vis_indx);
cv_com = com_sem(com_indx)./com_mu(com_indx);

grand_cv_ves = mean(cv_ves)
grand_cv_vis = mean(cv_vis)
grand_cv_com = mean(cv_com)

std_err_cv_ves = std(cv_ves)/sqrt(length(cv_ves))
std_err_cv_vis = std(cv_vis)/sqrt(length(cv_vis))
std_err_cv_com = std(cv_com)/sqrt(length(cv_com))
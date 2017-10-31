%% plot noise corr
r_noise_ves = [];
for i=[1:3 5:44], r_noise_ves = [r_noise_ves experiments(1).sessions(i).singleunits.clustering.ves.r_noise]; end
ecdf(r_noise_ves);
[F,X] = ecdf(r_noise_ves); plot(X,F,'Color','r','Linewidth',2);
hold on;

r_noise_vis = [];
for i=[1:3 5:44], r_noise_vis = [r_noise_vis experiments(1).sessions(i).singleunits.clustering.vis.r_noise]; end
[F,X] = ecdf(r_noise_vis); plot(X,F,'Color','g','Linewidth',2);
hold on;

r_noise_com = [];
for i=[1:3 5:44], r_noise_com = [r_noise_com experiments(1).sessions(i).singleunits.clustering.com.r_noise]; end
[F,X] = ecdf(r_noise_com); plot(X,F,'Color','b','Linewidth',2);

axis([-1 1 0 1]);
vline(0,'--k');
set(gca,'TickDir','Out','XTick',[-1 0 1],'YTick',[0 0.5 1]);

%% plot signal corr
figure;
r_sig_ves = [];
for i=[1:3 5:44], r_sig_ves = [r_sig_ves experiments(1).sessions(i).singleunits.clustering.ves.r_sig]; end
ecdf(r_sig_ves);
[F,X] = ecdf(r_sig_ves); plot(X,F,'Color','r','Linewidth',2);
hold on;

r_sig_vis = [];
for i=[1:3 5:44], r_sig_vis = [r_sig_vis experiments(1).sessions(i).singleunits.clustering.vis.r_sig]; end
[F,X] = ecdf(r_sig_vis); plot(X,F,'Color','g','Linewidth',2);
hold on;

r_sig_com = [];
for i=[1:3 5:44], r_sig_com = [r_sig_com experiments(1).sessions(i).singleunits.clustering.com.r_sig]; end
[F,X] = ecdf(r_sig_com); plot(X,F,'Color','b','Linewidth',2);

axis([-1 1 0 1]);
vline(0,'--k');
set(gca,'TickDir','Out','XTick',[-1 0 1],'YTick',[0 0.5 1]);
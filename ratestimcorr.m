indx1 = experiments(2).populations(1).ves.resp.indx;
indx2 = experiments(2).populations(1).ves.tuned.indx;
indx = indx1 & indx2;
for i=find(indx)
    mu = [experiments(2).singleunits(i).ves.rate_avg.mu];
    mu_ccw = mu(1:5);
    mu_cw = mu(6:end);
    [h(i),p_ves(i)] = ttest(mu_cw,mu_ccw);
end

indx1 = experiments(2).populations(1).vis.resp.indx;
indx2 = experiments(2).populations(1).vis.tuned.indx;
for i=find(indx)
    mu = [experiments(2).singleunits(i).vis.rate_avg.mu];
    mu_ccw = mu(1:5);
    mu_cw = mu(6:end);
    [h(i),p_vis(i)] = ttest(mu_cw,mu_ccw);
end

clear mu_cw mu_ccw;
indx = experiments(2).populations(1).ves.tuned.indx;
for i=find(indx)
    mu = [experiments(2).singleunits(i).ves.rate_avg.mu];
    mu_ccw(i) = mu(1);
    mu_cw(i) = mu(end);
end

clear mu_cw mu_ccw;
indx = experiments(2).populations(1).vis.resp.indx;
for i=find(indx)
    mu = [experiments(2).singleunits(i).vis.rate_avg.mu];
    mu_ccw(i) = mu(1);
    mu_cw(i) = mu(end);
end

clear mu_cw mu_ccw;
indx1 = experiments(2).populations(1).ves.tuned.indx;
indx2 = experiments(2).populations(1).ves.exc.indx;
indx = indx1 & indx2;
count = 0;
for i=find(indx)
    count = count+1;
    mu = [experiments(2).singleunits(i).ves.rate_avg.mu];
    mu_ccw(count,:) = mu(1:5); mu_ccw(count,:) = fliplr(mu_ccw(count,:));
    mu_cw(count,:) = mu(6:end);
    
end
mu = [mu_cw; mu_ccw];
[r_ves,p_ves] = corrcoef(repmat(1:5,[52 1]),mu);

clear mu_cw mu_ccw;
indx1 = experiments(2).populations(1).vis.tuned.indx;
indx2 = experiments(2).populations(1).vis.exc.indx;
indx = indx1 & indx2;
count = 0;
for i=find(indx)
    count = count+1;
    mu = [experiments(2).singleunits(i).vis.rate_avg.mu];
    mu_ccw(count,:) = mu(1:5); mu_ccw(count,:) = fliplr(mu_ccw(count,:));
    mu_cw(count,:) = mu(6:end);
    
end
mu = [mu_cw; mu_ccw];
[r_vis,p_vis] = corrcoef(repmat(1:5,[54 1]),mu)
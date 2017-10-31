function sig_pst = pst_sig(spk,t,ntrls,prs)

for i=1:ntrls
    rate_pst(i,:) = hist(spk(i).tspk,t)/(ntrls*prs.dt); % mean psth
    rate_pst(i,:) = smooth_pst(rate_pst(i,:),prs.dt,prs.tsmooth);
end
sig_pst = std(rate_pst);

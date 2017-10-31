function pst = smooth_pst(pst,dt,tsmooth)

fs = 1/dt;
sig = tsmooth*fs; %filter width
sz = fs; %filter size
t = linspace(-sz/2, sz/2, sz);
h = exp(-t.^2/(2*sig^2));
h = h/sum(h);
pst = conv(pst,h,'same');
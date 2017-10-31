function lfp = getlfp(fname,getch)

[fs, ~, ~, ~, lfp.v] = plx_ad_v(fname,getch-1);
t_beg = 0; t_end = length(lfp.v)/fs;
lfp.t = linspace(t_beg,t_end,length(lfp.v));
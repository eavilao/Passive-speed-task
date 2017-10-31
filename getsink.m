function [z_sink,t_sink] = getsink(csd)

[z_sink,t_sink]=find(csd==min(min(csd)));
csd_t = (csd(:,t_sink));
z_switch = find(diff(csd_t<0));
z1 = z_switch(find(z_switch<z_sink, 1, 'last' ));
z2 = z_switch(find(z_switch>z_sink, 1 ));
z_sink = z1+round(sum(csd_t(z1:z2).*(1:(z2-z1+1))')/sum(csd_t(z1:z2)));
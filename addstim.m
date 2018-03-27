function stim = addstim(stim1,tfix,tstim,prs)

% speed = azimuth for 1DAzi

stim.tfix = [tfix.on] - [tstim.on];
stim.tstim = [tstim.off] - [tstim.on];
speed = [stim1.speed]; speed(speed==-9999) = 0; stim.speed = speed;
modality = [stim1.modality]; modality(modality<0) = -1; stim.modality = modality;

tbeg = prs.tspk(1); tend = mean(stim.tstim) + prs.tspk(2);
time = linspace(tbeg,tend,(tend-tbeg)/prs.dt); stim.time = time;
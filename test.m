dataInMemory = ImportSMR('C:\Users\jlaks\Google Drive\Grad School\Angelaki lab\Thesis project\Path integration\Monkey\Data\test\test2.smr');
data = dataInMemory;
count = 0;
for i=[1 13:14]
    count = count+1;
    speed(count).wf = double(data(i).imp.adc)*data(i).hdr.adc.Scale;
end
fsamp = 5000/6;
dt = 1000*(1/fsamp); nsamp = length(speed(1).wf);
t = linspace(0, dt*nsamp, nsamp);
plot(t,speed(1).wf,'k'); hold on;
plot(t,speed(2).wf,'b');
plot(t,speed(3).wf,'r');
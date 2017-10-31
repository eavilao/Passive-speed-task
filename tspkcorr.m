function [r,p] = tspkcorr(trls,type)

switch type
    case 'true'
        for k=1:length(trls)
            [~,r1(k,:)] = spktimes2spktrain(trls(k).tspk1,0,1.5,0.05); % binwidth = 50ms
            [~,r2(k,:)] = spktimes2spktrain(trls(k).tspk2,0,1.5,0.05); % binwidth = 50ms
        end
    case 'shuffled'
        for k=1:length(trls)
            [~,r1(k,:)] = spktimes2spktrain(trls(k).tspk1,0,1.5,0.05); % binwidth = 50ms
            [~,r2(k,:)] = spktimes2spktrain(trls(k).tspk_rand,0,1.5,0.05); % binwidth = 50ms
        end
end

[ntrls, nt] = size(r1);
r1 = reshape(r1',[1 ntrls*nt]);
r2 = reshape(r2',[1 ntrls*nt]);
[r,p]=corr(r1',r2');
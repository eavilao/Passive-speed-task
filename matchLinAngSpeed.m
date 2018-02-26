% Match neurons that were responsive in linear and angular speed. 


% Only tuned neurons

neuron_match; 

compareInfo.labels = {'monk_id','session_id', 'ch_id', 'unitNum_linear', 'unitNum_angular', 'tuned_ves_linear', 'tuned_ves_angular', 'tuned_vis_linear', 'tuned_vis_angular', 'resp_ves_linear', 'resp_ves_angular', 'resp_vis_linear', 'resp_vis_angular'};

% linear
monk_ids1 = [experiments(1).singleunits.monk_id];
session_ids1 = [experiments(1).singleunits.session_id];
ch_nos1 = [experiments(1).singleunits.channel_no];
% angular
monk_ids2 = [experiments(2).singleunits.monk_id];
session_ids2 = [experiments(2).singleunits.session_id];
ch_nos2 = [experiments(2).singleunits.channel_no];
for i = 1:length(compareInfo.data)
    indx1 = find(monk_ids1 == compareInfo.data(i,1) & session_ids1 == compareInfo.data(i,2) & ch_nos1 == compareInfo.data(i,3)); %Linear speed
    indx1 = indx1(compareInfo.data(i,4)-1);
    indx2 = find(monk_ids2 == compareInfo.data(i,1) & session_ids2 == compareInfo.data(i,2) & ch_nos2 == compareInfo.data(i,3)); %Angular speed
    indx2 = indx2(compareInfo.data(i,5)-1);
    
    
    %%%% Tuning
    % Tuned to vestibular or not
    compareInfo.data(i,6)= experiments(1).singleunits(indx1).ves.stats.flags.tuning & experiments(1).singleunits(indx1).ves.stats.flags.exc;  % ves lin
    compareInfo.data(i,7)= experiments(2).singleunits(indx2).ves.stats.flags.tuning & experiments(2).singleunits(indx2).ves.stats.flags.exc; % ves ang
    % Tuned to visual or not
    compareInfo.data(i,8)= experiments(1).singleunits(indx1).vis.stats.flags.tuning & experiments(1).singleunits(indx1).vis.stats.flags.exc;   % vis lin
    compareInfo.data(i,9)= experiments(2).singleunits(indx2).vis.stats.flags.tuning & experiments(2).singleunits(indx2).vis.stats.flags.exc; % vis ang
    
    %%%%% Responsiveness
     % Responsive to vestibular or not
    compareInfo.data(i,10)= experiments(1).singleunits(indx1).ves.stats.flags.exc | experiments(1).singleunits(indx1).ves.stats.flags.sup;  % ves lin
    compareInfo.data(i,11)= experiments(2).singleunits(indx2).ves.stats.flags.exc | experiments(2).singleunits(indx2).ves.stats.flags.sup; % ves ang
    % Responsive to visual or not
    compareInfo.data(i,12)= experiments(1).singleunits(indx1).vis.stats.flags.exc | experiments(1).singleunits(indx1).vis.stats.flags.sup;   % vis lin
    compareInfo.data(i,13)= experiments(2).singleunits(indx2).vis.stats.flags.exc | experiments(2).singleunits(indx2).vis.stats.flags.sup; % vis ang
    
    
    %Separate 
    
    % Excitatory vestibular
    compareInfo.data(i,14)= experiments(1).singleunits(indx1).ves.stats.flags.exc ;  % ves lin Only exc
    compareInfo.data(i,15)=  experiments(2).singleunits(indx2).ves.stats.flags.exc; % ves ang Only exc
    
    % Suppresive vestibular
    compareInfo.data(i,16)= experiments(1).singleunits(indx1).ves.stats.flags.sup;   % ves lin Only sup
    compareInfo.data(i,17)= experiments(2).singleunits(indx2).ves.stats.flags.sup; % ves ang Only sup
    
    % Excitatory visual
    compareInfo.data(i,18)= experiments(1).singleunits(indx1).vis.stats.flags.exc ;  % vis lin Only exc
    compareInfo.data(i,19)=  experiments(2).singleunits(indx2).vis.stats.flags.exc; % vis ang Only exc
    
    % Suppresive visual
    compareInfo.data(i,20)= experiments(1).singleunits(indx1).vis.stats.flags.sup;   % vis lin Only sup
    compareInfo.data(i,21)= experiments(2).singleunits(indx2).vis.stats.flags.sup; % vis ang Only sup
    
    
    
    
end

% Check how many there are

%%% Tuned
vesMatch = sum(compareInfo.data(:,6) == 1 & compareInfo.data(:,7)==1)
visMatch = sum(compareInfo.data(:,8)==1 & compareInfo.data(:,9)==1)


% Responsive

vesMatchResp = sum(compareInfo.data(:,10) == 1 & compareInfo.data(:,11)==1)
visMatchResp = sum(compareInfo.data(:,12)==1 & compareInfo.data(:,13)==1)


% Responsive separate 

vesMatchExc= sum(compareInfo.data(:,14) == 1 & compareInfo.data(:,15)==1)
vesMatchSup= sum(compareInfo.data(:,16) == 1 & compareInfo.data(:,17)==1)

visMatchExc= sum(compareInfo.data(:,18) == 1 & compareInfo.data(:,19)==1)
visMatchSup= sum(compareInfo.data(:,20) == 1 & compareInfo.data(:,21)==1)



% for sanity check

% Linear speed
for m = 1:length(experiments(1).singleunits);
lin_tunedAll(m)= experiments(1).singleunits(m).vis.stats.flags.tuning; 
lin_excOnly(m) = experiments(1).singleunits(m).vis.stats.flags.exc;
end
lin_tunedAndExc = sum (lin_tunedAll & lin_excOnly);
lin_tunedAll_sum = sum(lin_tunedAll);
lin_excOnly_sum = sum(lin_excOnly);

% Angular speed
for m = 1:length(experiments(2).singleunits);
ang_tunedAll(m)= experiments(2).singleunits(m).vis.stats.flags.tuning; 
ang_excOnly(m) = experiments(2).singleunits(m).vis.stats.flags.exc;
end
ang_tunedAndExc = sum (ang_tunedAll &  ang_excOnly);
ang_tunedAll_sum = sum(ang_tunedAll);
ang_excOnly_sum = sum(ang_excOnly);




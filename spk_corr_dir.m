% code to compute from experiments_HD.mat spk corr to left or right. 

% singleunits
units = experiments.singleunits;
cond = {'ves','vis','com'};

for cell=1:length(units)
    for i=1:length(cond)
        headings= units(cell).(cond{i}).headings;
        choice = units(cell).(cond{i}).choice;
        for trialNum =1:length(headings)
            if choice(trialNum) == 0 & headings(trialNum)>0;
                direction_sua(cell,i,trialNum) = 1;
            elseif choice(trialNum) == 0 & headings(trialNum)<0;
                direction_sua(cell,i,trialNum) = -1;
            else choice(trialNum) == 0 & headings(trialNum)==0;
                direction_sua(cell,i,trialNum) = 0;
            end
        end
    end
end

% plot distribution of spks and heading per cond (requires manual intervention)
for cell=1:length(units)
    unique_head = unique(units(cell).(cond{i}).npsk_choice_heading(:,3));
    for i=1:length(cond)
        indx_head_pos= find(units(cell).(cond{i}).npsk_choice_heading(:,3) == unique_head(1));
        indx_head_neg = find(units(cell).(cond{i}).npsk_choice_heading(:,3) == unique_head(9));
        
        spks_pos(:,i) = units(cell).(cond{i}).npsk_choice_heading(indx_head_pos,1);
        spks_neg(:,i) = units(cell).(cond{i}).npsk_choice_heading(indx_head_neg,1);
    end
    
    % plot
    figure; hold on;
    ybounds = [min(min(spks_pos)) max(max(spks_pos))];  
    xbounds = [min(min(spks_neg)) max(max(spks_neg))]; 
    for i=1:length(cond)
        plot(spks_pos(:,i),spks_neg(:,i),'.', 'Color',[1 2 3]==i, 'MarkerSize', 18);
        xlim(xbounds); ylim(ybounds);
        title('spks vs heading')
    end
end





%% multiunits

units = experiments.multiunits;
cond = {'ves','vis','com'};

for cell=1:length(units)
    for i=1:length(cond)
        headings= units(cell).(cond{i}).headings;
        choice = units(cell).(cond{i}).choice;
        for trialNum =1:length(headings)
            if choice(trialNum) == 0 & headings(trialNum)>0;
                direction_mua(cell,i,trialNum) = 1;
            elseif choice(trialNum) == 0 & headings(trialNum)<0;
                direction_mua(cell,i,trialNum) = -1;
            else choice(trialNum) == 0 & headings(trialNum)==0;
                direction_mua(cell,i,trialNum) = 0;
            end
        end
    end
end
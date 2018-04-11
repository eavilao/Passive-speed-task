% plot fraction of excitatory and suppresive neurons for different windows
% load experiments_win_XXms.mat

%load
% fifty = load('experiments_win_50ms.mat');
% one_fifty = load('experiments_win_150ms.mat'); one_fifty.experiments.addpopulations('all', 'singleunits'); one_fifty.experiments.addpopulations('all', 'multiunits');
% two_fifty = load('experiments_win_250ms.mat'); two_fifty.experiments.addpopulations('all', 'singleunits'); two_fifty.experiments.addpopulations('all', 'multiunits');
% three_fifty = load('experiments_win_350ms.mat'); three_fifty.experiments.addpopulations('all', 'singleunits'); three_fifty.experiments.addpopulations('all', 'multiunits');
% four_fifty = load('experiments_win_450ms.mat'); four_fifty.experiments.addpopulations('all', 'singleunits'); four_fifty.experiments.addpopulations('all', 'multiunits');
% five_fifty = load('experiments_win_550ms.mat'); five_fifty.experiments.addpopulations('all', 'singleunits'); five_fifty.experiments.addpopulations('all', 'multiunits');
% original = load('experiments.mat');

% extract
cond = {'ves','vis', 'com'};
windows = [50 150 250 350 450 550];

for i=1:length(cond)
    % get fraction of excitatory neurons
    exc(i,:) = [fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; one_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; two_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; ...
        three_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; four_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; five_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons];
    
    % get fraction of suppressive neurons
    sup(i,:) = [fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; one_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; two_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; ...
        three_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; four_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; five_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons];
end

%identity
for i=1:length(cond)
    %exc
    for cellNum = 1:length(fifty.experiments(1).populations(1).(cond{i}).exc.indx);
        id_exc(cellNum,:,i) = [fifty.experiments.singleunits((cellNum)).(cond{i}).stats.flags.exc; one_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc; two_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc;...
            three_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc; four_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc; five_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc];
    end
    % sup
    for cellNum = 1:length(fifty.experiments(1).populations(1).(cond{i}).exc.indx);
        id_sup(cellNum,:,i) = [fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; one_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; two_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup;...
            three_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; four_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; five_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup];
    end
end


%% plot
%% fraction of neurons
subplot(2,1,1); hold on;
for i=1:length(cond)
    plot(windows,exc(i,:),'Color',[1 2 3]==i, 'LineWidth', 2)
end
set(gca,'xlim', [0 600], 'ylim', [0 0.5], 'xTick', [windows], 'yTick', [0 0.25 0.5], 'TickDir', 'out','FontSize', 18);
box off
title('Excitatory'); xlabel('window size (ms)'); ylabel('fraction of neurons')

% plot
subplot(2,1,2); hold on;
for i=1:length(cond)
    plot(windows,sup(i,:),'Color',[1 2 3]==i, 'LineWidth', 2)
end
set(gca,'xlim', [0 600], 'ylim', [0 0.1], 'xTick', [windows], 'yTick', [0 0.05 0.1], 'TickDir', 'out','FontSize', 18);
box off
title('Suppresive'); xlabel('window size (ms)'); ylabel('fraction of neurons')

%% neuron identity
map = [1 1 1
    0 0 0];

nunits = length(fifty.experiments(1).singleunits);  
% exc
for i=1:length(id_exc(1,1,:))
    figure; hold on;
    imagesc(windows,1:nunits,id_exc(:,:,i), [0 1]);
    colormap(map); set(gca, 'xtick', [windows]);
    ylabel('neuron'); xlabel('window size (ms)');
    title(cond{i}); box off
end

% sup
for i=1:length(id_sup(1,1,:))
    figure; hold on;
    imagesc(windows,1:nunits,id_sup(:,:,i), [0 1]);
    colormap(map); set(gca, 'xtick', [windows]);
    ylabel('neuron'); xlabel('window size (ms)');
    title(cond{i}); box off
end

%% cosine similarity index
for i=1:length(id_exc(1,1,:))
    for k = 1:length(id_exc(1,:,1))
        if k==6
            break
        else
            cos_sim_indx(k,i) = sum(id_exc(:,k,i).*id_exc(:,k+1,i))/(sqrt(sum(id_exc(:,k,i))*sum(id_exc(:,k+1,i))));
        end
    end
end

% plot 
figure; hold on; 
for i=1:length(cond)
    plot(cos_sim_indx(:,i),'Color',[1 2 3]==i, 'LineWidth', 2)
    set(gca,'xlim', [0 6], 'ylim', [0.3 1], 'xtick',[ 1 2 3 4 5],'ytick',[0.3 1], 'TickDir', 'out', 'FontSize', 18);
    xlabel ('window comparison'); ylabel('Cosine similarity index');
end






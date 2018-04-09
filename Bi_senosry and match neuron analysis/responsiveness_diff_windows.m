% plot fraction of excitatory and suppresive neurons for different windows
% load experiments_win_XXms.mat

%load
% fifty = load('experiments_win_50ms.mat');
% one_fifty = load('experiments_win_150ms.mat'); one_fifty.experiments.addpopulations('all', 'singleunits'); one_fifty.experiments.addpopulations('all', 'multiunits');
% two_fifty = load('experiments_win_250ms.mat'); two_fifty.experiments.addpopulations('all', 'singleunits'); two_fifty.experiments.addpopulations('all', 'multiunits');
% three_fifty = load('experiments_win_350ms.mat'); three_fifty.experiments.addpopulations('all', 'singleunits'); three_fifty.experiments.addpopulations('all', 'multiunits');
% four_fifty = load('experiments_win_450ms.mat'); four_fifty.experiments.addpopulations('all', 'singleunits'); four_fifty.experiments.addpopulations('all', 'multiunits');
% original = load('experiments.mat');

% extract
cond = {'ves','vis', 'com'};
windows = [50 150 250 350 450];

for i=1:length(cond)
    % get fraction of excitatory neurons
    exc(i,:) = [fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; one_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; two_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; ...
        three_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons; four_fifty.experiments(1).populations(1).(cond{i}).exc.f_neurons;];
    
    % get fraction of suppressive neurons
    sup(i,:) = [fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; one_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; two_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; ...
        three_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons; four_fifty.experiments(1).populations(1).(cond{i}).sup.f_neurons;];
end

%identity
for i=1:length(cond)
    %exc
    for cellNum = 1:length(fifty.experiments(1).populations(1).(cond{i}).exc.indx);
        id_exc(cellNum,:,i) = [fifty.experiments.singleunits((cellNum)).(cond{i}).stats.flags.exc; one_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc; two_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc;...
            three_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc; four_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.exc];
    end
    % sup
    for cellNum = 1:length(fifty.experiments(1).populations(1).(cond{i}).exc.indx);
        id_sup(cellNum,:,i) = [fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; one_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; two_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup;...
            three_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup; four_fifty.experiments.singleunits(cellNum).(cond{i}).stats.flags.sup];
    end
end


%% plot
%% fraction of neurons
subplot(2,1,1); hold on;
for i=1:length(cond)
    plot(windows,exc(i,:),'Color',[1 2 3]==i, 'LineWidth', 2)
end
set(gca,'xlim', [0 450], 'ylim', [0 0.5], 'xTick', [windows], 'yTick', [0 0.25 0.5], 'TickDir', 'out');
box off
title('Excitatory'); xlabel('window size (ms)'); ylabel('fraction of neurons')

% plot
subplot(2,1,2); hold on;
for i=1:length(cond)
    plot(windows,sup(i,:),'Color',[1 2 3]==i, 'LineWidth', 2)
end
set(gca,'xlim', [0 450], 'ylim', [0 0.1], 'xTick', [windows], 'yTick', [0 0.05 0.1], 'TickDir', 'out');
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

z=1;
% for i=1:length(id_exc(1,1,:))
% figure; hold on;
% for nunits= 1:length(id_exc(:,1,1))
% plot(windows,id_exc(nunits,:,i), 'ok'); plot(windows,id_exc(nunits,:,i), '.k');
% set(gca, 'xtick', [windows]);
% ylabel('neuron'); xlabel('window size (ms)');
% title(cond{i}); box off
% end
% end


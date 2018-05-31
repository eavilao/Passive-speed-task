function results=compute_psychometric(unit,prs)

cond = {'ves', 'vis', 'com'};
unique_headings = unit(1).ves.stim;

for u = 1:length(unit);
    for i=1:length(cond)
        all_headings{i,:} = unit(u).(cond{i}).headings;
        choice{i,:} = unit(u).(cond{i}).choice;
    end
    
    for i=1:length(cond)
        for j=1:length(unique_headings)
            id{j,:} = all_headings{i} == unique_headings(j);
            nTrials(i,j) = sum(id{j,:});
            nCorrect(i,j)= sum(all_headings{i} == unique_headings(j) & choice{i} == 0);
            pCorrect(i,j) = nCorrect(i,j)/nTrials(i,j);
            if unique_headings(j)<0
                pCorrect(i,j)= 1-pCorrect(i,j);
            end
            results(u).fit_data_psycho_cum{i}(j,1)=unique_headings(j);
            results(u).fit_data_psycho_cum{i}(j,2)=pCorrect(i,j);
            results(u).fit_data_psycho_cum{i}(j,3)=nTrials(i,j);
        end
    end
    
    % fit psychometric function using Wichman's MLE method to estimate threshold and bias(same as TEMPO GUI)
    for i=1:length(cond)
        wichman_psy = pfit(results(u).fit_data_psycho_cum{i},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        results(u).Thresh_psy{i} = wichman_psy.params.est(2);
        results(u).Bias_psy{i} = wichman_psy.params.est(1);
        results(u).psy_perf{i} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
    end
    
    %plot psychometric with fit
    h{1} = 'ro';  f{1} = 'r-';  g{1} = 'ro-';
    h{2} = 'gd';  f{2} = 'g-';  g{2} = 'gd-';
    h{3} = 'bs';  f{3} = 'b-';  g{3} = 'bs-';
    figure; hold on;
    for i=1:length(cond)
        xi = min(unique_headings) : 0.1 : max(unique_headings);
        beta = [0, 1.0];
        %   plot data in logarithmical space instead of linspace
        plot(unique_headings, pCorrect(i,:), h{i}, xi, cum_gaussfit(psy_perf{i}, xi), f{i}, 'MarkerSize', 6, 'Linewidth', 1.5);
        set(gca, 'TickDir', 'out', 'ylim',([0,1]), 'YTick', [0 0.5 1], 'FontSize', 16);
        xlabel('Heading Angles');
        ylabel('Rightward Choices');
        box off
    end
    title('Psychometric function')
    annotation('textbox', [0.55 0.33 0.4 0.1], 'string', ['Threshold (1=Ves, 2=Vis, 3=Com)= ' num2str([Thresh_psy{:}])])
    annotation('textbox',  [0.55 0.2 0.4 0.1], 'string', ['Bias (1=Ves, 2=Vis, 3=Com)= ' num2str([Bias_psy{:}])])
    
    % neurometric function
    %gather
    for i=1:length(cond)
        for j=1:length(unique_headings)
            r{i,j} = [unit(u).(cond{i}).nspk(j)'];
            r_avg(i,:) = [unit(u).(cond{i}).rate_avg.mu];
        end
    end
    
    % decide whether ves and vis is congruent tuning. Fit line by linear
    % regression first and compare the sign of each condition to decide whether
    % congruent or opposite, this is used to check whether congruent cells lead
    % to better neuronal performance in combined condition, and vice versa
    for i=1:length(cond)
        [rr,pp] = corrcoef(unique_headings,r_avg(i,:));
        line_re{i} = rr(1,2);
        line_p{i} = pp(1,2);
    end
    
    % calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
    for i=1:length(cond)
        for j=1:length(unique_headings)-1
            trials_n = all_headings{i} == unique_headings(j);
            fit_data_neuro_cum{i}(j,3) = sum(trials_n);
            if j < (1+length(unique_headings))/2
                results(u).Neuro_correct{i}(j) =  rocN(r{i,length(unique_headings)-j+1}{:},r{i,j}{:},100 );
            else
                results(u).Neuro_correct{i}(j) =  rocN(r{i,length(unique_headings)-j}{:}, r{i,(j+1)}{:},100);
            end
            if line_re{i} > 0
                results(u).Neuro_correct{i}(j) = 1 - Neuro_correct{i}(j);
            end
        end
    end
    
    for i=1:length(cond)
        fit_data_neuro_cum{i}(:,1) = unique_headings(unique_headings~=0);
        fit_data_neuro_cum{i}(:,2) = results(u).Neuro_correct{i}(:);
        wichman_neu = pfit(fit_data_neuro_cum{i}(2:7,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        Thresh_neu{i} = wichman_neu.params.est(2);
        % negative and positive infinite value means flat tuning
        if Thresh_neu{i}<0 | Thresh_neu{i}> 300
            Thresh_neu{i} = 300;
            wichman_neu.params.est(2) = 300;
        end
        results(u).Bias_neu{i} = wichman_neu.params.est(1);
        results(u).neu_perf{i} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
    end
    
    %plot neurometric function
    h{1} = 'ro';  f{1} = 'r-';  g{1} = 'ro-';
    h{2} = 'gd';  f{2} = 'g-';  g{2} = 'gd-';
    h{3} = 'bs';  f{3} = 'b-';  g{3} = 'bs-';
    
    figure; hold on;
    for i=1:length(cond)
        neu_heading = unique_headings(unique_headings~=0);
        xi = min(unique_headings) : 0.1 : max(unique_headings);
        plot(neu_heading, Neuro_correct{i}, h{i},xi,cum_gaussfit(neu_perf{i},xi),f{i},'MarkerSize', 6, 'Linewidth', 1.5);
        set(gca, 'TickDir', 'out', 'ylim',([0,1]), 'xlim',([-30,30]), 'YTick', [0 0.5 1], 'FontSize', 16);
        xlabel('Heading Angles');
        box off
        xlabel('Heading Angles');
        ylabel('Rightward Choices');
        title ('Neurometric function')
    end
end

end
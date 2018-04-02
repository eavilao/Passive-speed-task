function dist = computeDist(experiments,exp_name, unit_type)
% compute distance correlation between pairs of channels

switch exp_name
    case 'linearspeed'
        cond = {'ves','vis','com'};
        exp=1;
    case 'angularspeed'
        cond = {'ves','vis'};
        exp=2;
end

if unit_type == 'singleunits'
    unit_type = 1;
else unit_type == 'multiunits'
    unit_type =2;
end

numOfSessions = length(experiments(exp).sessions);

for condition = 1:length(cond)
    for sessNum = 1:length(experiments(exp).sessions)
        
        % load
        %         t_on = experiments(exp).populations(3).(cond{condition}).all.peak.t_on;
        t_on = experiments(exp).populations(3).(cond{condition}).all.SEP.v_peak(:,2);
        t_on = reshape(t_on,[16 numOfSessions]);  % number of channels (16ch UProbe)/ number of sessions
        
        % distance 1
        distance=1;
        reference=0;
        for i=1:length(t_on(:,sessNum))-1
            reference=reference+1; distance=distance+1;
            d_1{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_1(i,:,sessNum)= d_1{sessNum,reference}';
        end
        
        % distance 2
        distance=2;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_2{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_2(i,:,sessNum)= d_2{sessNum,reference}';
        end
        
        
        % distance 3
        distance=3;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_3{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_3(i,:,sessNum)= d_3{sessNum,reference}';
        end
        
        % distance 4
        distance=4;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_4{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_4(i,:,sessNum)= d_4{sessNum,reference}';
        end
        
        % distance 5
        distance=5;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_5{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_5(i,:,sessNum)= d_5{sessNum,reference}';
        end
        
        % distance 6
        distance=6;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_6{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_6(i,:,sessNum)= d_6{sessNum,reference}';
        end
        
        % distance 7
        distance=7;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_7{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_7(i,:,sessNum)= d_7{sessNum,reference}';
        end
        
        % distance 8
        distance=8;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_8{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_8(i,:,sessNum)= d_8{sessNum,reference}';
        end
        
        % distance 9
        distance=9;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_9{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_9(i,:,sessNum)= d_9{sessNum,reference}';
        end
        
        % distance 10
        distance=10;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_10{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_10(i,:,sessNum)= d_10{sessNum,reference}';
        end
        
        % distance 11
        distance=11;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_11{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_11(i,:,sessNum)= d_11{sessNum,reference}';
        end
        
        % distance 12
        distance=12;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_12{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_12(i,:,sessNum)= d_12{sessNum,reference}';
        end
        
        % distance 13
        distance=13;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_13{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_13(i,:,sessNum)= d_13{sessNum,reference}';
        end
        
        % distance 14
        distance=14;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_14{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_14(i,:,sessNum)= d_14{sessNum,reference}';
        end
        
        % distance 15
        distance=15;
        reference=0;
        for i=1:length(t_on(:,sessNum))-distance
            reference=reference+1; distance=distance+1;
            d_15{sessNum,reference}= t_on([reference distance],sessNum);
            dist.(cond{condition}).dist_15(i,:,sessNum)= d_15{sessNum,reference}';
        end
        
    end
end

%% compute correlation for every distance

for condition = 1:length(cond)
    
        %concatenate
        dist.(cond{condition}).distance_1 = [dist.(cond{condition}).dist_1(:,:,1)];
        dist.(cond{condition}).distance_2 = [dist.(cond{condition}).dist_2(:,:,1)];
        dist.(cond{condition}).distance_3 = [dist.(cond{condition}).dist_3(:,:,1)];
        dist.(cond{condition}).distance_4 = [dist.(cond{condition}).dist_4(:,:,1)];
        dist.(cond{condition}).distance_5 = [dist.(cond{condition}).dist_5(:,:,1)];
        dist.(cond{condition}).distance_6 = [dist.(cond{condition}).dist_6(:,:,1)];
        dist.(cond{condition}).distance_7 = [dist.(cond{condition}).dist_7(:,:,1)];
        dist.(cond{condition}).distance_8 = [dist.(cond{condition}).dist_8(:,:,1)];
        dist.(cond{condition}).distance_9 = [dist.(cond{condition}).dist_9(:,:,1)];
        dist.(cond{condition}).distance_10 = [dist.(cond{condition}).dist_10(:,:,1)];
        dist.(cond{condition}).distance_11 = [dist.(cond{condition}).dist_11(:,:,1)];
        dist.(cond{condition}).distance_12 = [dist.(cond{condition}).dist_12(:,:,1)];
        dist.(cond{condition}).distance_13= [dist.(cond{condition}).dist_13(:,:,1)];
        dist.(cond{condition}).distance_14 = [dist.(cond{condition}).dist_14(:,:,1)];
        dist.(cond{condition}).distance_15 = [dist.(cond{condition}).dist_15(:,:,1)];
        
        for sessNum = 2:length(experiments(exp).sessions)
        
        dist.(cond{condition}).distance_1 = [dist.(cond{condition}).distance_1; dist.(cond{condition}).dist_1(:,:,sessNum)];
        dist.(cond{condition}).distance_2 = [dist.(cond{condition}).distance_2; dist.(cond{condition}).dist_2(:,:,sessNum)];
        dist.(cond{condition}).distance_3 = [dist.(cond{condition}).distance_3; dist.(cond{condition}).dist_3(:,:,sessNum)];
        dist.(cond{condition}).distance_4 = [dist.(cond{condition}).distance_4; dist.(cond{condition}).dist_4(:,:,sessNum)];
        dist.(cond{condition}).distance_5 = [dist.(cond{condition}).distance_5; dist.(cond{condition}).dist_5(:,:,sessNum)];
        dist.(cond{condition}).distance_6 = [dist.(cond{condition}).distance_6; dist.(cond{condition}).dist_6(:,:,sessNum)];
        dist.(cond{condition}).distance_7 = [dist.(cond{condition}).distance_7; dist.(cond{condition}).dist_7(:,:,sessNum)];
        dist.(cond{condition}).distance_8 = [dist.(cond{condition}).distance_8; dist.(cond{condition}).dist_8(:,:,sessNum)];
        dist.(cond{condition}).distance_9 = [dist.(cond{condition}).distance_9; dist.(cond{condition}).dist_9(:,:,sessNum)];
        dist.(cond{condition}).distance_10 = [dist.(cond{condition}).distance_10; dist.(cond{condition}).dist_10(:,:,sessNum)];
        dist.(cond{condition}).distance_11 = [dist.(cond{condition}).distance_11; dist.(cond{condition}).dist_11(:,:,sessNum)];
        dist.(cond{condition}).distance_12 = [dist.(cond{condition}).distance_12; dist.(cond{condition}).dist_12(:,:,sessNum)];
        dist.(cond{condition}).distance_13 = [dist.(cond{condition}).distance_13; dist.(cond{condition}).dist_13(:,:,sessNum)];
        dist.(cond{condition}).distance_14 = [dist.(cond{condition}).distance_14; dist.(cond{condition}).dist_14(:,:,sessNum)];
        dist.(cond{condition}).distance_15 = [dist.(cond{condition}).distance_15; dist.(cond{condition}).dist_15(:,:,sessNum)];
        end
        
        % corr & std
        dist.(cond{condition}).corr.dist_1 = nancorr(dist.(cond{condition}).distance_1(:,1), dist.(cond{condition}).distance_1(:,2));
        dist.(cond{condition}).corr.dist_1_se = sqrt((1-dist.(cond{condition}).corr.dist_1)^2/length(dist.ves.distance_1));
        
        dist.(cond{condition}).corr.dist_2 = nancorr(dist.(cond{condition}).distance_2(:,1), dist.(cond{condition}).distance_2(:,2));
        dist.(cond{condition}).corr.dist_2_se = sqrt((1-dist.(cond{condition}).corr.dist_2)^2/length(dist.ves.distance_2));
        
        dist.(cond{condition}).corr.dist_3 = nancorr(dist.(cond{condition}).distance_3(:,1), dist.(cond{condition}).distance_3(:,2));
        dist.(cond{condition}).corr.dist_3_se = sqrt((1-dist.(cond{condition}).corr.dist_3)^2/length(dist.ves.distance_3));
        
        dist.(cond{condition}).corr.dist_4 = nancorr(dist.(cond{condition}).distance_4(:,1), dist.(cond{condition}).distance_4(:,2));
        dist.(cond{condition}).corr.dist_4_se = sqrt((1-dist.(cond{condition}).corr.dist_4)^2/length(dist.ves.distance_4));
        
        dist.(cond{condition}).corr.dist_5 = nancorr(dist.(cond{condition}).distance_5(:,1), dist.(cond{condition}).distance_5(:,2));
        dist.(cond{condition}).corr.dist_5_se = sqrt((1-dist.(cond{condition}).corr.dist_5)^2/length(dist.ves.distance_5));
        
        dist.(cond{condition}).corr.dist_6 = nancorr(dist.(cond{condition}).distance_6(:,1), dist.(cond{condition}).distance_6(:,2));
        dist.(cond{condition}).corr.dist_6_se = sqrt((1-dist.(cond{condition}).corr.dist_6)^2/length(dist.ves.distance_6));
        
        dist.(cond{condition}).corr.dist_7 = nancorr(dist.(cond{condition}).distance_7(:,1), dist.(cond{condition}).distance_7(:,2));
        dist.(cond{condition}).corr.dist_7_se = sqrt((1-dist.(cond{condition}).corr.dist_7)^2/length(dist.ves.distance_7));
        
        dist.(cond{condition}).corr.dist_8 = nancorr(dist.(cond{condition}).distance_8(:,1), dist.(cond{condition}).distance_8(:,2));
        dist.(cond{condition}).corr.dist_8_se = sqrt((1-dist.(cond{condition}).corr.dist_8)^2/length(dist.ves.distance_8));
        
        dist.(cond{condition}).corr.dist_9 = nancorr(dist.(cond{condition}).distance_9(:,1), dist.(cond{condition}).distance_9(:,2));
        dist.(cond{condition}).corr.dist_9_se = sqrt((1-dist.(cond{condition}).corr.dist_9)^2/length(dist.ves.distance_9));
        
        dist.(cond{condition}).corr.dist_10 = nancorr(dist.(cond{condition}).distance_10(:,1), dist.(cond{condition}).distance_10(:,2));
        dist.(cond{condition}).corr.dist_10_se = sqrt((1-dist.(cond{condition}).corr.dist_10)^2/length(dist.ves.distance_10));
        
        dist.(cond{condition}).corr.dist_11 = nancorr(dist.(cond{condition}).distance_11(:,1), dist.(cond{condition}).distance_11(:,2));
        dist.(cond{condition}).corr.dist_11_se = sqrt((1-dist.(cond{condition}).corr.dist_11)^2/length(dist.ves.distance_11));
        
        dist.(cond{condition}).corr.dist_12 = nancorr(dist.(cond{condition}).distance_12(:,1), dist.(cond{condition}).distance_12(:,2));
        dist.(cond{condition}).corr.dist_12_se = sqrt((1-dist.(cond{condition}).corr.dist_12)^2/length(dist.ves.distance_12));
        
        dist.(cond{condition}).corr.dist_13 = nancorr(dist.(cond{condition}).distance_13(:,1), dist.(cond{condition}).distance_13(:,2));
        dist.(cond{condition}).corr.dist_13_se = sqrt((1-dist.(cond{condition}).corr.dist_13)^2/length(dist.ves.distance_13));
        
        dist.(cond{condition}).corr.dist_14 = nancorr(dist.(cond{condition}).distance_14(:,1), dist.(cond{condition}).distance_14(:,2));
        dist.(cond{condition}).corr.dist_14_se = sqrt((1-dist.(cond{condition}).corr.dist_14)^2/length(dist.ves.distance_14));
        
        dist.(cond{condition}).corr.dist_15 = nancorr(dist.(cond{condition}).distance_15(:,1), dist.(cond{condition}).distance_15(:,2));
        dist.(cond{condition}).corr.dist_15_se = sqrt((1-dist.(cond{condition}).corr.dist_15)^2/length(dist.ves.distance_15));
    
    
    dist.(cond{condition}).corr.all = [dist.(cond{condition}).corr.dist_1;dist.(cond{condition}).corr.dist_2;dist.(cond{condition}).corr.dist_3;dist.(cond{condition}).corr.dist_4;...
        dist.(cond{condition}).corr.dist_5;dist.(cond{condition}).corr.dist_6;dist.(cond{condition}).corr.dist_7;dist.(cond{condition}).corr.dist_8;dist.(cond{condition}).corr.dist_9;...
        dist.(cond{condition}).corr.dist_10;dist.(cond{condition}).corr.dist_11;dist.(cond{condition}).corr.dist_12;dist.(cond{condition}).corr.dist_13;dist.(cond{condition}).corr.dist_14;dist.(cond{condition}).corr.dist_15];
    
    dist.(cond{condition}).corr.se_all = [dist.(cond{condition}).corr.dist_1_se;dist.(cond{condition}).corr.dist_2_se;dist.(cond{condition}).corr.dist_3_se;dist.(cond{condition}).corr.dist_4_se;...
        dist.(cond{condition}).corr.dist_5_se;dist.(cond{condition}).corr.dist_6_se;dist.(cond{condition}).corr.dist_7_se;dist.(cond{condition}).corr.dist_8_se;dist.(cond{condition}).corr.dist_9_se;...
        dist.(cond{condition}).corr.dist_10_se;dist.(cond{condition}).corr.dist_11_se;dist.(cond{condition}).corr.dist_12_se;dist.(cond{condition}).corr.dist_13_se;dist.(cond{condition}).corr.dist_14_se;dist.(cond{condition}).corr.dist_15_se];
    
    
end

dist = dist;

%% plot

% scatter ves
% plot(dist.ves.distance_1(:,1),dist.ves.distance_1(:,2), '.r'); hold on;
% plot(dist.ves.distance_2(:,1),dist.ves.distance_2(:,2), '.g');
interDist =[100:100:1500]; 
shadedErrorBar(interDist,dist.ves.corr.all,dist.ves.corr.se_all,'r'); hold on;
shadedErrorBar(interDist,dist.vis.corr.all,dist.vis.corr.se_all,'g');
shadedErrorBar(interDist,dist.com.corr.all,dist.com.corr.se_all,'b');
box off
set(gca, 'TickDir', 'out', 'ylim',([0.8 1]), 'ytick', [0.8 0.9 1], 'FontSize',18)
ylabel('corr coeff')
xlabel('Distance (mm)')

end

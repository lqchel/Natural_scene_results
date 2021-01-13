data = b2_g3_raw(2:end,:);
data = table2array(data(~isnan(data.response),:));
[unique_rows,index1,index2]= unique(data(:,1:9),'rows');
pilot.b2_g3 = data(index1,:);
colours = cbrewer('qual', 'Set1', 8); 



%% calculate and plot real & simulated catch trial accuracy
pilot_results = [pilot.b2_g1;pilot.b2_g2;pilot.b2_g3];
subject_id = unique(pilot_results(:,1));
% real catch accuracy
catch_accuracy = table(zeros(length(subject_id),1),zeros(length(subject_id),1),'VariableNames',{'Subject','Accuracy'});
for sub = 1:length(subject_id)
    current_subject = pilot_results(pilot_results(:,1) == subject_id(sub) & pilot_results(:,11) == 0 & pilot_results(:,2) >4,:);
    catch_accuracy.Accuracy(sub,:) = sum(current_subject(:,5) == current_subject(:,8)& current_subject(:,6) == current_subject(:,9))/...
        size(current_subject,1);
    catch_accuracy.Subject(sub,:) = subject_id(sub);
end

% catch accuracy simulation
nSub = 1000; 
nCatch = 12; 
nOption = 8 ; 

for iSub = 1:nSub 
for iTrial = 1: nCatch
if rand < 1/nOption 
corr(iTrial, iSub) = 1;
else 
corr(iTrial, iSub) = 0;
end
end
end

pCorr = sum(corr,1)/12;

% calculate cumulative percentages
[Y_s,edges_s] = histcounts(pCorr,1000); 
Ys_cumulative = cumsum(Y_s)/1000;

[Y,edges] = histcounts(catch_accuracy.Accuracy,22);
Y_cumulative  = cumsum(Y)/22;

% plot cumulative histogram
subplot(1,2,2),plot(edges,[Y_cumulative 1],'LineWidth',1.2);
xlabel('Catch trial accuracy'), ylabel('Cumulative percentage');
xlim([0 max(catch_accuracy.Accuracy)]), ylim([0 1]);
hold on
plot(edges_s,[Ys_cumulative 1],'LineWidth',1.2);
hold off

title('Simulated accuracy [0.01, 0.99] = [0, 0.38]');
set(gca,'FontName','Arial','FontSize',12);
legend({'Real','Simulated'});
legend('boxoff');
%% reaction time

% Exp 2
pilot_results(:,9) = pilot_results(:,8).*pilot_results(:,9);
exp_trials = pilot_results(pilot_results(:,11) ~= 0 & pilot_results(:,2) >4,:);
catch_trials = pilot_results(pilot_results(:,11) == 0 & pilot_results(:,2) >4,:);

for sub = 1:22
for i = 1:9
    exp_rt_matrix(sub,i) = mean(log(exp_trials(exp_trials(:,1) == subject_id(sub) & exp_trials(:,9) == i-5,10)));
    catch_rt_matrix(sub,i) = mean(log(catch_trials(catch_trials(:,1) == subject_id(sub) & catch_trials(:,9) == i-5,10)));
end
end

exp_rt_matrix = [exp_rt_matrix(:,1:4) exp_rt_matrix(:,6:9)];
catch_rt_matrix = [catch_rt_matrix(:,1:4) catch_rt_matrix(:,6:9)];



% Exp 1

exp1_data = importdata('Pooled Results 2.mat');
exp1_data(:,9) = exp1_data(:,8) .* exp1_data(:,9);
exp1_data(:,10) = exp1_data(:,10) .* 100;

for sub = 1:15
for i = 1:9
    exp1_rt_matrix(sub,i) = mean(log(exp1_data(exp1_data(:,1) == sub & exp1_data(:,9) == i-5,10)));
end
end
exp1_rt_matrix = [exp1_rt_matrix(:,1:4) exp1_rt_matrix(:,6:9)];

subplot(1,2,1)
errorbar([1:8],nanmean(exp_rt_matrix,1),nanstd(exp_rt_matrix,1)/sqrt(22),'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
hold on
errorbar([1:8],nanmean(catch_rt_matrix,1),nanstd(catch_rt_matrix,1)/sqrt(22),'.-','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta',...
    'MarkerFaceColor','magenta','LineWidth',1);

errorbar([1:8],nanmean(exp1_rt_matrix,1),nanstd(exp1_rt_matrix,1)/sqrt(22),'.-','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
xlabel('Decision x Confidence'), ylabel('Log-transformed RT');
xlim([1 8]);
hold off

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
legend({'Exp 2','Exp 2 Catch','Exp 1'},'Box','Off');
xlabel('Decision x Confidence'), ylabel('Log-transformed RT');
xlim([1 8]);


% comparing Exp 1 and 2 reaction time

for i = 1:8
    temp1 = [exp1_rt_matrix(:,i) i + zeros(15,1) ones(15,1)];
    temp2 = [exp_rt_matrix(:,i) i+zeros(22,1) 2+zeros(22,1)];
    temp_c = [temp1;temp2];
    if i == 1
        rt_table = temp_c;
    else
        rt_table = [rt_table;temp_c];
    end
    clear temp1
    clear temp2
    clear temp_c
end

rt_table = table(rt_table(:,1),rt_table(:,2),rt_table(:,3),'VariableNames',{'RT','DxC','Exp'});
lm_rt = fitlme(rt_table,'RT~DxC*Exp');
lm_DxC = fitlme(rt_table,'RT ~ DxC:Exp + Exp');
lm_Exp = fitlme(rt_table,'RT ~ DxC + DxC:Exp');
lm_interaction = fitlme(rt_table,'RT ~ DxC + Exp');

DxC_effect = compare(lm_DxC,lm_rt)
Exp_effect = compare(lm_Exp,lm_rt)
interaction_effect = compare(lm_interaction, lm_rt)



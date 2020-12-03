data = b2_g3_raw(2:end,:);
data = table2array(data(~isnan(data.response),:));
[unique_rows,index1,index2]= unique(data(:,1:9),'rows');
pilot.b2_g3 = data(index1,:);
colours = cbrewer('qual', 'Set1', 8); 
subject_id = unique(pilot_results(:,1));


%% calculate catch trial accuracy for each participant
pilot_results = [pilot.b2_g1;pilot.b2_g2;pilot.b2_g3];
catch_accuracy = table(zeros(length(subject_id),1),zeros(length(subject_id),1),'VariableNames',{'Subject','Accuracy'});
for sub = 1:length(subject_id)
    current_subject = pilot_results(pilot_results(:,1) == subject_id(sub) & pilot_results(:,11) == 0 & pilot_results(:,2) >4,:);
    catch_accuracy.Accuracy(sub,:) = sum(current_subject(:,5) == current_subject(:,8)& current_subject(:,6) == current_subject(:,9))/...
        size(current_subject,1);
    catch_accuracy.Subject(sub,:) = subject_id(sub);
end

[Y,edges] = histcounts(catch_accuracy.Accuracy,22);
Y_cumulative  = cumsum(Y);
subplot(1,2,2),plot(edges,[Y_cumulative 22],'LineWidth',1.2);
xlabel('Catch trial accuracy'), ylabel('Cumulative counts');
set(gca,'FontName','Arial','FontSize',12);
xlim([min(catch_accuracy.Accuracy) max(catch_accuracy.Accuracy)]), ylim([0 23]);


%% simulation of catch accuracy
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
histogram(pCorr);
xlabel('Accuracy'),ylabel('Count'),title('Quantile [0.01, 0.99] = [0, 0.38]');
set(gca,'FontName','Arial','FontSize',12);
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

subplot(1,2,1)
errorbar([1:8],nanmean(exp_rt_matrix,1),nanstd(exp_rt_matrix,1)/sqrt(22),'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
hold on
errorbar([1:8],nanmean(catch_rt_matrix,1),nanstd(catch_rt_matrix,1)/sqrt(22),'.-','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta',...
    'MarkerFaceColor','magenta','LineWidth',1);
hold off
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
legend({'Experiment','Catch'},'Box','Off');
xlabel('Decision x Confidence'), ylabel('Log-transformed RT');
xlim([1 8]);

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
errorbar([1:8],nanmean(exp1_rt_matrix,1),nanstd(exp1_rt_matrix,1)/sqrt(22),'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
xlabel('Decision x Confidence'), ylabel('Log-transformed RT');
xlim([1 8]);


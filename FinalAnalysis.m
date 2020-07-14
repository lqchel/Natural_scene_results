clear all
Results = importdata('Pooled Results.mat');

Results(:,9) = Results(:,8).*Results(:,9);
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;

% type 1 ROC curve analysis
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);

grandmatrix = zeros(2,15,2);

for condition = 1:2
    Results_NC = Results_N(Results_N(:,4)==condition - 1,:);
    Results_APC = Results_AP(Results_AP(:,4)==condition -1,:);
    
    
for sub = 1:15 
    Confidence_N = Results_NC(Results_NC(:,1)==sub,9);
    Confidence_AP = Results_APC(Results_APC(:,1)==sub,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    grandmatrix(1,sub,condition) = AUC;
    
    Confidence_N = [];
    Confidence_AP = [];
    Cumulative_NCounts = [];
    Cumulative_APCounts = [];
    Cumulative_Hit = [];
    Cumulative_FA = [];
    AUC = [];
end
end
%% hypo 2
grandmatrix = zeros(15,2);
for condition = 1:2
    
    for sub = 1:15
    if condition ==1
        Confidence_P = Results(Find_Congruent_CP & Results(:,1)==sub,9);
        Confidence_A = Results(Find_Congruent_IP & Results(:,1)==sub,9);
    else
        Confidence_P = Results(Find_Incongruent_IP & Results(:,1)==sub,9);
        Confidence_A = Results(Find_Incongruent_CP & Results(:,1)==sub,9);
    end

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_P == -i);
    Confidence_NCounts(i+5) = sum(Confidence_A == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_A);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_P);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_A);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_P);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

grandmatrix(sub,condition)= AUC;

clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end

% data distribution, plotting within-subject errors
% boxplot(grandmatrix),xticklabels({'Congruent condition','Incongruent condition'});
% ylabel('AUC');
% set(gca,'FontSize',12);

grandmean = mean(mean(grandmatrix));

for sub = 1:15
    indv_mean = mean(grandmatrix(sub,:));
    grandmatrix(sub,:) = grandmatrix(sub,:)- indv_mean + grandmean;
end
se1 = std(grandmatrix(:,1))/sqrt(15);
se2 = std(grandmatrix(:,2))/sqrt(15);


errorbar([1 2],[mean(grandmatrix(:,1)) mean(grandmatrix(:,2))],[se1 se2],'d-', 'MarkerSize',2, 'MarkerEdgeColor','red','MarkerFaceColor','red',...
    'LineWidth',1.2);
hold on
plot([0.5 2.5],[0.5 0.5],'k--');
hold off
set(gca,'FontSize',12);
xlim([0.5 2.5]),xticks([1 2]),ylim([0.3 0.8]), xticklabels({'Congruent condition','Incongruent condition'});



%% t tests
[h,p,ci,stats] = ttest(grandmatrix(1,:,1),0.5,'Alpha',0.05);
[incongruentAP,p,ci,stats] = ttest(grandmatrix(1,:,2),0.5,'Alpha',0.05);
[h,p,ci,stats] = ttest(grandmatrix(1,:,1),grandmatrix(1,:,2),'Alpha',0.05)


[h,p,ci,stats] = ttest(grandmatrix(:,1),0.5,'Alpha',0.025)
[h,p,ci,stats] = ttest(grandmatrix(:,2),0.5,'Alpha',0.025)
[h,p,ci,stats] = ttest(grandmatrix(:,1),matrix1,'Alpha',0.025)
[h,p,ci,stats] = ttest(grandmatrix(:,2),matrix1,'Alpha',0.025)
[h,p,ci,stats] = ttest(grandmatrix(:,1),grandmatrix(:,2))


%% for hypo one

%before getting into this stage, combine confidence (1-4) with response (yes =1, no= -1)
Results_NC = Results(Find_N,:); % all trials with absent test probes
Results_APC = Results(Find_IAP|Find_CAP,:); % all trials with present test probes

matrix1 = zeros(15,1);
for sub = 1:15 
    Results_NC = Results(Find_N,:);
    Confidence_N = Results_NC(Results_NC(:,1)==sub,9);
    Confidence_AP = Results_APC(Results_APC(:,1)==sub,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix1(sub,:)= AUC; % create matrix for individual AUCs



end

%plot boxplot of the AUCs
boxplot(matrix1)
set(gca,'FontSize',12,'XTickLabel','Present v.s. absent patches','FontSize',12);
ylabel('AUC'),ylim([0.5 1]);

%t-test of AUCs against chance (0.5)
se1 = std(matrix1(:,1))/sqrt(15);
[h,p,ci]= ttest(matrix1(:,1),0.5)


 fig_position = [200 200 600 400];
 f7 = figure('Position', fig_position);
h1 = raincloud_plot(matrix1, 'box_on', 1, 'color', colours(5,:), 'alpha', 0.6, 'bandwidth', .2,'density_type', 'ks',...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);

% sample analysis for individual

    Results_APC = Results(Find_CAP|Find_IAP,:);
    Results_NC = Results(Find_N,:);
    Confidence_N = Results_NC(Results_NC(:,1)==7,9);
    Confidence_AP = Results_APC(Results_APC(:,1)==7,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end
    colours = cbrewer('qual','Set3',8);
    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

    h = area(Cumulative_FA, Cumulative_Hit);
    h.FaceColor = colours(3,:);
    h.EdgeColor = colours(3,:);
    hold on
    plot(Cumulative_FA, Cumulative_Hit,'bo-','LineWidth',1);
    plot([0 1], [0 1],'k-');
    hold off
    xlabel('FA rate'),ylabel('Hit rate');
    set(gca,'FontSize',12)
    axis square;

for i = 1:9
    frequency_N(i) = size(Results(Results(:,1)==7& Find_N & Results(:,9)== i-5,:),1);
end

frequency_N = [frequency_N(1:4) frequency_N(6:9)];
%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});

%response value distribution for AP

for i = 1:9
    frequency_AP(i) = size(Results(Results(:,1)==7&(Find_IAP|Find_CAP) & Results(:,9)== i-5,:),1);
end
%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});

frequency_AP = [frequency_AP(1:4) frequency_AP(6:9)];
frequency_All = [reshape(frequency_AP,[8,1]) reshape(frequency_N,[8,1])];

b = bar(frequency_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
b(1).FaceColor = colours(6,:);
b(2).FaceColor = colours(5,:);
xlabel('Decision x confidence'),
ylabel('Frequency');
legend({'Patch present','Patch absent'});

%% for hypothesis 3
% calculate AUC on each location level
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;

Results = [Results eccentricity];
location = [0 6.48 9.16];


matrix3 = zeros(15, 3);
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_N,:);
    indvP = Results(Results(:,1)==sub & (Find_CAP|Find_IAP),:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i); % AUC calculation
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix3(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

% sort AUC, location and individual ID into a table, then fit lme model
AUC = reshape(matrix3,[45,1]);
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[3,1]);
location = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];


% data = table(AUC, subjects, location,'VariableName',{'AUC','subject','location'});
% lm1 = fitlme(data, 'AUC ~ location + (location|subject) + (1|subject)');
% lm2 = fitlme(data, 'AUC ~ 1 +(location|subject) + (1|subject) ');
% 
% compare(lm2,lm1);

% visualising the effect by plotting within-subjects errorbar

population_m = mean(mean(mean(matrix3))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix3,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix3_corr = matrix3 - repmat(indv_dff, [1,3]);

%calculate within-subject sem on each location level

for i = 1:3
    se(i)= std(matrix3_corr(:,i))/sqrt(15);
end

subplot(1,2,2),errorbar([0 6.48 9.16], mean(matrix3,1),se,'d-', 'MarkerSize',2,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
xlim([0 10]), xticks([0 6.48 9.16]),ylim([0.4 1]);
set(gca,'FontSize',12);
ylabel('AUC'),xlabel('Eccentricity (in dva)');
title('b)','FontWeight','normal');
hold on
plot([0 10], [0.5 0.5],'k--');
hold off

% post-hoc t tests
[h,p,ci] = ttest(matrix3(:,1),matrix3(:,2),0.025)
[h,p,ci] = ttest(matrix3(:,2),matrix3(:,3),0.025)
[h,p,ci] = ttest(matrix3(:,2)-matrix3(:,1), matrix3(:,3)-matrix3(:,2))
lm3 = step(lm1,'quadratic')


%% for hypothesis 4

%congruent
clear indvN
clear indvP
clear indvN_loc
clear indvP_loc

matrix4 = zeros(15, 3);
location = [0 6.48 9.16];
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Congruent_IP,:);
    indvP = Results(Results(:,1)==sub & Find_Congruent_CP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i); % AUC calculation
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix4(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

% sort AUC, location and individual ID into a table, then fit lme model
% AUC = reshape(matrix3,[45,1]);
% subjects(1:15,1)= 1:1:15;
% subjects = repmat(subjects,[3,1]);
% location = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];
% 
% 
% data = table(AUC, subjects, location,'VariableName',{'AUC','subject','location'});
% lm1 = fitlme(data, 'AUC ~ location + (location|subject) + (1|subject)');
% lm2 = fitlme(data, 'AUC ~ 1 +(location|subject) + (1|subject) ');
% 
% compare(lm2,lm1);


% visualising the effect by plotting within-subjects errorbar

population_m = mean(mean(mean(matrix4))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix4,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix4_corr = matrix4 - repmat(indv_dff, [1,3]);

%calculate within-subject sem on each location level

for i = 1:3
    se(i)= std(matrix4_corr(:,i))/sqrt(15);
end


%incongruent condition

clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff

matrix5 = zeros(15, 3);

for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Incongruent_CP,:);
    indvP = Results(Results(:,1)==sub & Find_Incongruent_IP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i); % AUC calculation
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix5(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

% plotting within-subjects errorbar

population_m = mean(mean(mean(matrix5))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix5,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix5_corr = matrix5 - repmat(indv_dff, [1,3]);

%calculate within-subject sem on each location level

for i = 1:3
    se1(i)= std(matrix5_corr(:,i))/sqrt(15);
end

% visualisation for congruent and incongruent condition


errorbar([0 6.48 9.16], mean(matrix4,1),se,'d-', 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
xlim([0 10]), xticks([0 6.48 9.16]),ylim([0.4 0.8]);
set(gca,'FontSize',12);
ylabel('AUC'),xlabel('Eccentricity (in dva)');
hold on
errorbar([0 6.48 9.16], mean(matrix5,1),se1,'d-', 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
plot([0 10],[0.5 0.5],'k--');
hold off
legend({'Congruent condition','Incongruent condition'})

AUC = [reshape(matrix4,[45,1]); reshape(matrix5,[45,1])];
clear subjects
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[6,1]);
eccentricity = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];
eccentricity = [eccentricity; eccentricity];
condition = [zeros(45,1);ones(45,1)];

data = table(AUC, subjects, eccentricity,condition,'VariableName',{'AUC','subject','eccentricity','condition'});
lm4 = fitlme(data, 'AUC ~ eccentricity*condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml1 = fitlme(data, 'AUC ~ eccentricity:condition + condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml2 = fitlme(data, 'AUC ~ eccentricity + condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml3 = fitlme(data, 'AUC ~ eccentricity:condition + eccentricity + (eccentricity|subject) + (condition|subject)+(1|subject)');

c_ecc = compare(lml1,lm4)
c_interaction = compare(lml2,lm4)
c_condition = compare(lml3,lm4)

%post-hoc tests
data2 = table(AUC(1:45,:),subjects(1:45,:),eccentricity(1:45,:), condition(1:45,:),'VariableName',{'AUC','subject','eccentricity','condition'});
lm5 = fitlme(data2, 'AUC ~ eccentricity + (eccentricity|subject)+(1|subject)')
lml4 = fitlme(data2, 'AUC ~ 1 + (eccentricity|subject)+(1|subject)')

data3 = table(AUC(46:90,:),subjects(46:90,:),eccentricity(46:90,:), condition(46:90,:),'VariableName',{'AUC','subject','eccentricity','condition'});
lm6 = fitlme(data3, 'AUC ~ eccentricity + (eccentricity|subject) + (1|subject)')
lml5 = fitlme(data3, 'AUC ~ 1 + (eccentricity|subject) + (1|subject)')

compare(lml4,lm5)
compare(lml5,lm6)

compare(lm2,lm1);

%% additional hit and cr analysis

% hypothesis 1


colours = cbrewer('qual', 'Pastel1', 8); 
grandmatrix = zeros(2,2);
grandmean1 = sum(Results(Find_CAP|Find_IAP|Find_N,8)==1)/length(Results(Find_CAP|Find_IAP|Find_N,:));
matrix1 = zeros(15,2);

for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
 
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
    end
    matrix1(sub,condition) = accuracy;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end
[h,p,ci]= ttest(matrix1(:,1),matrix1(:,2));

% calculate standard errors
grandmean = repmat(mean(mean(matrix1)),15,1);
individualmean = mean(matrix1,2);
indi_diff = individualmean - grandmean;
indi_diff = [indi_diff indi_diff];
matrix1_corr = matrix1 - indi_diff;

for i = 1:2
    grandmatrix(1,i) = mean(matrix1_corr(:,i));
    grandmatrix(2,i) = std(matrix1_corr(:,i))/sqrt(15);
end

subplot(1,2,1),d = bar(grandmatrix(1,:),'BarWidth',0.7);
ylabel('Percentage of accurate judgment');
ylim([0 1]);
title('c)','FontWeight','normal');
set(gca,'XTickLabel',{'Present','Absent'},'FontSize',12);
d.FaceColor = colours(2,:);
H = sigstar({{'Present','Absent'}},0.001);
hold on
errorbar(grandmatrix(1,:),grandmatrix(2,:),'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off




%hypothesis 2

colours = cbrewer('qual', 'Pastel1', 8); 
grandmatrix = zeros(2,2);
matrix1 = zeros(15,2);



for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP|Find_Incongruent_IP),:);
 
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)==3;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==2;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    end
    matrix1(sub,condition) = accuracy;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

% calculate standard errors
grandmean = repmat(mean(mean(matrix1)),15,1);
individualmean = mean(matrix1,2);
indi_diff = individualmean - grandmean;
indi_diff = [indi_diff indi_diff];
matrix1_corr = matrix1 - indi_diff;

for i = 1:2
    grandmatrix(1,i) = mean(matrix1_corr(:,i));
    grandmatrix(2,i) = std(matrix1_corr(:,i))/sqrt(15);
end

subplot(1,2,1),d = bar(grandmatrix(1,:),'BarWidth',0.7);
ylabel('Endorsement rate');
ylim([0 1]);
title('c)','FontWeight','normal');
set(gca,'XTickLabel',{'Present','Absent'},'FontSize',12);
d.FaceColor = colours(2,:);
hold on
errorbar(grandmatrix(1,:),grandmatrix(2,:),'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off



% hypothesis 3
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];

colours = cbrewer('qual', 'Set1', 8); 
grandmatrix = zeros(2,2);
grandmean1 = sum(Results(Find_CAP|Find_IAP|Find_N,8)==1)/length(Results(Find_CAP|Find_IAP|Find_N,:));

matrix1 = zeros(90,4);
location = [0 6.48 9.16];

for loc = 1:3
for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N)& Results(:,end)==location(loc),:);

for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
    end
    rowindex = size(matrix1,1)-sum(matrix1(:,1)==0)+1;
    
    matrix1(rowindex,1) = sub;
    matrix1(rowindex,2) = condition - 1;
    matrix1(rowindex,3) = location(loc);
    matrix1(rowindex,4) = accuracy;
    
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

end


% % fit lme
% data1 = table(matrix1(matrix1(:,2)==0,1),matrix1(matrix1(:,2)==0,3),matrix1(matrix1(:,2)==0,4),...
%     'VariableName',{'Subjects','Location','Accuracy'});
% lm1 = fitlme(data1,'Accuracy ~ Location  + (1|Subjects) + (Location|Subjects) ');
% lm2 = fitlme(data1,'Accuracy ~ 1+ (1|Subjects) + (Location|Subjects)');
% compare(lm2,lm1)
% 
% data2 = table(matrix1(matrix1(:,2)==1,1),matrix1(matrix1(:,2)==1,3),matrix1(matrix1(:,2)==1,4),...
%     'VariableName',{'Subjects','Location','Accuracy'});
% lm3 = fitlme(data2,'Accuracy ~ Location  + (1|Subjects) + (Location|Subjects) ');
% lm4 = fitlme(data2,'Accuracy ~ 1+ (1|Subjects) + (Location|Subjects)');
% compare(lm4,lm3)

% calculate standard errors for hit rate
matrix2 = [matrix1(matrix1(:,2)==0 & matrix1(:,3)==0,4) matrix1(matrix1(:,2)==0 & matrix1(:,3)==6.48,4)...
          matrix1(matrix1(:,2)==0 & matrix1(:,3)==9.16,4)];

grandmean1= repmat(mean(mean(matrix2)),15,1);
individualmean1 = mean(matrix2,2);
indi_diff1 = individualmean1 - grandmean1;
indi_diff1 = [indi_diff1 indi_diff1 indi_diff1];
matrix2_corr = matrix2 - indi_diff1;

for i = 1:3
    grandmatrix1(1,i) = mean(matrix2_corr(:,i));
    grandmatrix1(2,i) = std(matrix2_corr(:,i))/sqrt(15);
end

% calculate standard errors for CR
matrix3 = [matrix1(matrix1(:,2)==1 & matrix1(:,3)==0,4) matrix1(matrix1(:,2)==1 & matrix1(:,3)==6.48,4)...
          matrix1(matrix1(:,2)==1 & matrix1(:,3)==9.16,4)];

grandmean2= repmat(mean(mean(matrix3)),15,1);
individualmean2 = mean(matrix3,2);
indi_diff2 = individualmean2 - grandmean2;
indi_diff2 = [indi_diff2 indi_diff2 indi_diff2];
matrix3_corr = matrix3 - indi_diff2;

for i = 1:3
    grandmatrix2(1,i) = mean(matrix3_corr(:,i));
    grandmatrix2(2,i) = std(matrix3_corr(:,i))/sqrt(15);
end

subplot(1,2,1),b = errorbar([0 6.48 9.16],grandmatrix1(1,:),grandmatrix1(2,:),'o-',...
    'MarkerSize',0.8,'LineWidth',1.2);
b.Color = colours(2,:);
hold on
r = errorbar([0 6.48 9.16],grandmatrix2(1,:),grandmatrix2(2,:),'o-',...
    'MarkerSize',0.8,'LineWidth',1.2);
r.Color = colours(1,:);
hold off

set(gca,'FontSize',12);
ylim([0.4 1]), xticks([0 6.48 9.16]),ylabel('Percentage of accurate judgment'),xlabel('Eccentricity (in dva)'),...
    legend({'Patch present','Patch absent'}),title('a)','FontWeight','normal');
    

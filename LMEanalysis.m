clear all
%%% exp 1
Results = importdata('Pooled Results.mat');
sub_AUC = AUC_matrix(Results,1);

%%% exp 2
run('get_data2.m');
sub_AUC = AUC_matrix(Results,2);

%% analysis on criterion-free measures - hypo1, type 1

% exp 1 table
data_h1_t1 = table(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)~= -1,3),sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)~= -1,4),...
    sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)~= -1,1),...
    'VariableNames',{'eccentricity','AUC','subjects'});

%%% exp 2 table
% data_h1_t1 = table(categorical(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)~= -1,2)),sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)~= -1,4),...
%     sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)~= -1,1),...
%     'VariableNames',{'eccentricity','AUC','subjects'});

lme_h1_t1_full = fitlme(data_h1_t1,'AUC ~ eccentricity + (1|subjects)');
lme_h1_t1_reduced = fitlme(data_h1_t1,'AUC ~ 1+ (1|subjects)')
compare(lme_h1_t1_reduced,lme_h1_t1_full)


[h1,p1] = ttest(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)== 0,4),0.5,'alpha',0.05/3);
[h2,p2] = ttest(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)== 1,4),0.5,'alpha',0.05/3);
[h3,p3] = ttest(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2)== 2,4),0.5,'alpha',0.05/3);


%% hypo1, type 2
% exp 1 table
data_h1_t2 = table(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)~= -1,3),sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)~= -1,4),sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)~= -1,1),...
    'VariableNames',{'eccentricity','AUC','subjects'});

% % exp2 table
% data_h1_t2 = table(categorical(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)~= -1,2)),sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)~= -1,4),sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)~= -1,1),...
%     'VariableNames',{'eccentricity','AUC','subjects'});

lme_h1_t2_full = fitlme(data_h1_t2,'AUC ~ eccentricity + (1|subjects)');
lme_h1_t2_reduced = fitlme(data_h1_t2,'AUC ~ 1+ (1|subjects)')
compare(lme_h1_t2_reduced,lme_h1_t2_full)
% 
% [h,p] = ttest(sub_AUC.hypo1_Type2(:,1),0.5);
% LME.overall_h1_t2 = p;
% LME.quadratic_h1_t2 = compare(lme_h1_t2_quadratic,lme_h1_t2_full);
mean(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)== -1,4))
std(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2)== -1,4))/sqrt(240)
[h4,p4] = ttest(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type2(:,2)== 0,4),0.5,'alpha',0.05/3);
[h5,p5] = ttest(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type2(:,2)== 1,4),0.5,'alpha',0.05/3);
[h6,p6] = ttest(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type2(:,2)== 2,4),0.5,'alpha',0.05/3);

%% hypo2, type 1
%%% exp 1 table
data_h2_t1 = table(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,4),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,1),...
    categorical(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,2)),'VariableNames',{'eccentricity','AUC','subjects','condition'});

%%% exp 2 table
% data_h2_t1 = table(cagtegorical(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,3)),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,1),...
%     categorical(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1,2)),'VariableNames',{'eccentricity','AUC','subjects','condition'});
 
lme_h2_t1_full = fitlme(data_h2_t1,'AUC~eccentricity*condition + (1|subjects)');
lme_h2_t1_ecc = fitlme(data_h2_t1,'AUC~eccentricity:condition + condition + (1|subjects)');
lme_h2_t1_interaction = fitlme(data_h2_t1,'AUC~eccentricity + condition + (1|subjects)');
lme_h2_t1_condition = fitlme(data_h2_t1,'AUC~eccentricity + eccentricity:condition + (1|subjects)');

ecc_h2_t1 = compare(lme_h2_t1_ecc,lme_h2_t1_full)
interaction_h2_t1 = compare(lme_h2_t1_interaction, lme_h2_t1_full)
condition_h2_t1 = compare(lme_h2_t1_condition,lme_h2_t1_full)

mean(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 1,5)), std(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 1,5))/sqrt(240)
[h1,p1] = ttest(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 1,5),0.5,'alpha', 0.025)

mean(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 2,5)), std(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 2,5))/sqrt(240)
[h1,p1] = ttest(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 2,5),0.5,'alpha', 0.025)
[h2,p2] = ttest(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== -1 & sub_AUC.hypo2_Type1(:,2)== 2,5))

[h3,p3] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 0 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 0 & sub_AUC.hypo2_Type1(:,2)== 2,5),'alpha', 0.0167)
[h4,p4] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 1 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 1 & sub_AUC.hypo2_Type1(:,2)== 2,5),'alpha', 0.0167)
[h5,p5] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 2 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 2 & sub_AUC.hypo2_Type1(:,2)== 2,5),'alpha', 0.0167)

[h3,p3] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 0 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 1 & sub_AUC.hypo2_Type1(:,2)== 1,5),'alpha', 0.0167)
[h4,p4] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 1 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 2 & sub_AUC.hypo2_Type1(:,2)== 1,5),'alpha', 0.0167)
[h5,p5] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 0 & sub_AUC.hypo2_Type1(:,2)== 1,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 2 & sub_AUC.hypo2_Type1(:,2)== 1,5),'alpha', 0.0167)

[h6,p6] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 0 & sub_AUC.hypo2_Type1(:,2)== 2,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 1 & sub_AUC.hypo2_Type1(:,2)== 2,5),'alpha', 0.0167)
[h7,p7] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 1 & sub_AUC.hypo2_Type1(:,2)== 2,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 2 & sub_AUC.hypo2_Type1(:,2)== 2,5),'alpha', 0.0167)
[h8,p8] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 0 & sub_AUC.hypo2_Type1(:,2)== 2,5),sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)== 2 & sub_AUC.hypo2_Type1(:,2)== 2,5),'alpha', 0.0167)

%% hypo 2, type 2

%%% exp 1 table
data_h2_t2 = table(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,4),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,5),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,1),...
    categorical(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,2)),'VariableNames',{'eccentricity','AUC','subjects','condition'});

%%% exp 2 table
% data_h2_t2 = table(categorical(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,3)),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,5),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,1),...
%     categorical(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)~= -1,2)),'VariableNames',{'eccentricity','AUC','subjects','condition'});
 
lme_h2_t2_full = fitlme(data_h2_t2,'AUC~eccentricity*condition + (1|subjects)');
lme_h2_t2_ecc = fitlme(data_h2_t2,'AUC~eccentricity:condition + condition + (1|subjects)');
lme_h2_t2_interaction = fitlme(data_h2_t2,'AUC~eccentricity + condition + (1|subjects)');
lme_h2_t2_condition = fitlme(data_h2_t2,'AUC~eccentricity + eccentricity:condition + (1|subjects)');

ecc_h2_t2 = compare(lme_h2_t2_ecc,lme_h2_t2_full)
interaction_h2_t2 = compare(lme_h2_t2_interaction, lme_h2_t2_full)
condition_h2_t2 = compare(lme_h2_t2_condition,lme_h2_t2_full)

cong_h2t2 = mean(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 1,5))
se_h2t2 = std(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 1,5))/sqrt(240)
incong_h2t2 = mean(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 2,5))
incongse_h2t2 = std(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 2,5))/sqrt(240)

[h1,p1] = ttest(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 1,5),0.5,'alpha',0.025)
[h2,p2] = ttest(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 2,5),0.5,'alpha',0.025)
[h3,p3] = ttest(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 2,5),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== -1 & sub_AUC.hypo2_Type2(:,2)== 1,5))
[h4,p4] = ttest2(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== 0 & sub_AUC.hypo2_Type2(:,2)== 2,5),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== 0 & sub_AUC.hypo2_Type2(:,2)== 1,5),...
    'alpha',0.05/4)
[h5,p5] = ttest2(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== 1 & sub_AUC.hypo2_Type2(:,2)== 2,5),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== 1 & sub_AUC.hypo2_Type2(:,2)== 1,5),...
    'alpha',0.05/3)
[h6,p6] = ttest2(sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== 2 & sub_AUC.hypo2_Type2(:,2)== 2,5),sub_AUC.hypo2_Type2(sub_AUC.hypo2_Type2(:,3)== 2 & sub_AUC.hypo2_Type2(:,2)== 1,5),...
    'alpha',0.05/3)


%% Comparing Experiment 1 & 2 present vs. null differentation performance

% Type 1 AUC, across and on each eccentricity level
[h1,p1] = kstest2(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2) == -1, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type1(:,2) == -1, 4),'alpha', 0.0125)
[h2,p2] = kstest2(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2) == 0, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type1(:,2) == 0, 4),'alpha', 0.0125)
[h3,p3] = kstest2(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2) == 1, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type1(:,2) == 1, 4),'alpha', 0.0125)
[h4,p4] = kstest2(sub_AUC.hypo1_Type1(sub_AUC.hypo1_Type1(:,2) == 2, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type1(:,2) == 2, 4),'alpha', 0.0125)


% Type 2 AUC, across and on each eccentricity level
[h1,p1] = kstest2(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2) == -1, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type2(:,2) == -1, 4),'alpha', 0.0125)
[h2,p2] = kstest2(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2) == 0, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type2(:,2) == 0, 4),'alpha', 0.0125)
[h3,p3] = kstest2(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2) == 1, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type2(:,2) == 1, 4),'alpha', 0.0125)
[h4,p4] = kstest2(sub_AUC.hypo1_Type2(sub_AUC.hypo1_Type2(:,2) == 2, 4),sub_AUC_1.hypo1_Type1(sub_AUC_1.hypo1_Type2(:,2) == 2, 4),'alpha', 0.0125)


% Type 1 AUC, across and on each eccentricity level
[h1,p1] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 1 & sub_AUC.hypo2_Type1(:,3) == -1, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 1 &...
    sub_AUC_1.hypo2_Type1(:,3) == -1, 5),'alpha', 0.0125)
[h2,p2] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 1 & sub_AUC.hypo2_Type1(:,3) == 0, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 1 &...
    sub_AUC_1.hypo2_Type1(:,3) == 0, 5),'alpha', 0.0125)
[h3,p3] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 1 & sub_AUC.hypo2_Type1(:,3) == 1, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 1 &...
    sub_AUC_1.hypo2_Type1(:,3) == 1, 5),'alpha', 0.0125)
[h4,p4] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 1 & sub_AUC.hypo2_Type1(:,3) == 2, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 1 &...
    sub_AUC_1.hypo2_Type1(:,3) == 2, 5),'alpha', 0.0125)

[h1,p1] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 2 & sub_AUC.hypo2_Type1(:,3) == -1, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 2 &...
    sub_AUC_1.hypo2_Type1(:,3) == -1, 5),'alpha', 0.0125)
[h2,p2] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 2 & sub_AUC.hypo2_Type1(:,3) == 0, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 2 &...
    sub_AUC_1.hypo2_Type1(:,3) == 0, 5),'alpha', 0.0125)
[h3,p3] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 2 & sub_AUC.hypo2_Type1(:,3) == 1, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 2 &...
    sub_AUC_1.hypo2_Type1(:,3) == 1, 5),'alpha', 0.0125)
[h4,p4] = kstest2(sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,2) == 2 & sub_AUC.hypo2_Type1(:,3) == 2, 5),sub_AUC_1.hypo2_Type1(sub_AUC_1.hypo2_Type1(:,2) == 2 &...
    sub_AUC_1.hypo2_Type1(:,3) == 2, 5),'alpha', 0.0125)
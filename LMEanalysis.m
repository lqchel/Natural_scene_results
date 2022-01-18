clear all
Results = get_data2(2);
sub_AUC = AUC_matrix(Results);

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

% interaction between condition and eccentricity significant, so fit an lme
% model to congruent and incongruent separately

% congruent
data = sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1 & sub_AUC.hypo2_Type1(:,2) == 1,:);
AUC_lme([data(:,1) data(:,3:4)],1,1)

% incongruent
data = sub_AUC.hypo2_Type1(sub_AUC.hypo2_Type1(:,3)~= -1 & sub_AUC.hypo2_Type1(:,2) == 2,:);
AUC_lme([data(:,1) data(:,3:4)],1,1)

%% check no responses no modifid patches only

subject_id = unique(Results(:,1));

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
index{1} = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1; % congruent, modified
index{2} = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1; % incongruent, original



location = [0 1 2];
accurate_judgement = [1 -1 1 -1]; % yes to original, no to modified

% congruent
for in = 1:2
    condition_data = Results(index{in},:);
    b = 1;
    for loc = 1:3
        for sub = 1:length(subject_id)
                current_sub_data = condition_data(condition_data(:,1) == subject_id(sub) & condition_data(:,13) == location(loc),:);
                mt(b,:)= [sub location(loc)  sum(current_sub_data(:,8) == accurate_judgement(in))/length(current_sub_data(:,8))]; 
                b = b+1;
        end
    end
    AUC_lme(mt,1)
    clear mt
end



% condition_data = Results(index{2},:);
% b = 1;
% for loc = 1:3
%     for sub = 1:length(subject_id)
%         current_sub_data = condition_data(condition_data(:,1) == subject_id(sub) & condition_data(:,13) == location(loc),:);
%         mt(b,:)= [sub location(loc) sum(current_sub_data(:,8) == -1)/length(current_sub_data(:,8))];
%         b = b+1;
%     end
% end
% AUC_lme(mt,1)

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



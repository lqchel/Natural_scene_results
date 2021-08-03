clear all

sub_AUC = AUC_matrix([b1_g2;b1_g4;b2_g1;b2_g4]);

%% analysis on criterion-free measures - hypo1, type 1

eccentricity = [zeros(22,1); 6.5 + zeros(22,1); 9.2 + zeros(22,1)];
sub_number = reshape(1:22,[22,1]);
subjects = [sub_number; sub_number; sub_number];
data_h1_t1 = table(eccentricity,[sub_AUC.hypo1_Type1(:,2);sub_AUC.hypo1_Type1(:,3);sub_AUC.hypo1_Type1(:,4)],subjects,'VariableNames',{'eccentricity','AUC','subjects'});
lme_h1_t1_full = fitlme(data_h1_t1,'AUC ~ eccentricity + (1|subjects)');
lme_h1_t1_reduced = fitlme(data_h1_t1,'AUC ~ (1|subjects)')
% lme_h1_t1_quadratic = fitlme(data_h1_t1,'AUC ~ 1 + (1|subjects)');


[h,p] = ttest(sub_AUC.hypo1_Type1(:,1),0.5);
LME.overall_h1_t1 = p;
LME.quadratic_h1_t1 = compare(lme_h1_t1_quadratic,lme_h1_t1_full);
[h1,p1] = ttest(sub_AUC.hypo1_Type1(:,2),sub_AUC.hypo1_Type1(:,3),'alpha',0.05/3);
[h2,p2] = ttest(sub_AUC.hypo1_Type1(:,3),sub_AUC.hypo1_Type1(:,4),'alpha',0.05/3);
[h3,p3] = ttest(sub_AUC.hypo1_Type1(:,2),sub_AUC.hypo1_Type1(:,4),'alpha',0.05/3);


%% hypo1, type 2

data_h1_t2 = table(eccentricity,[sub_AUC.hypo1_Type2(:,2);sub_AUC.hypo1_Type2(:,3);sub_AUC.hypo1_Type2(:,4)],subjects,'VariableNames',{'eccentricity','AUC','subjects'});
lme_h1_t2_full = fitlme(data_h1_t2,'AUC ~ eccentricity^2 + (1|subjects)');
lme_h1_t2_quadratic = fitlme(data_h1_t2,'AUC ~ 1 + (1|subjects)');

[h,p] = ttest(sub_AUC.hypo1_Type2(:,1),0.5);
LME.overall_h1_t2 = p;
LME.quadratic_h1_t2 = compare(lme_h1_t2_quadratic,lme_h1_t2_full);

[h4,p4] = ttest(sub_AUC.hypo1_Type2(:,2),sub_AUC.hypo1_Type2(:,3),'alpha',0.05/3);
[h5,p5] = ttest(sub_AUC.hypo1_Type2(:,3),sub_AUC.hypo1_Type2(:,4),'alpha',0.05/3);
[h6,p6] = ttest(sub_AUC.hypo1_Type2(:,2),sub_AUC.hypo1_Type2(:,4),'alpha',0.05/3);

%% hypo2, type 1

[h1,p1] = ttest(sub_AUC.hypo2_Type1_cong(:,1),0.5,'alpha',0.05/3);
[h2,p2] = ttest(sub_AUC.hypo2_Type1_incong(:,1),0.5,'alpha',0.05/3);
[h3,p3] = ttest(sub_AUC.hypo2_Type1_incong(:,1),sub_AUC.hypo2_Type1_cong(:,1),'alpha',0.05/3);
LME.overall_h2_t1 = [p1 p2 p3];

data_h2_t1 = table([eccentricity;eccentricity],[reshape(sub_AUC.hypo2_Type1_cong(:,2:4),[66,1]); reshape(sub_AUC.hypo2_Type1_incong(:,2:4),...
    [66,1])],[subjects;subjects],[ones(66,1);1+ones(66,1)],'VariableNames',{'eccentricity','AUC','subjects','condition'});

lme_h2_t1_full = fitlme(data_h2_t1,'AUC~eccentricity*condition + (1|subjects)');
lme_h2_t1_ecc = fitlme(data_h2_t1,'AUC~eccentricity:condition + condition + (1|subjects)');
lme_h2_t1_interaction = fitlme(data_h2_t1,'AUC~eccentricity + condition + (1|subjects)');
lme_h2_t1_condition = fitlme(data_h2_t1,'AUC~eccentricity + eccentricity:condition + (1|subjects)');

LME.ecc_h2_t1 = compare(lme_h2_t1_ecc,lme_h2_t1_full)
LME.interaction_h2_t1 = compare(lme_h2_t1_interaction, lme_h2_t1_full)
LME.condition_h2_t1 = compare(lme_h2_t1_condition,lme_h2_t1_full)





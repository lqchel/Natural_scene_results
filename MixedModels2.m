clear all
Results = importdata('Pooled & Cleaned Results.mat');
Results(:,9) = Results(:,9)-0.5;
Results(:,9) = Results(:,8).*Results(:,9);
% signal present and absent trial classification

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
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
% 

%% lme for hypothesis 3-1
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];


Results_H1 = Results(Find_CAP|Find_IAP|Find_N,:);
APN = Results_H1(:,5)~= 1;
congruency_H1 = categorical(Results_H1(:,4));
eccentricity_H1 = Results_H1(:,end);
Response_H1 = Results_H1(:,9);
subjects_H1 = Results_H1(:,1);
dataset_H1 = table(APN, eccentricity_H1, Response_H1, subjects_H1,congruency_H1, Results_H1(:,2),'VariableName',{'PatchType', 'Eccentricity','Response','Participants','Congruency','Image'});
lm1 = fitlme(dataset_H1,'Response~ PatchType * Eccentricity + (1 | Participants) + (1 | Image)');

%% likelihood ratio test
lm2 = fitlme(dataset_H1,'Response~ PatchType + Eccentricity + (1 | Participants) + (1 | Image)');
lm3 = fitlme(dataset_H1,'Response~ PatchType + PatchType:Eccentricity + (1 | Participants) + (1 | Image)');
lm4 = fitlme(dataset_H1,'Response~ Eccentricity + PatchType*Eccentricity + (1 | Participants) + (1 | Image)');
lm5 = fitlme(dataset_H1,'Response~ PatchType * Eccentricity + (1|Image)');
lm6 = fitlme(dataset_H1,'Response~ PatchType * Eccentricity + (1|Participants)');

test_interaction = compare(lm2,lm1)
test_eccentricity = compare(lm3,lm1)
test_APN = compare(lm4,lm1)
test_randsubject = compare(lm5,lm1)
test_randimage = compare(lm6,lm1)

%% assumption checks
histogram(residuals(lm1));
plot(fitted(lm1),residuals(lm1),'bx'), xlabel('Predicted values'),ylabel('Residuals'),title('Residuals against predicted values'); % violated
histogram(randomEffects(lm1));

%% visualisation of effects
yfit1 = fitted(lm1,'Conditional',false);
yfit1r = fitted(lm1);
gscatter(dataset_H1.Eccentricity,yfit1,dataset_H1.PatchType);

%% lme for hypothesis 3-2

Results_H2 = Results(Find_Congruent_CP|Find_Incongruent_IP|Find_Incongruent_CP|Find_Congruent_IP,:);
eccentricity_H2 = Results_H2(:,end);
PN = (Results_H2(:,4)==0 & Results_H2(:,5)==2)|(Results_H2(:,4)==1 & Results_H2(:,5)==3);
dataset_H2 = table(PN, eccentricity_H2, Results_H2(:,1), Results_H2(:,2),Results_H2(:,4),Results_H2(:,9),'VariableName',{'PatchType','Eccentricity','Subjects','Image','Congruency','Response'});
lme1 = fitlme(dataset_H2,'Response ~ PatchType*Eccentricity*Congruency + (1|Subjects)+(1|Image)');
lme2 = fitlme(dataset_H2, 'Response ~ PatchType*Eccentricity + (1|Congruency)+(1|Subjects) + (1|Image)')

%% testing the model

lme_mPatch = fitlme(dataset_H2,'Response ~ Eccentricity + Congruency + PatchType:Eccentricity + PatchType:Congruency + Eccentricity:Congruency + PatchType:Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_mEccentricity = fitlme(dataset_H2,'Response ~ PatchType + Congruency + PatchType:Eccentricity + PatchType:Congruency + Eccentricity:Congruency + PatchType:Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_mCongruency = fitlme(dataset_H2,'Response ~ PatchType + Eccentricity + PatchType:Eccentricity + PatchType:Congruency + Eccentricity:Congruency + PatchType:Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_mPE = fitlme(dataset_H2,'Response ~ PatchType + Eccentricity + Congruency + PatchType:Congruency + Eccentricity:Congruency + PatchType:Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_mPC = fitlme(dataset_H2,'Response ~ PatchType + Eccentricity + Congruency + PatchType:Eccentricity + Eccentricity:Congruency + PatchType:Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_mEC = fitlme(dataset_H2,'Response ~ PatchType + Eccentricity + Congruency + PatchType:Eccentricity + PatchType:Congruency + PatchType:Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_3intr = fitlme(dataset_H2,'Response ~ PatchType + Eccentricity + Congruency + PatchType:Eccentricity + PatchType:Congruency + Eccentricity:Congruency + (1|Subjects) + (1|Image)');
lme_rsubject = fitlme(dataset_H2,'Response ~ PatchType*Eccentricity*Congruency +(1|Image)');
lme_rimage = fitlme(dataset_H2,'Response ~ PatchType*Eccentricity*Congruency +(1|Subjects)' );

test_mPatch = compare(lme_mPatch,lme1)
test_mEccentricity = compare(lme_mEccentricity,lme1)
test_mCongruency = compare(lme_mCongruency,lme1)
test_mPE = compare(lme_mPE, lme1)
test_mPC = compare(lme_mPC, lme1)
test_3intr = compare(lme_3intr,lme1)











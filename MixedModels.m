%% LME model analysis
clear all
Results = importdata('Pooled & Cleaned Results.mat');
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
%% compute accuracy, I don't think it makes sense so should be discarded?
% Find_AP = Find_CAP | Find_IAP;
% CorrectY = Find_AP & Results(:,8)==1;
% Find_AllN = Find_N|Find_Congruent_IP|Find_Incongruent_CP;
% CorrectN = Find_AllN & Results(:,8)==-1;
% All_Correct = CorrectY | CorrectN;
% All_false = -1*~All_Correct;

% accuracy = All_Correct + All_false;
% conf_accuracy = Results(:,9).*accuracy;
% histogram(conf_accuracy)

%% patch difference analysis

PD_level = categorical(Find_N(~(Find_CAP|Find_IAP)));
dataset3 = table(Results(~(Find_CAP|Find_IAP),11),PD_level,'VariableNames',{'PatchDifference','PatchLevel'});
lm_PD = fitlme(dataset3,'PatchDifference ~ PatchLevel');
plotResiduals(lm_PD,'histogram')
histogram(PD_level, Results(~(Find_CAP|Find_IAP),11));

%% building LME model for patch difference predicting responses

central = Results(:,7)==5;
periphery_1= 2.*(Results(:,7)==2|Results(:,7)==4|Results(:,7)==6|Results(:,7)==8);
periphery_2 = 3.* (Results(:,7)==1|Results(:,7)==3|Results(:,7)==7|Results(:,7)==9);
location = categorical(central + periphery_1 + periphery_2,[1 2 3],{'central','periphery_1','periphery_2'});

subject = Results(~(Find_CAP|Find_IAP),1);
image = Results(~(Find_CAP|Find_IAP),2);
congruency = categorical(Results(~(Find_CAP|Find_IAP),4));

PD1 = Find_Incongruent_CP;
PD2 = 2.* Find_N;
PD_level = PD1 + PD2;
PD_level = categorical(PD_level(~(Find_CAP|Find_IAP)));


Response = Results(~(Find_CAP|Find_IAP),8)==-1;
dataset2 = table(Response,Results(~(Find_CAP|Find_IAP),11),congruency,subject,image,PD_level,'VariableNames',{'Response','PatchDifference',...
    'Congruency','Participant','Image','PatchType'});


glm1 = fitglme(dataset2,'Response~PatchDifference * PatchType +(1|Participant)+ (1|Image)','Distribution','Binomial','FitMethod','MPL',...
    'Link','logit');
glm1



glm2 = fitglme(dataset2,'Response~PatchDifference +(1|Participant)+ (1|Image)','Distribution','Binomial','FitMethod','MPL',...
    'Link','logit');
glm2



glm3 = fitglme(dataset2,'Response~PatchDifference + PatchType +(1|Participant)+ (1|Image)','Distribution','Binomial','FitMethod','REMPL',...
    'Link','logit');
glm3

glm4 = fitglme(dataset2,'Response~PatchDifference + Congruency +(PatchDifference|Participant)+ (1|Image)','Distribution','Binomial','FitMethod','MPL',...
    'Link','logit');
glm4

%% lme, N patches only
subject = Results(Find_N,1);
image = Results(Find_N,2);
dataset = table(Results(Find_N,9),Results(Find_N,11),subject,image,'VariableNames',{'Response','PatchDifference',...
    'Participant','Image'});

lm3 = fitlme(dataset,'Response~PatchDifference +(1|Participant)+ (1|Image)');
lm3

yfit3 = fitted(lm3, 'Conditional', false);
yfit3r = fitted(lm3);
s = scatter(dataset.PatchDifference,yfit3r,'b.');
title('Response x confidence across patch difference values'),...
xlabel('Patch difference'),ylabel('Response x confidence'),ylim([-4 4]),yticks([-4:1:4]);
s.MarkerFaceAlpha = 0.3;
hold on
plot(dataset.PatchDifference,yfit3,'r-','LineWidth',2.5);
hold off
%% lme, incongruent and congruent objects only
subject = Results(Find_Incongruent_CP|Find_Congruent_IP,1);
image = Results(Find_Incongruent_CP|Find_Congruent_IP,2);
dataset2 = table(Results(Find_Incongruent_CP|Find_Congruent_IP,9),Results(Find_Incongruent_CP|Find_Congruent_IP,11),subject,image,'VariableNames',{'Response','PatchDifference',...
    'Participant','Image'});

lm4 = fitlme(dataset2,'Response~PatchDifference + (1|Participant)+ (1|Image)');
lm4

yfit4 = fitted(lm4, 'Conditional', false);
yfit4r = fitted(lm4);
s = scatter(dataset2.PatchDifference,yfit4r,'b.');
title('Response x confidence across patch difference values for critical objects'),...
xlabel('Patch difference'),ylabel('Response x confidence');
s.MarkerFaceAlpha = 0.3;
hold on
plot(dataset2.PatchDifference,yfit4,'r-','LineWidth',2.5);
hold off

r = Results(Find_Incongruent_CP|Find_Congruent_IP,9);

%% incongruent object
subject = Results(Find_Incongruent_CP,1);
image = Results(Find_Incongruent_CP,2);
dataset3 = table(Results(Find_Incongruent_CP,9),Results(Find_Incongruent_CP,11),subject,image,'VariableNames',{'Response','PatchDifference',...
    'Participant','Image'});

lm5 = fitlme(dataset3,'Response~PatchDifference + (1|Participant)+ (1|Image)');
lm5

yfit5 = fitted(lm5, 'Conditional', false);
yfit5r = fitted(lm5);
s = scatter(dataset3.PatchDifference,yfit5r,'b.');
title('Response x confidence across patch difference values for incongruent objects'),...
xlabel('Patch difference'),ylabel('Response x confidence');
s.MarkerFaceAlpha = 0.3;
hold on
plot(dataset3.PatchDifference,yfit5,'r-','LineWidth',2.5);
hold off


%% lme, for both N and incongruent/congruent
subject = Results(~(Find_CAP|Find_IAP),1);
image = Results(~(Find_CAP|Find_IAP),2);
congruency = categorical(Results(~(Find_CAP|Find_IAP),4));

PD1 = Find_Incongruent_CP;
PD2 = 2.* Find_N;
PD_level = PD1 + PD2;
PD_level = categorical(PD_level(~(Find_CAP|Find_IAP)));

dataset = table(Results(~(Find_CAP|Find_IAP),9),Results(~(Find_CAP|Find_IAP),11),subject,image,PD_level,'VariableNames',{'Response','PatchDifference',...
    'Participant','Image','PatchLevel'});

lm1 = fitlme(dataset,'Response~PatchDifference*PatchLevel +(1|Participant)+ (1|Image)');
lm1

yfit1 = fitted(lm1,'Conditional',false);
gscatter(dataset.PatchDifference,yfit1,dataset.PatchLevel,'rgb','.',12,'off'),ylim([-3 2]),...
xlabel('Patch difference in sum RGB values'), ylabel('Predicted response x confidence values'), title('Fixed-effect prediction of Model');

yfit2=(fitted(lm1));
gscatter(dataset.PatchDifference,yfit2,dataset.PatchLevel,'rgb','.',3,'off'),...
xlabel('Patch difference in sum RGB values'), ylabel('Predicted response x confidence values'), title('Model prediction accounting for random effects');

lm1_res = residuals(lm1);
plot(yfit2,lm1_res,'bx'), xlabel('Predicted values'),ylabel('Residuals'),title('Residuals against predicted values');

lm3 = fitlme(dataset,'Response~PatchDifference +(1|Participant)+ (1|Image)');
lm3

lm4 = fitlme(dataset,'Response~PatchDifference + PatchLevel +(1|Participant)+ (1|Image)');
lm4


%% lme for hypothesis 1
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
dataset_H1 = table(APN, eccentricity_H1, Response_H1, subjects_H1,congruency_H1,'VariableName',{'PatchType', 'Eccentricity','Response','Participants','Congruency'});
lm1 = fitlme(dataset_H1, 'Response ~ PatchType + PatchType:Congruency + (PatchType|Participants)+ (1|Participants)');

%% model comparison for hypothesis 1
lm2 = fitlme(dataset_H1, 'Response ~ PatchType  + (PatchType|Participants)+ (1|Participants)');
lm3 = fitlme(dataset_H1,'Response ~ PatchType:Congruency + (PatchType|Participants)+ (1|Participants)');
lm4 = fitlme(dataset_H1,'Response ~ PatchType + PatchType:Congruency + (1|Participants)');
lm5 = fitlme(dataset_H1,'Response ~ PatchType + PatchType:Congruency + (PatchType|Participants)');

test_APN = compare(lm3,lm1);
test_Congruency = compare(lm2,lm1);
test_randslope = compare(lm4,lm1);
test_randintcpt = compare(lm5,lm1);

%% assumption checks
histogram(residuals(lm1));
plot(fitted(lm1),residuals(lm1),'bx'), xlabel('Predicted values'),ylabel('Residuals'),title('Residuals against predicted values'); % violated

%% lm3 for hypothesis 3-1
lm6 = fitlme(dataset_H1,'Response~ PatchType * Eccentricity + (PatchType|Participants) + (Eccentricity|Participants) + (1|Eccentricity)');

%% model comparison for hypothesis 3-1
lm7 = fitlme(dataset_H1,'Response~ PatchType + Eccentricity + (PatchType|Participants)+ (Eccentricity|Participants) + (1|Eccentricity)');
lm8 = fitlme(dataset_H1,'Response~ PatchType + PatchType:Eccentricity + (PatchType|Participants)+ (Eccentricity|Participants) + (1|Eccentricity)');
lm9 = fitlme(dataset_H1,'Response~ PatchType * Eccentricity + (PatchType|Participants) + (1|Eccentricity)');
lm10 = fitlme(dataset_H1,'Response~ PatchType * Eccentricity + (PatchType|Participants) + (Eccentricity|Participants)');
test_interaction = compare(lm7,lm6)
test_eccentricity = compare(lm8,lm6)
test_erandslope = compare(lm9,lm6)
test_erandincpt = compare(lm10,lm6)

%% assumption checks
histogram(residuals(lm6));
plot(fitted(lm6),residuals(lm6),'bx'), xlabel('Predicted values'),ylabel('Residuals'),title('Residuals against predicted values'); % violated

lme = glmfit(dataset_H1)

%% hypothesis 3-1
 
Response_H31 = (Results_H1(:,8)==1).*1;
dataset_H31 = table(APN, eccentricity_H1, Response_H31, subjects_H1,congruency_H1,Results_H1(:,2),'VariableName',{'PatchType', 'Eccentricity','Response','Participants','Congruency','Image'});
glm1 = fitglme(dataset_H31,'Response ~ PatchType*Eccentricity +(1|Participants)+ (1|Image)','Distribution','Binomial','FitMethod','MPL');

%% model comparison
glm2 = fitglme(dataset_H31,'Response ~ PatchType + (1|Participants)','Distribution','Binomial','FitMethod','MPL');
glm3=fitglme(dataset_H31,'Response ~ (1|Participants)','Distribution','Binomial','FitMethod','MPL');







%% data processing
clear all;
load('pilot_ql_01.mat');

PresentResponse = response > 0;
AbsentResponse_Index = response<0;
AbsentResponse = -1.* AbsentResponse_Index;
Response = PresentResponse + AbsentResponse;
Confidence = response;

NPatch = PatchFeature == 'n';
PPatch = PatchFeature == 'p';
PPatch = PPatch .* 3;
APatch = PatchFeature == 'a';
APatch = APatch .* 2;
PatchFeature = NPatch + PPatch + APatch;
PatchFeature = PatchFeature(1:1698,:);

PatchSource_S = PatchSource == 'S';
PatchSource_C = PatchSource == 'C';
PatchSource_C = PatchSource_C .* 2;
PatchSource_I = PatchSource == 'I';
PatchSource_I = PatchSource_I .* 3;
PatchSource = PatchSource_S + PatchSource_C + PatchSource_I;
PatchSource = PatchSource(1:1698,:);

ImageFeature = ImageFeature == "Incongruent";
ImageFeature = 1.*ImageFeature(1:1698,:);

CurrentTrialNumber = CurrentTrialNumber(1:1698,:);
QuestionNumber = QuestionNumber(1:1698,:);

Patch_location = Patch_location(1:1698,:);
for i = 1:length(Patch_location)
    Patch_Location(i,:) = str2num(Patch_location(i,:));
end

%PatchSource: shinji == 1, congruent ==2, incongruent == 3
%ImageFeature for the current trial: congruent ==0, incongruent == 1;
%PatchFeature: shiji == 1, critical object absent == 2, critical object
%present == 3
%Response: Present = 1, absent = -1
%

Results = [CurrentTrialNumber QuestionNumber ImageFeature PatchSource PatchFeature Patch_Location Response Confidence ResponseTime];
% load('pilotresult_RMG_01_01.mat');
load('RMG Alternative.mat');
Results = Trialresponse;
Results(:,8) = Results(:,8).*Results(:,7);

%by_congruency =input('Analyse by congruency of trials? (0/1) \n','s');

%% load data files
for sub = 1:15
%        for ses = 1:2
%            for block = 1:2
             folderE = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\results and analysis\exp results";
            filePatternE = fullfile(folderI, '*.mat');
            theFilesE = struct2cell(dir(filePatternI));
            selectnameE(1:60,1) = string(theFilesI(1,:));
            
            folderT = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\results and analysis\exp results";
            filePatternT = fullfile(folderI, '*.mat');
            theFilesT = struct2cell(dir(filePatternI));
            selectnameT(1:60,1) = string(theFilesI(1,:));
            
            ExpRst_1 = importdata(selectname((sub-1).*4 + 1,1));
            ExpRst_1 = ExpRst_1(ExpRst_1(:,2)~= 1 && ExpRst_1(:,2) ~= 2,:); 
            ExpRst_2 = importdata(selectname((sub-1).*4 + 2,1));
            ExpRst_3 = importdata(selectname((sub-1).*4 + 3,1));
            ExpRst_3 = ExpRst_3(ExpRst_3(:,2)~= 1 && ExpRst_3(:,2) ~= 2,:);
            ExpRst_4 = importdata(selectname((sub-1).*4 + 4,1));
            
               
%            end
%        end
end

%% signal present and absent trial classification

%Find trials that presented N patches -- signal absent for hypo 1
Find_N = Results(:,4) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,3) == 0 & Results(:,4) == 2; 
Find_IAP = Results(:,3) == 1 & Results(:,4) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,3) == 0 & Results(:,4) == 2 & Results(:,5) == 3;
Find_Incongruent_IP = Results(:,3) == 1 & Results(:,4) == 3 & Results(:,5) == 3;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP = Results(:,3) == 0 & Results(:,4) == 3 & Results(:,5) == 3;
Find_Incongruent_CP = Results(:,3) == 1 & Results(:,4) == 2 & Results(:,5) == 3;

%% hypo 1, not by congruency
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);
% type 1 ROC curve analysis

Confidence_N = Results_N(:,8);
Confidence_AP = Results_AP(:,8);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

%     subplot(2,2,1),plot([0:1:8],Confidence_NCounts,'b-','LineWidth',1.5),xlabel('Criterion (Strict to lenient)'),ylabel('Frequency'),title('Critical objects discrimination')
%     hold on
%     plot([0:1:8],Confidence_APCounts,'r-','LineWidth',1.5)
%     legend({'N', 'AP'});
%     hold off
% 
%     subplot(2,2,2),plot([0:1:8],Cumulative_NCounts,'LineWidth',1.5),xlabel('Criterion (strict to lenient)'),ylabel('Cumulative Frequency'),title('Critical objects discrimination')
%     hold on
%     plot([0:1:8],Cumulative_APCounts,'LineWidth',1.5)
%     legend({'N Patch','A and P Patches'});
%     hold off
% 
%     subplot(2,2,3),plot([0:1:8],Cumulative_FA,'LineWidth',1.5),xlabel('Criterion (strict to lenient)'),ylabel('Cumulative Probability'),title('Critical objects discrimination')
%     hold on
%     plot([0:1:8],Cumulative_Hit,'LineWidth',1.5);
%     legend({'N Patch','A and P Patches'})
%     hold off

subplot(3,2,1),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC curve for A/P patch discrimination AUC = ', num2str(AUC)]);
axis square

Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];

%% hypo 2, not by congruency
Results_N = Results(Find_Congruent_IP | Find_Incongruent_CP,:);
Results_AP = Results(Find_Congruent_CP | Find_Incongruent_IP,:);

Confidence_N = Results_N(:,8);
Confidence_AP = Results_AP(:,8);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(3,2,2),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['RMG ROC curve for critical objects discrimination AUC = ', num2str(AUC)]);
axis square


Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];


%% hypo 1, by congruency
%Congruent trials
Results_N = Results(Find_N & Results(:,3)== 0,:);
Results_AP = Results(Find_CAP,:);

Confidence_N = Results_N(:,8);
Confidence_AP = Results_AP(:,8);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(3,2,3),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['A/P discrimination in congruent trials AUC = ', num2str(AUC)]);
axis square


Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];

%incongruent trials
Results_N = Results(Find_N & Results(:,3)==1,:);
Results_AP = Results(Find_IAP,:);

Confidence_N = Results_N(:,8);
Confidence_AP = Results_AP(:,8);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(3,2,4),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['A/P discrimination in incongruent trials AUC = ', num2str(AUC)]);
axis square


Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];

%% hypo 2, by congruency
%congruent trials
Results_AP = Results(Find_Congruent_CP,:);
Results_N =  Results(Find_Congruent_IP,:);

Confidence_N = Results_N(:,8);
Confidence_AP = Results_AP(:,8);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(3,2,5),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['Critical objects discrimination in congruent trials AUC = ', num2str(AUC)]);
axis square


Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];

%Incongruent trials
Results_N = Results(Find_Incongruent_CP,:);
Results_AP = Results(Find_Incongruent_IP,:);

Confidence_N = Results_N(:,8);
Confidence_AP = Results_AP(:,8);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(3,2,6),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['Critical objects discrimination in incongruent trials AUC = ', num2str(AUC)]);
axis square

%% confidence ratings for hypo 1
%not by trial type
Find_CAP_Cor = Results_AP(:,3) == 0 & Results_AP(:,4) == 2 & Results_AP(:,7) == 1;
Find_IAP_Cor = Results_AP(:,3) == 1 & Results_AP(:,4) == 3 & Results_AP(:,7) == 1;
accuracy_collapsed1 = [Find_CAP_Cor Find_IAP_Cor];
for i = 1: length(Results_AP)
    if any(accuracy_collapsed1(i,:)==1)
        Trial_Accuracy1(i,:) = 1;
    else
        Trial_Accuracy1(i,:) = -1;
    end
end

Confidence_AP = Results_AP(:,8) .* Trial_Accuracy1;


Find_N_Cor = Results_N(:,7)== -1;
for a = 1:length(Results_N)
    if Find_N_Cor(a,:) == 1
        Trial_Accuracy2(a,:) = 1;
    else
        Trial_Accuracy2(a,:) = -1;
    end
end

Confidence_N = Results_N(:,8) .* Trial_Accuracy2;

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == i);
    Confidence_NCounts(i+5) = sum(Confidence_N == i);
end
plot([-4:1:4],Confidence_NCounts,'b-','LineWidth',1.5),xlim([-4 4]), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for absent patches')

plot([-4:1:4],Confidence_APCounts,'r-','LineWidth',1.5),xlim([-4 4]), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for present patches')

%by trial type
%congruent trials
Results_CAP = Results(Find_CAP,:);
Find_Cong_CAP = Results_CAP(:,7) == 1;
for b = 1:length(Results_CAP)
    if Find_Cong_CAP(b,:) == 1
        Trial_Accuracy3(b,:) = 1;
    else
        Trial_Accuracy3(b,:) = -1;
    end
end
Confidence_CAP = Results_CAP(:,8).*Trial_Accuracy3;

Results_CN = Results(Results(:,3) == 0 & Results(:,4) == 1,:);
Find_Cong_CN = Results_CN(:,7) == -1;
for c = 1:length(Results_CN)
    if Find_Cong_CN(c,:) == 1
        Trial_Accuracy4(c,:) = 1;
    else
        Trial_Accuracy4(c,:) = -1;
    end
end
Confidence_CN = Results_CN(:,8).*Trial_Accuracy4;

%incongruent trials
Results_IAP = Results(Find_IAP,:);
Find_Incong_IAP = Results_IAP(:,7) == 1;
for d = 1:length(Results_IAP)
    if Find_Incong_IAP(d,:) == 1
        Trial_Accuracy5(d,:) = 1;
    else
        Trial_Accuracy5(d,:) = -1;
    end
end
Confidence_IAP = Results_IAP(:,8).*Trial_Accuracy5;

Results_IN = Results(Results(:,3) == 1 & Results(:,4) == 1,:);
Find_Incong_IN = Results_IN(:,7) == -1;
for e = 1:length(Find_Incong_IN)
    if Find_Incong_IN(e,:) == 1
        Trial_Accuracy6(e,:) = 1;
    else
        Trial_Accuracy6(e,:) = -1;
    end
end
Confidence_IN = Results_IN(:,8).*Trial_Accuracy6;

subplot(2,2,1),histogram(Confidence_CAP,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for present patches in congruent trials')
subplot(2,2,2),histogram(Confidence_CN,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combinedfor absent patches in congruent trials')
subplot(2,2,3),histogram(Confidence_IAP,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combinedfor present patches in incongruent trials')
subplot(2,2,4),histogram(Confidence_IN,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for absent patches in incongruent trials')

%testing hypo 2
% x = categorical({'Correct responses','Incorrect responses'});
% d = [mean(All_Correct) mean(All_False)];
% error = [std(All_Correct) std(All_False)];
% bar(x,d,0.60), ylabel('Mean Confidence Ratings', 'FontSize', 14), title('Mean confidence ratings for correct and incorrect decisions')
% hold on
% errorbar(x,d,error,'k.','linewidth',0.75)
% hold off
% 
% subplot(2,1,1),hist(All_Correct), ylim([0 1550]),xticks(1:1:4), xlabel('Confidence level score'), ylabel('Frequency'), title('Confidence ratings distribution for correct responses');
% subplot(2,1,2),hist(All_False), ylim([0 60]), xticks(1:1:4), xlabel('Confidence level score'), ylabel('Frequency'), title('Confidence ratings distribution for incorrect responses');

%% confidence for hypo 2
%Congruent Trial, discrimination between IP and CP

for f = 1:length(Congruent_CP)
    if Congruent_CP(f,7) == 1
        Trial_Accuracy7(f,:) = 1;
    else
        Trial_Accuracy7(f,:) = -1;
    end
end
Confidence_CCP = Congruent_CP(:,8).*Trial_Accuracy7;

for g = 1:length(Congruent_IP)
    if Congruent_IP(g,7) == -1
        Trial_Accuracy8(g,:) = 1;
    else
        Trial_Accuracy8(g,:) = -1;
    end
end
Confidence_CIP = Congruent_IP(:,8).*Trial_Accuracy8;

%incongruent trial, discrimination between IP and CP
for h = 1:length(Incongruent_IP)
    if Incongruent_IP(h,7) == 1
        Trial_Accuracy9(h,:) = 1;
    else
        Trial_Accuracy9(h,:) = -1;
    end
end
Confidence_IIP = Incongruent_IP(:,8).*Trial_Accuracy9;

for i = 1:length(Incongruent_CP)
    if Incongruent_CP(i,7) == -1
        Trial_Accuracy10(i,:) = 1;
    else
        Trial_Accuracy10(i,:) = -1;
    end
end
Confidence_ICP = Incongruent_CP(:,8).*Trial_Accuracy10;

subplot(2,2,1),histogram(Confidence_CCP,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for congruent object in congruent trials')
subplot(2,2,2),histogram(Confidence_CIP,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for incongruent object in congruent trials')
subplot(2,2,3),histogram(Confidence_IIP,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for congruent object in incongruent trials')
subplot(2,2,4),histogram(Confidence_ICP,'BinWidth',1),xlim([-4 4]),xticks(-4:1:4), xlabel('Confidence level score'),ylabel('Frequency'), title('Confidence and decision accuracy combined for incongruent object in incongruent trials')


%% simulated data analysis
load('pilotresult_cop_02_01.mat')
TrialResponse = [CurrentTrialNumber QuestionNumber ImageFeature PatchSource PatchFeature Patch_location Response Confidence ResponseTime];
sum(TrialResponse(:,3))/length(TrialResponse); %%0.4851; not exactly 50 50 cuz this is calculated based on the number of patches were tested for congruent images
sum(TrialResponse(:,4)==1)/length(TrialResponse); %%0.7345; about 75% of nishimoto patches

subplot(3,2,1),histogram(TrialResponse(TrialResponse(:,4)==1,6)),xlabel('Location'), ylabel('Frequency'),title('Nishimoto patch location distribution');
subplot(3,2,2),histogram(TrialResponse(TrialResponse(:,4)==1,7)),xlabel('Confidence Scores'), ylabel('Frequency'),title('Nishimoto patch confidence ratings');
subplot(3,2,3),histogram(TrialResponse(TrialResponse(:,3)==0,8)),xlabel('Confidence Scores'),ylabel('Frequency'),title('Congruent Trial Confidence'); %change to 0 to see congruent
subplot(3,2,4),histogram(TrialResponse(TrialResponse(:,3)==1,8)),xlabel('Confidence Scores'),ylabel('Frequency'),title('incongruent Trial Confidence')
subplot(3,2,5),histogram(TrialResponse(TrialResponse(:,3)==0,7)),xlabel('Response'),ylabel('Frequency'),title('Congruent Trial Responses');
subplot(3,2,6),histogram(TrialResponse(TrialResponse(:,3)==1,7)),xlabel('Response'),ylabel('Frequency'),title('Incongruent Trial Responses');




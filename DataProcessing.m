 clear all            
% %% filepaths of patches and data files
% 
%     folderE = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\results and analysis\exp results";
%     filePatternE = fullfile(folderE, '*.mat');
%     theFilesE = struct2cell(dir(filePatternE));
%     selectnameE(1:60,1) = string(theFilesE(1,:));
%  
%     folderT = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\results and analysis\trial lists";
%     filePatternT = fullfile(folderT, '*.mat');
%     theFilesT = struct2cell(dir(filePatternT));
%     selectnameT(1:60,1) = string(theFilesT(1,:));
% 
%     folderI = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\incongruent patch";
%     filePatternI = fullfile(folderI, '*.jpg');
%     theFilesI = struct2cell(dir(filePatternI));
%     selectnameI(1:1260,1) = string(theFilesI(1,:));
% 
%     folderC = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\congruent patch";
%     filePatternC = fullfile(folderC, '*.jpg');
%     theFilesC = struct2cell(dir(filePatternC));
%     selectnameC(1:1260,1) = string(theFilesC(1,:));
% 
%     folderS = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\nishimoto patch';
%     filePatternS = fullfile(folderS, '*.jpg');
%     theFilesS = struct2cell(dir(filePatternS));
%     selectnameS(1:7044,1) = string(theFilesS(1,:));
% 
%     testPatch_dir = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\patches';
%     
% %% load data, and make a big matrix for each participant
% 
%    
%     ExpRst_1 = importdata(selectnameE((sub-1).*4 + 1,1));
%     ExpRst_1 = ExpRst_1(ExpRst_1(:,1)~= 1 & ExpRst_1(:,1) ~= 2,:); 
%     ExpRst_2 = importdata(selectnameE((sub-1).*4 + 2,1));
%     ExpRst_3 = importdata(selectnameE((sub-1).*4 + 3,1));
%     ExpRst_3 = ExpRst_3(ExpRst_3(:,1)~= 1 & ExpRst_3(:,1) ~= 2,:);
%     ExpRst_4 = importdata(selectnameE((sub-1).*4 + 4,1));
% 
% disp(selectnameE((sub-1).*4 + 1,1))
% disp(selectnameE((sub-1).*4 + 3,1))
% 
% 
%     TrialList_1 = importdata(selectnameT((sub-1).*4 + 1,1));
%     TrialList_2 = importdata(selectnameT((sub-1).*4 + 2,1));
%     TrialList_3 = importdata(selectnameT((sub-1).*4 + 3,1));
%     TrialList_4 = importdata(selectnameT((sub-1).*4 + 4,1));
% 
% disp(selectnameT((sub-1).*4 + 1,1))
% 
%     ExpRst = [ExpRst_1; ExpRst_2; ExpRst_3; ExpRst_4];
%     TrialList = [TrialList_1; TrialList_2; TrialList_3; TrialList_4];
% 
%     Image_num = zeros(length(TrialList),1);
%     for i = 1: length(TrialList)
%         Image_num(i,1)= str2num(TrialList(i,2));
%     end
% 
%     ExpRst = [Image_num ExpRst(:,2:9)];
%             
%% calculate patch difference for each patch
% ExpRst = importdata('Pooled Results.mat');
%      PatchDiff = zeros(length(ExpRst),1);
% for sub = 1:15     
%     
%      for a = 1: length(ExpRst)
%      
%         Ori_name = ['_' num2str(ExpRst(a,1)) '_' num2str(ExpRst(a,6)) '_'];
%          
%         if contains(char(TrialList(a,4)), '.jpg')
%             TestPatch_name = char(fullfile(testPatch_dir,TrialList(a,4)));
%             TestPatch = imread(TestPatch_name);
%             
%         else
%             
%             patchname = char(selectnameS(str2num(TrialList(a,4)),1));
%             patchname = fullfile(folderS, patchname);
%             TestPatch = imread(patchname);
%             
%         end
%         
%        
%          if ExpRst(a,3)== 0
%              Ori_fullname = selectnameC(contains(selectnameC,Ori_name));
%              Ori_fullname = char(fullfile(folderC,Ori_fullname));
%              Ori_patch = imread(Ori_fullname);
%          else
%              Ori_fullname = selectnameI(contains(selectnameI, Ori_name));
%              Ori_fullname = char(fullfile(folderI, Ori_fullname));
%              Ori_patch = imread(Ori_fullname);
%          end
%        
%        TestPatch = imresize(TestPatch, [147 147]);
%        Ori_patch = imresize(Ori_patch, [147 147]);
%        PatchDiff(a,1) =  sum(sum(sum(TestPatch - Ori_patch)));
%          
%      end
%      
% %      ExpRst2 = [ExpRst PatchDiff];
% %      
% %      if sub == 1
% %          ExpRst_pooled = ExpRst2;
% %      else
% %          ExpRst_pooled = [ExpRst_pooled; ExpRst2];
% %      end
% %      
% %      
% %     if sub <10
% %         save_path = ['Expresult_0' num2str(sub) '.mat'];
% %     else
% %         save_path  = ['Expresult_' num2str(sub) '.mat'];
% %     end
% 
%     %save(save_path, 'ExpRst2','-mat');   
% 
% end

% for sub = 1:15
%     
%     if sub < 10
%         RstName = ['Expresult_0' num2str(sub) '.mat'];
%     else
%         RstName = ['Expresult_' num2str(sub) '.mat'];
%     end
%    
%     ExpRst = importdata(RstName);
%     
%     subnumber = zeros(length(ExpRst), 1) + sub;
%     
%     ExpRst = [subnumber ExpRst];
%     
%     save(RstName,'ExpRst','-mat');
%     
%      if sub == 1
%          ExpRst_pooled = ExpRst;
%      else
%          ExpRst_pooled = [ExpRst_pooled; ExpRst];
%      end
% end
% 
% % 
% ExpRst_pooled = importdata('Pooled & Cleaned Results.mat');
% ResponseTime = ExpRst_pooled(:,10);
% histogram(ResponseTime);
% quantileRT = quantile(ResponseTime, [0 0.025 0.5 0.975 1]);
% histogram(ResponseTime(ResponseTime <= quantileRT(4) & ResponseTime >= quantileRT(2))),xlabel('Response Time by Seconds'), ylabel('Frequency'), title('Distribution of Response Time');
% % boxplot(ResponseTime(ResponseTime <= quantileRT(4) & ResponseTime >= quantileRT(2)));
% 
% PooledRst_cleaned = ExpRst_pooled(ExpRst_pooled(:,10) <= quantileRT(4) & ExpRst_pooled(:,10) >= quantileRT(2),:);
% 
% %save('Pooled Results.mat', 'ExpRst_pooled', '-mat');
% save('Pooled & Cleaned Results.mat', 'PooledRst_cleaned', '-mat');

%% data cleaning

% results = importdata('Pooled Results.mat');
% results = results(:,2:12);
% 
%     folderT = "C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\results and analysis\trial lists";
%     filePatternT = fullfile(folderT, '*.mat');
%     theFilesT = struct2cell(dir(filePatternT));
%     selectnameT(1:60,1) = string(theFilesT(1,:));
% 
% for i = 1:60
%     TrialList_sub = importdata(selectnameT(i, 1));
%     if i == 1
%         TrialList = TrialList_sub(:,4);
%     else
%         TrialList = [TrialList; TrialList_sub(:,4)];
%     end
% end
% 
% TrialList_pooled = zeros(length(TrialList),1);
% 
% for a = 1:length(TrialList)
%     if contains(TrialList(a,:),'.jpg')
%         TrialList(a,:)= "0";
%         TrialList_pooled(a,:)= str2num(TrialList(a,:));
%     else
%         TrialList_pooled(a,:)= str2num(TrialList(a,:));
%     end
%    
% end
% 
% results = [results TrialList_pooled];
% 
% save('Pooled Results.mat','results','-mat');
% 
% LogRT = log10(results(:,10));
% results_cleaned = results(LogRT <= 0.63 & LogRT >= -0.70,:);
% 
% % histogram(LogRT), xlim([-1 1.5]), xlabel('Log-transformed response time'), ylabel('Frequency'), title('Distribution of log-transformed response time');
% % histogram(LogRT(LogRT <= 0.63 & LogRT >= -0.70)), xlabel('Log-transformed response time'), ylabel('Frequency'), title('Distribution of cleaned, log-transformed response time');
% % histogram(results(:,10)), xlim([0 50]),xlabel('Log-transformed response time'), ylabel('Frequency'), title('Distribution of log-transformed response time');
% 
% save('Pooled & Cleaned Results.mat', 'results_cleaned', '-mat');



%% data cleaning by accuracy: find participants that perform significantly worse than the others
clear all;
Results_pooled = importdata('Pooled & Cleaned Results.mat');
AUC_table = [1:1:15; zeros(2,15)];
Results_pooled(:,9) = Results_pooled(:,8).*Results_pooled(:,9);

for sub = 1:15

Results = Results_pooled(Results_pooled(:,1) == sub,:);

% signal present and absent trial classification

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

%% hypo 1, not by congruency
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);
% type 1 ROC curve analysis

Confidence_N = Results_N(:,9);
Confidence_AP = Results_AP(:,9);

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    
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
AUC_table(2,sub) = AUC;

% subplot(1,2,1),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC curve for A/P patch discrimination AUC = ', num2str(AUC)]);
% axis square

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

Confidence_N = Results_N(:,9);
Confidence_AP = Results_AP(:,9);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4));
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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]),2);
AUC_table(3,sub) = AUC; 
% 
% subplot(1,2,2),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC curve for critical objects discrimination AUC = ', num2str(AUC)]);
% axis square


Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];
end
% 
% disp(AUC_table);
% subplot(2,1,1), histogram(AUC_table(2,:),15);
% subplot(2,1,2), histogram(AUC_table(3,:),15);

%% calculate patch difference for each patch
ExpRst = importdata('Pooled Results.mat');
ExpRst = ExpRst(:,2:end);
PatchDiff = zeros(length(ExpRst),1);

folderT = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\results and analysis\trial lists";
filePatternT = fullfile(folderT, '*.mat');
theFilesT = struct2cell(dir(filePatternT));
selectnameT(1:60,1) = string(theFilesT(1,:));

folderI = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\incongruent patch";
filePatternI = fullfile(folderI, '*.jpg');
theFilesI = struct2cell(dir(filePatternI));
selectnameI(1:1260,1) = string(theFilesI(1,:));

folderC = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\congruent patch";
filePatternC = fullfile(folderC, '*.jpg');
theFilesC = struct2cell(dir(filePatternC));
selectnameC(1:1260,1) = string(theFilesC(1,:));

folderS = 'D:\New photos\nishimoto patch';
filePatternS = fullfile(folderS, '*.jpg');
theFilesS = struct2cell(dir(filePatternS));
selectnameS(1:7044,1) = string(theFilesS(1,:));

testPatch_dir = 'C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\patches';


for sub = 1:15   
    TrialList_1 = importdata(selectnameT((sub-1).*4 + 1,1));
    TrialList_2 = importdata(selectnameT((sub-1).*4 + 2,1));
    TrialList_3 = importdata(selectnameT((sub-1).*4 + 3,1));
    TrialList_4 = importdata(selectnameT((sub-1).*4 + 4,1));
    TrialList_indv = [TrialList_1; TrialList_2; TrialList_3; TrialList_4];
    if sub == 1
    TrialList = TrialList_indv;
    else
    TrialList = [TrialList; TrialList_indv];
    end
end

 for a = 1:length(ExpRst)
     disp(a);

    Ori_name = ['_' num2str(ExpRst(a,2)) '_' num2str(ExpRst(a,7)) '_'];

    if ExpRst(a,5)~= 1
        patchname = char(fullfile(testPatch_dir,TrialList(a,4)));
        TestPatch = imread(patchname);
    else
        patchname = char(selectnameS(str2num(TrialList(a,4)),1));
        patchname = fullfile(folderS, patchname);
        TestPatch = imread(patchname);

    end

     if ExpRst(a,4) == 0
         Ori_fullname = selectnameC(contains(selectnameC,Ori_name));
         Ori_fullname = char(fullfile(folderC,Ori_fullname));
         Ori_patch = imread(Ori_fullname);
     else
         Ori_fullname = selectnameI(contains(selectnameI, Ori_name));
         Ori_fullname = char(fullfile(folderI, Ori_fullname));
         Ori_patch = imread(Ori_fullname);
     end

    TestPatch = double(imresize(TestPatch, [147 147]));
    Ori_patch = double(imresize(Ori_patch, [147 147]));
    PatchDiff(a,1) =  sum(sum(sum(TestPatch - Ori_patch)));

    
 end
 
 ExpRst(:,11)= PatchDiff;
 save('Pooled Trial List.mat','TrialList','-mat');
 
 TrialList_pooled = zeros(length(TrialList),1);
 
 for a = 1:length(TrialList)
    if contains(TrialList(a,4),'.jpg')
        TrialList_pooled(a,:)= 0;
    else
        TrialList_pooled(a,:)= str2num(TrialList(a,4));
    end
   
 end

ExpRst = [ExpRst TrialList_pooled];

save('Pooled Results.mat','ExpRst','-mat')


% 
ExpRst_P = importdata('Pooled Results.mat');
ExpRst_p = ExpRst_P(ExpRst_P(:,10)>=0.4 & ExpRst_P(:,10)<=4,:);
%LogRT = log10(ExpRst(:,10));
%ExpRst = [ExpRst LogRT];

for condition = 1:2
    ExpRst = ExpRst_p(ExpRst_p(:,4)== condition-1, :);
for sub = 1:15
    Results_indv = ExpRst(ExpRst(:,1)==sub,:);
    mean_indv = mean(Results_indv(:,10));
    sd_indv = 3.*std(Results_indv(:,10));
    min = mean_indv - sd_indv;
    max = mean_indv + sd_indv;
    Results_cleaned_indv = Results_indv(Results_indv(:,10)<= max & Results_indv(:,10)>= min,:); 
    
    disp(length(Results_cleaned_indv)/length(Results_indv));
    
    if sub == 1
        results_cleaned = Results_cleaned_indv;
    else
        results_cleaned = [results_cleaned; Results_cleaned_indv];
    end
    
end
if condition == 1
    results_final_cleaned = results_cleaned;
else
    results_final_cleaned = [results_final_cleaned; results_cleaned];
end
end
(length(ExpRst_P)-length(results_final_cleaned))/length(ExpRst_P);
%  m_pooled = mean(results_cleaned(:,end));
%  sd_pooled = 3.*std(results_cleaned(:,end));
%  pmin = m_pooled- sd_pooled;
%  pmax = m_pooled + sd_pooled;
%  results_cleaned_final = results_cleaned(results_cleaned(:,end)<= pmax & results_cleaned(:,end)>= pmin,:);

subplot(1,2,1),histogram(LogRT),xlim([-1 1.5]), xlabel('Log-transformed response time'), ylabel('Frequency'), title('Distribution of log-transformed response time');
axis square
subplot(1,2,2), histogram(results_cleaned_final(:,end)),xlabel('Log-transformed response time'), ylabel('Frequency'), title('Distribution of log-transformed response time after cleaning');
axis square

results_cleaned_final = results_cleaned_final(:,1:12);

size(ExpRst_P,1)-size(results_cleaned_final)

% save('Pooled & Cleaned Results.mat','results_cleaned_final','-mat');

results = importdata('Pooled & Cleaned Results.mat');
max(results(:,10));

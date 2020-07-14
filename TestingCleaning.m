ExpRst_P = importdata('Pooled Results.mat');
ExpRst_p = ExpRst_P(ExpRst_P(:,10)>=0.4 & ExpRst_P(:,10)<=4,:);
%LogRT = log10(ExpRst(:,10));
%ExpRst = [ExpRst LogRT];

% this filters reaction time

for condition = 1:2
    ExpRst = ExpRst_p(ExpRst_p(:,4)== condition-1, :);
for sub = 1:15
    Results_indv = ExpRst(ExpRst(:,1)==sub,:);
    mean_indv = mean(Results_indv(:,10));
    sd_indv = 3.*std(Results_indv(:,10));
    min = mean_indv - sd_indv;
    max = mean_indv + sd_indv;
    Results_cleaned_indv = Results_indv(Results_indv(:,10)<= max & Results_indv(:,10)>= min,:); 
    
  
    
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
(size(ExpRst_P,1)-size(results_final_cleaned,1))/length(ExpRst_P)

%% hypo 1, not by congruency
Results = results_final_cleaned;
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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2)

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
Results_N = Results(Find_Incongruent_CP,:);
Results_AP = Results(Find_Incongruent_IP,:);

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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2)

clear Confidence_N;
clear Confidence_AP;
clear Results_N;
clear Results_AP;
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC ;
clear results_final_cleaned;
clear Results;

%% method 2
disp('method 2');
ExpRst_P = importdata('Pooled Results.mat');
LogRT = log10(ExpRst_P(:,10));
ExpRst_P = [ExpRst_P LogRT];

for condition = 1:2
    ExpRst = ExpRst_p(ExpRst_p(:,4)== condition-1, :);
for sub = 1:15
    Results_indv = ExpRst(ExpRst(:,1)==sub,:);
    mean_indv = mean(Results_indv(:,end));
    sd_indv = 3.*std(Results_indv(:,end));
    min = mean_indv - sd_indv;
    max = mean_indv + sd_indv;
    Results_cleaned_indv = Results_indv(Results_indv(:,end)<= max & Results_indv(:,end)>= min,:); 
    
    
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
(length(ExpRst_P)-length(results_final_cleaned))/length(ExpRst_P)

Results = results_final_cleaned;
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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2)

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
Results_N = Results(Find_Incongruent_CP,:);
Results_AP = Results(Find_Incongruent_IP,:);

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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2)

clear Confidence_N;
clear Confidence_AP;
clear Results_N;
clear Results_AP;
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC ;
clear results_final_cleaned;
clear Results;


%% method 3
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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2)

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
Results_N = Results(Find_Incongruent_CP,:);
Results_AP = Results(Find_Incongruent_IP,:);

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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2)

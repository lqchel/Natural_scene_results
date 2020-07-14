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

%% for hypo one

Results_NC = Results(Find_N,:);
Results_APC = Results(Find_IAP|Find_CAP,:);
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
    matrix1(sub,:)= AUC;



end

%% for hypo 2
grandmatrix = zeros(15,21);

    
    for sub = 1:15

    Confidence_P = Results((Find_Congruent_CP | Find_Incongruent_IP) & Results(:,1)==sub,9);
    Confidence_A = Results((Find_Congruent_IP | Find_Incongruent_CP)& Results(:,1)==sub,9);

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

grandmatrix(sub)= AUC;

clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
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

boxplot(matrix1)
set(gca,'FontSize',12,'XTickLabel','Present v.s. absent patches','FontSize',12);
ylabel('AUC'),ylim([0.5 1]);
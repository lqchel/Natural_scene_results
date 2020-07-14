% type 2 analysis for the experiment'

clear all
Results = importdata('Pooled Results.mat');
%% location
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];
location = [0 6.5 9.2];

%% trial index
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

%% raw confidence distribution

for sub = 1:15
    m(sub) = mean(Results(Results(:,1)==sub,9));
    sd(sub) = std(Results(Results(:,1)==sub,9));
end
subplot(1,2,1), histogram(m),title('Subject mean confidence'),ylim([0 8]);
subplot(1,2,2), histogram(sd),title('Individual confidence SD'), ylim([0 8]);

%% z-transform confidence
conf_array = zeros(length(Results),1);
for sub = 1:15
    Confidence = Results(Results(:,1)==sub,9);
    for s = 1:length(Confidence)
        z_conf(s) = (Confidence(s) - m(sub))/sd(sub);
    end
    if sub == 1
    conf_array(1:length(Confidence),:)= z_conf;
    row_id = length(Confidence);
    else
    conf_array(row_id+1:row_id+length(Confidence),:) = z_conf;
    row_id = row_id+length(Confidence);
    end
    clear Confidence
end
conf_array = [Results(:,1) conf_array];
Results = [Results conf_array(:,2)];

%% confidence for original patches, hit and miss
correct = Results_AP(Results_AP(:,8)==1,:);
incorrect = Results_AP(Results_AP(:,8)==-1,:);
for loc = 1:3
    for sub = 1:15
        hit(sub,loc) = mean(correct(correct(:,1)==sub & correct(:,end) == location(loc),9));
        miss(sub,loc) = mean(incorrect(incorrect(:,1)==sub & incorrect(:,end) == location(loc),9));
    end
end
se_hit = within_se(hit,15,3);
se_miss = within_se(miss,15,3);

errorbar(location,mean(hit,1),se_hit);
hold on
errorbar(location,mean(miss,1),se_miss);
hold off

table1 = table(correct(:,9),correct(:,end),correct(:,1),'VariableNames',{'Confidence','Eccentricity','Subjects'});
lm_hit = fitlme(table1, 'Confidence ~ Eccentricity + (1|Subjects)');
table2 = table(incorrect(:,9),incorrect(:,end),incorrect(:,1),'VariableNames',{'Confidence','Eccentricity','Subjects'});
lm_miss = fitlme(table2, 'Confidence ~ Eccentricity + (1|Subjects)');

%% hypo 1, separating congruent/incongruent conditions
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

%% hypo 1, not separating conditions
Results = Results(Results(:,1)==7,:);
Results_NC = Results(Find_N,:); % all trials with absent test probes
Results_APC = Results(Find_IAP|Find_CAP,:); % all trials with present test probes
Results_Correct = [Results_NC(Results_NC(:,8)==-1,:); Results_APC(Results_APC(:,8)==1,:)];
Results_Incorrect = [Results_NC(Results_NC(:,8)==1,:); Results_APC(Results_APC(:,8)==-1,:)];

matrix1 = zeros(15,1,3);
hit_matrix = zeros(15,5,3);
FA_matrix = zeros(15,5,3);

for a = 1:3
    for sub = 1:15 

    Confidence_Incorrect = Results_Incorrect(Results_Incorrect(:,1)==sub & Results_Incorrect(:,13)==location(a),9);
    Confidence_Correct = Results_Correct(Results_Correct(:,1)==sub& Results_Correct(:,13)==location(a),9);
%     Confidence_Incorrect = Results_Incorrect(:,9);
%     Confidence_Correct = Results_Correct(:,9);
    
    for i = -4:-1
        Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
        Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
    end
    for i = 1:4
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i) = Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit];
    Cumulative_FA = [0 Cumulative_FA];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

    matrix1(sub,:,a)= AUC; % create matrix for individual AUCs
    hit_matrix(sub,:,a) = Cumulative_Hit;
    FA_matrix(sub,:,a) = Cumulative_FA;
    
    clear Cumulative_Hit
    clear Cumulative_FA
    clear AUC
 
    end
end
final_AUC = mean(mean(mean(matrix1)));
final_AUC = round(final_AUC,2);
colours = cbrewer('qual','Pastel2',3);
fill([mean(mean(FA_matrix,3)) 1],[mean(mean(hit_matrix,3)) 0],colours(3,:));
hold on
plot(mean(mean(FA_matrix,3)),mean(mean(hit_matrix,3)),'bo-','LineWidth',1.2);
plot([0 1],[0 1], '--k');
hold off
axis square
title(['Type 2 ROC for present vs N, AUC = ' num2str(final_AUC)]);

% across eccentricity
AUC_ecc = reshape(mean(matrix1,1),[1,3]);
plot([0 6.48 9.16], AUC_ecc,'bo-','LineWidth',1.2);
xticks([0 6.48 9.16]),ylim([0.5 0.8]);

%% hypo 2
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

%% type 2 of participant 7
Results = Results(Results(:,1)==7,:);
Results_NC = Results(Find_N,:); % all trials with absent test probes
Results_APC = Results(Find_IAP|Find_CAP,:); % all trials with present test probes
Results_Correct = [Results_NC(Results_NC(:,8)==-1,:); Results_APC(Results_APC(:,8)==1,:)];
Results_Incorrect = [Results_NC(Results_NC(:,8)==1,:); Results_APC(Results_APC(:,8)==-1,:)];

Confidence_Correct = Results_Correct(:,9);
Confidence_Incorrect = Results_Incorrect(:,9);

    for i = -4:-1
        Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
        Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
    end
    for i = 1:4
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i) = Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit];
    Cumulative_FA = [0 Cumulative_FA];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);



colours = cbrewer('qual','Pastel2',3);
fill([1 Cumulative_FA],[0 Cumulative_Hit],colours(3,:));
hold on
plot(Cumulative_FA,Cumulative_Hit,'bo-','LineWidth',1.2);
plot([0 1],[0 1], 'k--');

hold off
axis square
title('b)');
set(gca,'FontName','Arial','FontSize',12)
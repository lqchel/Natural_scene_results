%% manuscript first draft analysis
Results = importdata('OnlinePilotData.mat');

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];

location = [0 6.5 9.2];

colours = cbrewer('qual','Set1',8);

Results(:,9) = Results(:,8).*Results(:,9);
Find_N = Results(:,5) ==1;

% % present patch trials -- signal present for hypo 1
% Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
% Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 
% 
% %Congruent trial with Congruent object, and incongruent trial with
% %incongruent object -- signal present for hypo 2
% Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
% Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;
% 
% %Incongruent trial with congruent object, congruent trial with incongruent
% %object -- signal absent for hypo 2
% Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
% Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1;
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2

Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;


%% hypothesis 1 AUC

Results_NC = Results(Find_N,:); % find absent patchese
Results_APC = Results(Find_IAP|Find_CAP,:); % find present patches
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

se1 = std(matrix1)/sqrt(15);

% AUC across eccentricities

matrix3 = zeros(15, 3);
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_N,:);
    indvP = Results(Results(:,1)==sub & (Find_CAP|Find_IAP),:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i); % AUC calculation
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
    matrix3(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

se2 = within_se(matrix3,15,3);

% plot graph

subplot(2,2,1),errorbar(-5,mean(matrix1),se1,'.','MarkerSize',12,...
    'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
ylabel('Objective Type 1 AUC');
xlabel('Eccentricity (dva)');
xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar([0 6.5 9.2],mean(matrix3),se2,'.-','MarkerSize',12,...
    'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
title('present vs. null','FontName','Arial');

grandmatrix = zeros(15,2);

%% hypothesis 2 AUC
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
    Confidence_APCounts(i+5) = nansum(Confidence_P == -i);
    Confidence_NCounts(i+5) = nansum(Confidence_A == -i);
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

% across eccentricities

%congruent

matrix4 = zeros(15, 3);
location = [0 6.5 9.2];
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Congruent_IP,:);
    indvP = Results(Results(:,1)==sub & Find_Congruent_CP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = nansum(Confidence_AP == -i); % AUC calculation
        Confidence_NCounts(i+5) = nansum(Confidence_N == -i);
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
    matrix4(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC
    
end
end



%incongruent condition

clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff

matrix5 = zeros(15, 3);

for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Incongruent_CP,:);
    indvP = Results(Results(:,1)==sub & Find_Incongruent_IP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = nansum(Confidence_AP == -i); % AUC calculation
        Confidence_NCounts(i+5) = nansum(Confidence_N == -i);
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
    matrix5(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

% t tests

[h1,p1,CI,stats] = ttest(matrix4(:,1),0.5,'alpha',0.05/3)
[h2,p2,CI,stats] = ttest(matrix4(:,2),0.5,'alpha',0.05/3)
[h3,p3,CI,stats] = ttest(matrix4(:,3),0.5,'alpha',0.05/3)

[h4,p4,CI,stats] = ttest(matrix5(:,1),0.5,'alpha',0.05/3)
[h5,p5,CI,stats] = ttest(matrix5(:,2),0.5,'alpha',0.05/3)
[h6,p6,CI,stats] = ttest(matrix5(:,3),0.5,'alpha',0.05/3)


subplot(2,2,2),d = errorbar(-5,mean(grandmatrix(:,1)),std(grandmatrix(:,1))/sqrt(15),'d','MarkerSize',6,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
ylabel('Objective Type 1 AUC');
xlabel('Eccentricity (dva)');
xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
hold on
errorbar(-5,mean(grandmatrix(:,2)),std(grandmatrix(:,2))/sqrt(15),'d','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.5 9.2],nanmean(matrix4),within_se(matrix4,15,3),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.5 9.2],nanmean(matrix5),within_se(matrix5,15,3),'d--','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
plot([-7 11],[0.5 0.5],'k--');
title('original vs. modified','FontName','Arial');
legend({'Congruent','Incongruent'},'Box','off');
hold off

% end


%% hypothesis 2 hit and FA
    matrix1 = zeros(15,3,2);
    location = [0 6.5 9.2];

    for loc = 1:3
    for sub = 1:15
        R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        end

        matrix1(sub,loc,condition) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end

    end

    %incongruent --------------------------------------------------------------------------------
     matrix2 = zeros(15,3,2);

    for loc = 1:3
    for sub = 1:15
        R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP|Find_Incongruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        end

        matrix2(sub,loc,condition) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end

    end

%% for data table
%congruent
for condition = 1:2
    for loc = 1:3
        current_accuracy = matrix1(:,loc,condition);
        patch_type = zeros(15,1) + condition;
        eccentricity = zeros(15,1) + location(loc);
        con_loc_matrix = [num current_accuracy eccentricity patch_type];
        
        if loc == 1
            condition_matrix = con_loc_matrix;
        else
            condition_matrix = [condition_matrix; con_loc_matrix];
        end
        
        clear current_accuracy
        clear patch_type
        clear eccentricity
        clear con_loc_matrix
    end
    
    if condition == 1
        all_matrix_cong = condition_matrix;
    else
        all_matrix_cong = [all_matrix_cong; condition_matrix];
    end
    
end
%incongruent
for condition = 1:2
    for loc = 1:3
        current_accuracy = matrix2(:,loc,condition);
        patch_type = zeros(15,1) + condition;
        eccentricity = zeros(15,1) + location(loc);
        con_loc_matrix = [num current_accuracy eccentricity patch_type];
        
        if loc == 1
            condition_matrix = con_loc_matrix;
        else
            condition_matrix = [condition_matrix; con_loc_matrix];
        end
        
        clear current_accuracy
        clear patch_type
        clear eccentricity
        clear con_loc_matrix
    end
    
    if condition == 1
        all_matrix_incong = condition_matrix;
    else
        all_matrix_incong = [all_matrix_cong; condition_matrix];
    end
    
end
 data4 = table(all_matrix_cong(:,1),all_matrix_cong(:,2),all_matrix_cong(:,3),all_matrix_cong(:,4),'VariableNames',{'Subject','Accuracy',...
     'Eccentricity','PatchType'});
lme5 = fitlme(data4,'Accuracy ~ Eccentricity*PatchType + (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject)' );
lme4 = fitlme(data4,'Accuracy ~ Eccentricity + PatchType+ (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject) ');
lme3 = fitlme(data4,'Accuracy ~ Eccentricity + Eccentricity:PatchType+ (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject) ');
lme1 = fitlme(data4,'Accuracy ~ PatchType + Eccentricity:PatchType+ (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject) ');
congruent_interaction = compare(lme4,lme5)
congruent_eccentricity = compare(lme1,lme5)
congruent_patch = compare(lme3,lme5)
 
 
 
data5 = table(all_matrix_incong(:,1),all_matrix_incong(:,2),all_matrix_incong(:,3),all_matrix_incong(:,4),'VariableNames',{'Subject','Accuracy',...
     'Eccentricity','PatchType'});
 
lme6 = fitlme(data5,'Accuracy ~ Eccentricity*PatchType + (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject)' );
lme7 = fitlme(data5,'Accuracy ~ Eccentricity + PatchType+ (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject) ');
lme8 = fitlme(data5,'Accuracy ~ Eccentricity + Eccentricity:PatchType+ (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject) ');
lme9 = fitlme(data5,'Accuracy ~ PatchType + Eccentricity:PatchType+ (1|Subject) + (PatchType|Subject)+ (Eccentricity|Subject) ');
incongruent_interaction = compare(lme7,lme6)
incongruent_eccentricity = compare(lme9,lme6)
incongruent_patch = compare(lme8,lme6)

%% hypothesis 1 type 2 AUC

Results(:,9) = abs(Results(:,9));
Results_NC = Results(Find_N,:); % all trials with absent test probes
Results_APC = Results(Find_IAP|Find_CAP,:); % all trials with present test probes
Results_Correct = [Results_NC(Results_NC(:,8)==-1,:); Results_APC(Results_APC(:,8)==1,:)];
Results_Incorrect = [Results_NC(Results_NC(:,8)==1,:); Results_APC(Results_APC(:,8)==-1,:)];

matrix6 = zeros(15,3);


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

    matrix6(sub,a)= AUC; % create matrix for individual AUCs
%     hit_matrix(sub,:,a) = Cumulative_Hit;
%     FA_matrix(sub,:,a) = Cumulative_FA;
    
    clear Cumulative_Hit
    clear Cumulative_FA
    clear AUC
 
    end
end



% final_AUC = round(final_AUC,2);
se7 = std(mean(matrix6,2))/sqrt(15);

% AUC on each eccentricity levels
AUC_ecc = reshape(mean(matrix6,1),[1,3]);
se4 = within_se(matrix6,15,3);

% t test
[h,p,CI,stats] = ttest(mean(matrix6,2),0.5)
[h,p,CI,stats] = ttest(matrix6(:,1),0.5,'alpha',0.05/3)
[h,p,CI,stats] = ttest(matrix6(:,2),0.5,'alpha',0.05/3)
[h,p,CI,stats] = ttest(matrix6(:,3),0.5,'alpha',0.05/3)

% lme model

for i = 1:15
    num(i,1) = i;
end
sub_num = [num; num; num];
location_num = [zeros(15,1); zeros(15,1)+6.5; zeros(15,1)+ 9.2];
AUC_lme = [matrix6(:,1);matrix6(:,2);matrix6(:,3)];
data1 = table(sub_num,location_num,AUC_lme,'VariableNames',{'Subject','location','AUC'});

lm1 = fitlme(data1,'AUC~location + (1|Subject)');
lm0 = fitlme(data1,'AUC~ 1+ (1|Subject)');
compare(lm0,lm1)

subplot(2,2,3),errorbar(-5,mean2(matrix6),se7,'.','MarkerSize',14,...
    'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
ylabel('Subjective Type 2 AUC');


xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 0.8]);
xlabel('Eccentricity (dva)');

legend('off');
hold on
errorbar([0 6.5 9.2],AUC_ecc,se4,'.-','MarkerSize',14,...
    'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
title('present vs. null','FontName','Arial');

clear Confidence_Correct
clear Confidence_Incorrect

%% hypothesis 2 type 2 AUC 
matrix7 = zeros(15,2);


for condition = 1:2
    
    for sub = 1:15
    if condition ==1
        Results_N = Results(Find_Congruent_IP,:);
        Results_A = Results(Find_Congruent_CP,:);
        Confidence_Correct = [Results_N(Results_N(:,1)== sub & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== sub &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==-1,9)];
    else
        Results_N = Results(Find_Incongruent_CP,:);
        Results_A = Results(Find_Incongruent_IP,:);
        Confidence_Correct = [Results_N(Results_N(:,1)== sub &Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== sub &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==-1,9)];
    end

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

    matrix7(sub,condition)= AUC; % create matrix for individual AUCs

clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end




% across eccentricities

%congruent

matrix9 = zeros(15, 3);
location = [0 6.5 9.2];
for condition = 1:3
    
    for sub = 1:15
  
        Results_N = Results(Find_Congruent_IP & Results(:,end)== location(condition),:);
        Results_A = Results(Find_Congruent_CP & Results(:,end) == location(condition),:);
        Confidence_Correct = [Results_N(Results_N(:,1)== sub & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== sub &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==-1,9)];
 

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

    matrix9(sub,condition)= AUC; % create matrix for individual AUCs

clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end


%incongruent condition
clear Results_N
clear Results_A


matrix0 = zeros(15, 3);

for condition = 1:3
    
    for sub = 1:15
  
        Results_N = Results(Find_Incongruent_CP & Results(:,end)== location(condition),:);
        Results_A = Results(Find_Incongruent_IP & Results(:,end) == location(condition),:);
        Confidence_Correct = [Results_N(Results_N(:,1)== sub & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== sub &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==sub & Results_A(:,8)==-1,9)];
 

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

    matrix0(sub,condition)= AUC; % create matrix for individual AUCs

clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end

% t test
[h,p,CI,stats] = ttest(matrix7(:,1),0.5,'alpha',0.025)
[h,p,CI,stats] = ttest(matrix7(:,2),0.5,'alpha',0.025)

[h1,p1,CI,stats] = ttest(matrix9(:,1),0.5,'alpha',0.05/3)
[h2,p2,CI,stats] = ttest(matrix9(:,2),0.5,'alpha',0.05/3)
[h3,p3,CI,stats] = ttest(matrix9(:,3),0.5,'alpha',0.05/3)

[h,p,CI,stats] = ttest(matrix0(:,1),0.5,'alpha',0.05/3)
[h,p,CI,stats] = ttest(matrix0(:,2),0.5,'alpha',0.05/3)
[h,p,CI,stats] = ttest(matrix0(:,3),0.5,'alpha',0.05/3)

[h,p,CI,stats] = ttest(matrix7(:,1),matrix7(:,2),'alpha',0.05/3)

% for data table
clear matrix1
matrix1(:,:,1) = matrix9;
matrix1(:,:,2) = matrix0;
for condition = 1:2
    for loc = 1:3
        current_accuracy = matrix1(:,loc,condition);
        patch_type = zeros(15,1) + condition;
        eccentricity = zeros(15,1) + location(loc);
        con_loc_matrix = [num current_accuracy eccentricity patch_type];
        
        if loc == 1
            condition_matrix = con_loc_matrix;
        else
            condition_matrix = [condition_matrix; con_loc_matrix];
        end
        
        clear current_accuracy
        clear patch_type
        clear eccentricity
        clear con_loc_matrix
    end
    
    if condition == 1
        all_matrix = condition_matrix;
    else
        all_matrix = [all_matrix; condition_matrix];
    end
    
end

% lme models with condition x eccentricity interaction 
data6 = table(all_matrix(:,1),all_matrix(:,2),all_matrix(:,3),all_matrix(:,4),'VariableNames',{'Subject','AUC','Eccentricity','Condition'});
lme1 = fitlme(data6, 'AUC~Eccentricity*Condition + (1|Subject) + (Eccentricity|Subject) + (Condition|Subject)');
lme2 = fitlme(data6, 'AUC~Eccentricity + Condition + (1|Subject) + (Eccentricity|Subject) + (Condition|Subject)');
lme3 = fitlme(data6, 'AUC~Eccentricity:Condition + Condition + (1|Subject) + (Eccentricity|Subject) + (Condition|Subject)');
lme4 = fitlme(data6, 'AUC~Eccentricity:Condition + Eccentricity + (1|Subject) + (Eccentricity|Subject) + (Condition|Subject)');

ecc_effect = compare(lme3,lme1)
con_effect = compare(lme2,lme1)
interaction_effect = compare(lme4,lme1)

% lme models

for i = 1:15
    num(i,1) = i;
end
sub_num = [num; num; num];
location_num = [zeros(15,1); zeros(15,1)+6.5; zeros(15,1)+ 9.2];
%congruent
AUC_lme1 = [matrix9(:,1);matrix9(:,2);matrix9(:,3)];
data2 = table(sub_num,location_num,AUC_lme1,'VariableNames',{'Subject','location','AUC'});

lm2 = fitlme(data2,'AUC~location + (1|Subject)');
lm3 = fitlme(data2,'AUC~ 1+ (1|Subject)');
compare(lm3,lm2)

%incongruent
AUC_lme2 = [matrix0(:,1);matrix0(:,2);matrix0(:,3)];
data3 = table(sub_num,location_num,AUC_lme2,'VariableNames',{'Subject','location','AUC'});

lm4 = fitlme(data3,'AUC~location + (1|Subject)');
lm5 = fitlme(data3,'AUC~ 1+ (1|Subject)');
compare(lm5,lm4)


% 
% subplot(2,2,4),d = errorbar(-5,mean(matrix7(:,1)),std(matrix7(:,1))/sqrt(15),'d','MarkerSize',6,...
%     'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
% ylabel('Subjective Type 2 AUC');
% xlabel('Eccentricity (dva)');
% xlim([-7 11]),xticks([-5 0 6.5 9.2]);
% set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
% ylim([0.4 0.8]);
% legend('off');
% hold on
% errorbar(-5,mean(matrix7(:,2)),std(matrix7(:,2))/sqrt(15),'d','MarkerSize',6,...
%     'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
% errorbar([0 6.5 9.2],nanmean(matrix9),within_se(matrix9,15,3),'d-','MarkerSize',6,...
%     'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
% errorbar([0 6.5 9.2],nanmean(matrix0),within_se(matrix0,15,3),'d--','MarkerSize',6,...
%     'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
% plot([-7 11],[0.5 0.5],'k--');
% title('original vs. modified','FontName','Arial');
% legend({'Congruent','Incongruent'},'Box','off');
% hold off

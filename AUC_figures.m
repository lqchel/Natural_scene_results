%%% this function plots Type 1 and Type 2 AUC results for Qianchen's experiment. data =
%%% data matrix, num_sub = number of subjects in the data matrix.
%%% for this function to work, please make sure to install PsychToolbox (http://psychtoolbox.org/)


function [out] = AUC_figures(data)
% Results = data(data(:,11)~=0,:); %% exclude catch trials
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results'));
colours = cbrewer('qual','Set1',8);
for exp = 1:2
    
Results = data{exp};
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));

% Results(:,9) = Results(:,8).*Results(:,9);
Find_N = Results(:,5) ==1;

location = [0 1 2]; % 0 = fovea, 1 = peripheral (2,4,6,8), 2 = para-peripheral (1,3,7,9)

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1; %% Results(:,6) == 1 for exp 2
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;
end

%% hypothesis 1 AUC

Results_NC = Results(Find_N,:); % find absent patches
Results_APC = Results(Find_IAP|Find_CAP,:); % find present patches
matrix1 = zeros(num_sub,1);

for sub = 1:num_sub 

    Results_NC = Results(Find_N,:);
    Confidence_N = Results_NC(Results_NC(:,1)==subject_id(sub),9);
    Confidence_AP = Results_APC(Results_APC(:,1)==subject_id(sub),9);
    
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

if num_sub >= 3
se1 = std(matrix1);
end

% AUC across eccentricities

matrix3 = zeros(num_sub, 3);
for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_N,:);
    indvP = Results(Results(:,1)==subject_id(sub) & (Find_CAP|Find_IAP),:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,13)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,13) == location(a),:);
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

if num_sub >= 3
se2 = std(matrix3);
end
% plot graph
out = figure;

if num_sub < 3
    subplot(2,2,1),plot(-5,mean(matrix1),'.','MarkerSize',14,...
    'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    hold on
    plot([0 6.5 9.2],mean(matrix3,1),'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    plot([-7 11],[0.5 0.5],'k--');
    hold off
    title('present vs. null','FontName','Arial');
else
    subplot(2,2,1),errorbar(-5,mean(matrix1),se1,'.','MarkerSize',14,...
        'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
    hold on
    errorbar([0 6.5 9.2],mean(matrix3),se2,'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
    plot([-7 11],[0.5 0.5],'k--');
    hold off
end  
    ylabel('Objective Type 1 AUC');
    xlabel('Patch location');
    xlim([-7 11]),xticks([-5 0 6.5 9.2]);
    set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
    ylim([0.4 1]);
    legend('off');
    title('present vs. null','FontName','Arial');




%% hypothesis 2 AUC
grandmatrix = zeros(num_sub,2);
for condition = 1:2
    
    for sub = 1:num_sub
    if condition ==1
        Confidence_P = Results(Find_Congruent_CP & Results(:,1)==subject_id(sub),9);
        Confidence_A = Results(Find_Congruent_IP & Results(:,1)==subject_id(sub),9);
    else
        Confidence_P = Results(Find_Incongruent_IP & Results(:,1)==subject_id(sub),9);
        Confidence_A = Results(Find_Incongruent_CP & Results(:,1)==subject_id(sub),9);
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

matrix4 = zeros(num_sub, 3);
for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_Congruent_IP,:);
    indvP = Results(Results(:,1)==subject_id(sub) & Find_Congruent_CP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,13)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,13) == location(a),:);
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

matrix5 = zeros(num_sub, 3);

for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_Incongruent_CP,:);
    indvP = Results(Results(:,1)==subject_id(sub) & Find_Incongruent_IP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,13)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,13) == location(a),:);
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

figure(out);

if num_sub<3
    subplot(2,2,2),plot(-5,nanmean(grandmatrix(:,1)),'d','MarkerSize',6,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    legend('off');
    hold on
    plot(-5,mean(grandmatrix(:,2)),'d','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([0 6.5 9.2],nanmean(matrix4,1),'-d','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([0 6.5 9.2],nanmean(matrix5,1),'--d','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([-7 11],[0.5 0.5],'k--');
    legend({'Congruent','Incongruent'},'Box','off');
    hold off
else
    subplot(2,2,2),errorbar(-5,mean(grandmatrix(:,1)),std(grandmatrix(:,1))/sqrt(num_sub),'d','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    legend('off');
    hold on
    errorbar(-5,mean(grandmatrix(:,2)),std(grandmatrix(:,2))/sqrt(num_sub),'d','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    errorbar([0 6.5 9.2],nanmean(matrix4),within_se(matrix4,num_sub,3),'d-','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    errorbar([0 6.5 9.2],nanmean(matrix5),within_se(matrix5,num_sub,3),'d--','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([-7 11],[0.5 0.5],'k--');
    legend({'Congruent','Incongruent'},'Box','off');
    hold off
end

ylabel('Objective Type 1 AUC');
xlabel('Patch location');
xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
ylim([0.4 1]);
title('original vs. modified','FontName','Arial');


%% hypothesis 1 type 2 AUC

Results(:,9) = abs(Results(:,9));
Results_NC = Results(Find_N,:); % all trials with absent test probes
Results_APC = Results(Find_IAP|Find_CAP,:); % all trials with present test probes
Results_Correct = [Results_NC(Results_NC(:,8)==-1,:); Results_APC(Results_APC(:,8)==1,:)];
Results_Incorrect = [Results_NC(Results_NC(:,8)==1,:); Results_APC(Results_APC(:,8)==-1,:)];

matrix6 = zeros(num_sub,3);


for a = 1:3
    for sub = 1:num_sub 

    Confidence_Incorrect = Results_Incorrect(Results_Incorrect(:,1)==subject_id(sub) & Results_Incorrect(:,13)==location(a),9);
    Confidence_Correct = Results_Correct(Results_Correct(:,1)==subject_id(sub)& Results_Correct(:,13)==location(a),9);
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


% AUC on each eccentricity levels
AUC_ecc = reshape(nanmean(matrix6,1),[1,3]);

if num_sub >= 3
se4 = within_se(matrix6,num_sub,3);
se7 = nanstd(mean(matrix6,2))/sqrt(num_sub);
end

figure(out);
if num_sub<3
    
    subplot(2,2,3),plot(-5,nanmean(nanmean(matrix6)),'.','MarkerSize',14,...
    'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    hold on
    plot([0 6.5 9.2],AUC_ecc,'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    plot([-7 11],[0.5 0.5],'k--');
    hold off

else
    subplot(2,2,3),errorbar(-5,nanmean(nanmean(matrix6)),se7,'.','MarkerSize',14,...
        'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
    hold on
    errorbar([0 6.5 9.2],AUC_ecc,se4,'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
    plot([-7 11],[0.5 0.5],'k--');
    hold off

end

ylabel('Subjective Type 2 AUC');
xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
ylim([0.4 0.9]);
xlabel('Patch location');
title('present vs. null','FontName','Arial');
legend('off');
clear Confidence_Correct
clear Confidence_Incorrect

%% hypothesis 2 type 2 AUC 
matrix7 = zeros(num_sub,2);


for condition = 1:2
    
    for sub = 1:num_sub
    if condition ==1
        Results_N = Results(Find_Congruent_IP,:);
        Results_A = Results(Find_Congruent_CP,:);
        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
    else
        Results_N = Results(Find_Incongruent_CP,:);
        Results_A = Results(Find_Incongruent_IP,:);
        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
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

matrix9 = zeros(num_sub, 3);
for condition = 1:3
    
    for sub = 1:num_sub
  
        Results_N = Results(Find_Congruent_IP & Results(:,13)== location(condition),:);
        Results_A = Results(Find_Congruent_CP & Results(:,13) == location(condition),:);
        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
 
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


matrix0 = zeros(num_sub, 3);

for condition = 1:3
    
    for sub = 1:num_sub
  
        Results_N = Results(Find_Incongruent_CP & Results(:,13)== location(condition),:);
        Results_A = Results(Find_Incongruent_IP & Results(:,13) == location(condition),:);
        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
 

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

figure(out);
disp(matrix0)

if num_sub < 3
    subplot(2,2,4),plot(-5,mean(matrix7(:,1)),'d','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    legend('off');
    hold on
    plot(-5,mean(matrix7(:,2)),'d','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([0 6.5 9.2],nanmean(matrix9,1),'-d','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([0 6.5 9.2],nanmean(matrix0,1),'--d','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([-7 11],[0.5 0.5],'k--');
    legend({'Congruent','Incongruent'},'Box','off');
    hold off
else
    subplot(2,2,4),errorbar(-5,mean(matrix7(:,1)),std(matrix7(:,1))/sqrt(num_sub),'d','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    legend('off');
    hold on
    errorbar(-5,mean(matrix7(:,2)),std(matrix7(:,2))/sqrt(15),'d','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    errorbar([0 6.5 9.2],nanmean(matrix9),within_se(matrix9,num_sub,3),'d-','MarkerSize',6,...
        'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    errorbar([0 6.5 9.2],nanmean(matrix0),within_se(matrix0,num_sub,3),'d--','MarkerSize',6,...
        'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
    plot([-7 11],[0.5 0.5],'k--');
    legend({'Congruent','Incongruent'},'Box','off');
    hold off
end

title('original vs. modified','FontName','Arial');
ylabel('Subjective Type 2 AUC');
xlabel('Patch location');
xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
ylim([0.4 0.9]);
end

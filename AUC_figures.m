%%% this function plots Type 1 and Type 2 AUC results for Qianchen's experiment. data =
%%% data matrix, num_sub = number of subjects in the data matrix.
%%% for this function to work, please make sure to install PsychToolbox (http://psychtoolbox.org/)


function [out] = AUC_figures(data)
% Results = data(data(:,11)~=0,:); %% exclude catch trials
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results'));
colour_palatte = cbrewer('qual','Set1',8);
colours = [colour_palatte(2,:); colour_palatte(5,:)];
for exp = 1:2
Results = data{exp};
Find_N{exp} = Results(:,5) ==1;
location = [0 1 2]; % 0 = fovea, 1 = peripheral (2,4,6,8), 2 = para-peripheral (1,3,7,9)
% present patch trials -- signal present for hypo 1
Find_CAP{exp} = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP{exp} = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP{exp} = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1; %% Results(:,6) == 1 for exp 2
Find_Incongruent_IP{exp} = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP{exp} = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
Find_Incongruent_CP{exp} = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;

clear Results
end

%% hypothesis 1 AUC
out = figure;
subplot(2,3,1);
for exp = 1:2
Results = data{exp};
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));

Results_NC = Results(Find_N{exp},:); % find absent patches
Results_APC = Results(Find_IAP{exp}|Find_CAP{exp},:); % find present patches
matrix1 = zeros(num_sub,1);

for sub = 1:num_sub 

    Results_NC = Results(Find_N{exp},:);
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
se1 = std(matrix1)/sqrt(num_sub);
end

% AUC across eccentricities

matrix3 = zeros(num_sub, 3);
for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_N{exp},:);
    indvP = Results(Results(:,1)==subject_id(sub) & (Find_CAP{exp}|Find_IAP{exp}),:); % trial classification
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
se2 = std(matrix3)/sqrt(num_sub);
end
% plot graph
errorbar(-5,mean(matrix1),se1,'.','MarkerSize',14,...
'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
if exp == 1
    hold on
end
errorbar([0 6.5 9.2],mean(matrix3),se2,'.-','MarkerSize',14,...
'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');

if exp == 2
hold off
ylabel('Objective Type 1 AUC');
xlim([-7 11]),xticks([-5 0 6.5 9.2]);
set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
ylim([0.4 1]);
legend('off');
title('P+O vs. N','FontName','Arial');
end
clear matrix1
clear matrix3

end

%% hypothesis 1 type 2 AUC
figure(out);
subplot(2,3,4);

for exp = 1:2
Results = data{exp};
Results(:,9) = abs(Results(:,9));
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));
Results_NC = Results(Find_N{exp},:); % all trials with absent test probes
Results_APC = Results(Find_IAP{exp}|Find_CAP{exp},:); % all trials with present test probes
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
    
    clear Cumulative_Hit
    clear Cumulative_FA
    clear AUC
 
    end
end

% AUC on each eccentricity levels
AUC_ecc = reshape(nanmean(matrix6,1),[1,3]);

if num_sub >= 3
se4 = nanstd(matrix6)/sqrt(num_sub);
se7 = nanstd(mean(matrix6,2))/sqrt(num_sub);
end


errorbar(-5,nanmean(nanmean(matrix6)),se7,'.','MarkerSize',14,...
    'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
if exp == 1
    hold on
end
errorbar([0 6.5 9.2],AUC_ecc,se4,'.-','MarkerSize',14,...
    'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');

if exp == 2
    hold off
    ylabel('Subjective Type 2 AUC');
    xlim([-7 11]),xticks([-5 0 6.5 9.2]);
    set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
    ylim([0.4 0.9]);
    xlabel('Patch location');
    legend('off');
end
clear Confidence_Correct
clear Confidence_Incorrect

end

%% hypothesis 2, Congruent Type 1 AUC
figure(out);
subplot(2,3,2);

for exp = 1:2
Results = data{exp};
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));    
grandmatrix = zeros(num_sub,1);
    
    for sub = 1:num_sub
        Confidence_P = Results(Find_Congruent_CP{exp} & Results(:,1)==subject_id(sub),9);
        Confidence_A = Results(Find_Congruent_IP{exp} & Results(:,1)==subject_id(sub),9);

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

        grandmatrix(sub)= AUC;

        clear Confidence_P
        clear Confidence_A
        clear Cumulative_NCounts;
        clear Cumulative_APCounts;
        clear Cumulative_Hit;
        clear Cumulative_FA;
        clear AUC;
    end

% across eccentricities

%congruent

    matrix4 = zeros(num_sub, 3);
    for sub = 1:num_sub
        indvN = Results(Results(:,1)==subject_id(sub) & Find_Congruent_IP{exp},:);
        indvP = Results(Results(:,1)==subject_id(sub) & Find_Congruent_CP{exp},:); % trial classification
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

    errorbar(-5,mean(grandmatrix),nanstd(grandmatrix)/sqrt(num_sub),'.','MarkerSize',14,...
            'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1);
    if exp == 1
        hold on
    end
    errorbar([0 6.5 9.2],nanmean(matrix4),nanstd(matrix4)/sqrt(num_sub),'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1);
    if exp == 2
        plot([-7 11],[0.5 0.5],'k--');
        hold off
        legend('off');
        xlim([-7 11]),xticks([-5 0 6.5 9.2]);
        set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
        ylim([0.4 1]);
        title('O vs. M, Congruent','FontName','Arial');
    end
    clear grandmatrix
    clear matrix4
end

%% Hypothesis 2, Type 1, Incongruent

clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff

figure(out);
subplot(2,3,3);

for exp = 1:2
Results = data{exp};
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));    
grandmatrix = zeros(num_sub,1);
matrix5 = zeros(num_sub, 3);

% collapsed across eccentricities
    for sub = 1:num_sub
        Confidence_P = Results(Find_Incongruent_IP{exp} & Results(:,1)==subject_id(sub),9);
        Confidence_A = Results(Find_Incongruent_CP{exp} & Results(:,1)==subject_id(sub),9);

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

        grandmatrix(sub)= AUC;

        clear Confidence_P
        clear Confidence_A
        clear Cumulative_NCounts;
        clear Cumulative_APCounts;
        clear Cumulative_Hit;
        clear Cumulative_FA;
        clear AUC;
    end


% on each eccentricity level
        for sub = 1:num_sub
            indvN = Results(Results(:,1)==subject_id(sub) & Find_Incongruent_CP{exp},:);
            indvP = Results(Results(:,1)==subject_id(sub) & Find_Incongruent_IP{exp},:); % trial classification
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
        clear AUC
        end
        end
       
       errorbar(-5,mean(grandmatrix),nanstd(grandmatrix)/sqrt(num_sub),'.','MarkerSize',14,...
                'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1);
        if exp == 1
            hold on
        end
        errorbar([0 6.5 9.2],nanmean(matrix5),nanstd(matrix5)/sqrt(num_sub),'.-','MarkerSize',14,...
            'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1);
        if exp == 2
            plot([-7 11],[0.5 0.5],'k--');
            hold off
            xlim([-7 11]),xticks([-5 0 6.5 9.2]);
            set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
            ylim([0.4 1]);
            title('O vs. M, Incongruent','FontName','Arial');
        end
        legend('','Experiment 1','','Experiment 2','','Box','off');
        clear grandmatrix
        clear matrix5
end

%% Hypothesis 2 Type 2 AUC, Congruent
figure(out);
subplot(2,3,5);
for exp = 1:2
Results = data{exp};
Results(:,9) = abs(Results(:,9));
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));
matrix7 = zeros(num_sub,1);
    
    for sub = 1:num_sub

        Results_N = Results(Find_Congruent_IP{exp},:);
        Results_A = Results(Find_Congruent_CP{exp},:);
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

    matrix7(sub)= AUC; % create matrix for individual AUCs

clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

% across eccentricities

%congruent

matrix9 = zeros(num_sub, 3);
for condition = 1:3
    
    for sub = 1:num_sub
  
        Results_N = Results(Find_Congruent_IP{exp} & Results(:,13)== location(condition),:);
        Results_A = Results(Find_Congruent_CP{exp} & Results(:,13) == location(condition),:);
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

    errorbar(-5,nanmean(matrix7),nanstd(matrix7)/sqrt(num_sub),'.','MarkerSize',14,...
        'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
    if exp == 1
        hold on
    end
    errorbar([0 6.5 9.2],nanmean(matrix9),nanstd(matrix9)/sqrt(num_sub),'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
    plot([-7 11],[0.5 0.5],'k--');

    if exp == 2
        hold off
        xlim([-7 11]),xticks([-5 0 6.5 9.2]);
        set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
        ylim([0.4 0.9]);
        xlabel('Patch location');
        legend('off');
    end
    clear Confidence_Correct
    clear Confidence_Incorrect
    clear matrix7
    clear matrix9
end

%% Hypothesis 2, Type 2, Incongruent
clear Results_N
clear Results_A
figure(out);
subplot(2,3,6);

for exp = 1:2
Results = data{exp};
Results(:,9) = abs(Results(:,9));
num_sub = length(unique(Results(:,1)));
subject_id = unique(Results(:,1));
matrix7 = zeros(num_sub,1);
matrix9 = zeros(num_sub, 3);

   for sub = 1:num_sub

        Results_N = Results(Find_Incongruent_CP{exp},:);
        Results_A = Results(Find_Incongruent_IP{exp},:);
        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
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

    matrix7(sub)= AUC; % create matrix for individual AUCs

clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

for condition = 1:3
    
    for sub = 1:num_sub
  
        Results_N = Results(Find_Incongruent_CP{exp} & Results(:,13)== location(condition),:);
        Results_A = Results(Find_Incongruent_IP{exp} & Results(:,13) == location(condition),:);
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
    errorbar(-5,nanmean(matrix7),nanstd(matrix7)/sqrt(num_sub),'.','MarkerSize',14,...
        'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
    if exp == 1
        hold on
    end
    errorbar([0 6.5 9.2],nanmean(matrix9),nanstd(matrix9)/sqrt(num_sub),'.-','MarkerSize',14,...
        'MarkerFaceColor',colours(exp,:),'MarkerEdgeColor',colours(exp,:),'Color',colours(exp,:),'LineWidth',1.2,'Capsize',10);
    plot([-7 11],[0.5 0.5],'k--');

    if exp == 2
        hold off
        xlim([-7 11]),xticks([-5 0 6.5 9.2]);
        set(gca,'XTickLabel',{'All','F','P','P-P'},'FontSize',12,'Box','off');
        ylim([0.4 0.9]);
        xlabel('Patch location');
        legend('off');
    end
    clear matrix7
    clear matrix9
end
end

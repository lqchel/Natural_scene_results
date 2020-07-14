
% function [out] = massive_analysis(data,num_sub)

Results = importdata('Pooled Results 2.mat'); %uncomment for lab testing
%results
% Results = data;
num_sub = 15;
Results(:,9) = Results(:,8).*Results(:,9);

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
Find_N = Results(:,5) ==1;

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

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];

%% hypothesis 1

colours = cbrewer('qual', 'Set1', 8); 
grandmatrix = zeros(2,2);
matrix1 = zeros(num_sub,2);

for sub = 1:num_sub
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
 
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
    end
    matrix1(sub,condition) = accuracy;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end
% [h,p,ci]= ttest(matrix1(:,1),matrix1(:,2));

% calculate standard errors
grandmatrix(1,:) = mean(matrix1);

for condition = 1:2
grandmatrix(2,condition) = std(matrix1(:,condition))/sqrt(num_sub);
end

%% hypothesis 1 with eccentricity


colours = cbrewer('qual', 'Set1', 8); 
grandmatrix1 = zeros(2,3,2);

% hit and CR rate
matrix1 = zeros(num_sub,3,2);
location = [0 6.48 9.16];

for loc = 1:3
for sub = 1:num_sub
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N)& Results(:,end)==location(loc),:);

for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
    end
 
    matrix1(sub,loc,condition) = accuracy;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

end

for condition = 1:2 
    
    grandmatrix1(1,:,condition) = mean(matrix1(:,:,condition));
    grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
    
end

out = figure;
subplot(2,4,1),errorbar(-5,grandmatrix(1,1),grandmatrix(2,1),'d','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8);
hold on
errorbar(-5,grandmatrix(1,2),grandmatrix(2,2),'o','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
    'k','Color','k','LineWidth',0.8);
hold off
ylabel('%Accurate judgments');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
legend('off');

hold on
errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
    'k','Color','k','LineWidth',0.8,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
legend({'Hit','CR'});
title('a)','FontName','Arial');
clear grandmatrix1
clear matrix1
clear R_indv
clear grandmatrix
clear Find_patch




%% hypothesis 1 AUC

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];

Results_NC = Results(Find_N,:); % find absent patchese
Results_APC = Results(Find_IAP|Find_CAP,:); % find present patches
matrix1 = zeros(num_sub,1);

for sub = 1:num_sub 
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

se1 = std(AUC)/sqrt(num_sub);

% AUC across eccentricities

matrix3 = zeros(num_sub, 3);
for sub = 1:num_sub
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

se2 = within_se(matrix3,num_sub,3);

% plot graph
figure(out)
subplot(2,4,2),errorbar(-5,mean(matrix1),se1,'d','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
ylabel('AUC');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar([0 6.48 9.16],mean(matrix3),se2,'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
title('b)','FontName','Arial');



%% hypothesis 2
clear grandmatrix
clear matrix1

% hypothesis 2 eccentricity

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];
    % congruent ---------------------------------------------------------------------------------------------------
    matrix1 = zeros(num_sub,3,2);
    location = [0 6.48 9.16];

    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
        end

        matrix1(sub,loc,condition) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end

    end

    for condition = 1:2 

        grandmatrix1(1,:,condition) = nanmean(matrix1(:,:,condition));
        grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
        
    end
    disp(grandmatrix1)
    disp(matrix1)


    %incongruent --------------------------------------------------------------------------------
     matrix2 = zeros(num_sub,3,2);

    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP|Find_Incongruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
        end

        matrix2(sub,loc,condition) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end

    end

    for condition = 1:2 

        grandmatrix2(1,:,condition) = nanmean(matrix2(:,:,condition));
        grandmatrix2(2,:,condition)= within_se(matrix2(:,:,condition),num_sub,3);

    end

    % calculate hit and CR collapsed across eccentricity levels
    all_cong(1,:) = [nanmean(grandmatrix1(1,:,1)) nanmean(grandmatrix1(1,:,2))]; 
    all_incong(1,:) = [nanmean(grandmatrix2(1,:,1)) nanmean(grandmatrix2(1,:,2))];

    for c = 1:2
        all_cong(2,c) = std(mean(matrix1(:,:,c),2))/sqrt(num_sub);
        all_incong(2,c) = std(mean(matrix2(:,:,c),2))/sqrt(num_sub);
    end
    
    % plot figures
    % congruent
    
figure(out);

for a = 1:2
    ax(2.*a - 1)= a - 0.14;
    ax(2.*a) = a + 0.14;
end
    m_hypo1 = [all_cong(1,:) all_incong(1,:)];
    se_hypo1 = [all_cong(2,:) all_incong(2,:)];

subplot(2,4,3), s = bar([all_cong(1,:);all_incong(1,:)],'grouped');
s(1).FaceColor = colours(2,:);
s(2).FaceColor = colours(1,:);

hold on
errorbar(ax,m_hypo1,se_hypo1,'k.','MarkerSize',1,'LineWidth',1);
hold off
ylim([0 1]);
ylabel('%Accurate judgments');
set(gca,'XTickLabel',{'Congruent','Incongruent'},'FontSize',12);
title('c)','FontName','Arial');


% errorbar(-5,all_cong(1,1),all_cong(2,1),'d','MarkerSize',8,...
%     'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
% hold on
% errorbar(-5, all_cong(1,2), all_cong(2,2),'o','MarkerSize',6,...
%     'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(1,:),'LineWidth',1);
% errorbar(-5,all_incong(1,1),all_incong(2,1),'d','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
%     colours(2,:),'Color',colours(2,:),'LineWidth',1);
% errorbar(-5,all_incong(1,2),all_incong(2,2),'o','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
%     colours(1,:),'Color',colours(1,:),'LineWidth',1);
% hold off
% ylabel('%Accurate judgments');
% xlim([-7 11]),xticks([-5 0 6.48 9.16]);
% set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
% ylim([0 1]);
% legend('off');


% hold on
% errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',6,...
%     'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',6,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',...
%     colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix2(1,:,1),grandmatrix2(2,:,1),'d--','MarkerSize',7,...
%     'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix2(1,:,2),grandmatrix2(2,:,2),'o--','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
%     colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
% plot([-7 11],[0.5 0.5],'k--');
% hold off
% legend({'Congruent Hit','Congruent CR','Incongruent Hit','Incongruent CR'});




%% hypothesis 2 AUC
grandmatrix = zeros(num_sub,2);


for condition = 1:2
    
    for sub = 1:num_sub
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

matrix4 = zeros(num_sub, 3);
location = [0 6.48 9.16];
for sub = 1:num_sub
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

matrix5 = zeros(num_sub, 3);

for sub = 1:num_sub
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

figure(out);
subplot(2,4,4),d = errorbar(-5,mean(grandmatrix(:,1)),std(grandmatrix(:,1))/sqrt(num_sub),'d','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
ylabel('AUC');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar(-5,mean(grandmatrix(:,2)),std(grandmatrix(:,2))/sqrt(num_sub),'d','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.48 9.16],nanmean(matrix4),within_se(matrix4,num_sub,3),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.48 9.16],nanmean(matrix5),within_se(matrix5,num_sub,3),'d--','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
plot([-7 11],[0.5 0.5],'k--');
legend({'Congruent','Incongruent'});
title('d)','FontName','Arial');

hold off

% end



%% hypothesis 1 confidence
% 
grandmatrix = zeros(2,2);
matrix1 = zeros(num_sub,2);

for sub = 1:num_sub
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
    R_indv(:,9) = R_indv(:,9).* R_indv(:,8);
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1 & R_indv(:,8) ==1;
        Results_P = R_indv(Find_patch,:);
        confidence = sum(Results_P(:,9))/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1 & R_indv(:,8) == -1;
        Results_P = R_indv(Find_patch,:);
        confidence = sum(Results_P(:,9))/size(Results_P,1);
    end
    
    matrix1(sub,condition) = confidence;
    
    clear Results_P
    clear condition_mean 
    clear confidence
    clear Find_patch
end
end

% calculate standard errors
grandmatrix(1,:) = mean(matrix1);
grandmatrix(2,:) = within_se(matrix1,num_sub,2);



% across eccetricities
grandmatrix1 = zeros(2,3,2);
matrix1 = zeros(num_sub,3,2);
location = [0 6.48 9.16];

for loc = 1:3
for sub = 1:num_sub
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N)& Results(:,end)==location(loc),:);
    R_indv(:,9) = R_indv(:,9).* R_indv(:,8);
for condition = 1:2
   if condition == 1
        Find_patch = R_indv(:,5)~=1 & R_indv(:,8) ==1;
        Results_P = R_indv(Find_patch,:);
        confidence = sum(Results_P(:,9))/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1 & R_indv(:,8) == -1;
        Results_P = R_indv(Find_patch,:);
        confidence = sum(Results_P(:,9))/size(Results_P,1);
    end
     
    matrix1(sub,loc,condition) = confidence;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

end

for condition = 1:2 
    
    grandmatrix1(1,:,condition) = mean(matrix1(:,:,condition));
    grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
    
end

figure(out);
subplot(2,4,5),errorbar(-5,grandmatrix(1,1),grandmatrix(2,1),'d','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8);
hold on
errorbar(-5,grandmatrix(1,2),grandmatrix(2,2),'o','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
    'k','Color','k','LineWidth',0.8);
hold off
ylabel('Confidence level');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([1 4]),yticks([1:1:4]);
legend('off');

hold on
errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
    'k','Color','k','LineWidth',0.8,'Capsize',10);
hold off
legend({'Hit','CR'});
title('e)','FontName','Arial');

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

final_AUC = mean(mean(matrix6));
final_AUC = round(final_AUC,2);
se3 = std(mean(mean(matrix6)))/sqrt(num_sub);

% AUC on each eccentricity levels
AUC_ecc = reshape(mean(matrix6,1),[1,3]);
se4 = within_se(matrix6,num_sub,3);
figure(out)

subplot(2,4,6),errorbar(-5,final_AUC,se3,'d','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
ylabel('AUC');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar([0 6.48 9.16],AUC_ecc,se4,'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
title('f)','FontName','Arial');

clear Confidence_Correct
clear Confidence_Incorrect

%% hypothesis 2 confidence

Results(:,9) = abs(Results(:,9));

    % congruent ---------------------------------------------------------------------------------------------------
    matrix7 = zeros(num_sub,3,2);
    location = [0 6.48 9.16];

    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==2 & R_indv(:,8)==1;
            Results_P = R_indv(Find_patch,:);
            confidence = sum(Results_P(:,9))/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==3 & R_indv(:,8)==-1;
            Results_P = R_indv(Find_patch,:);
            confidence = sum(Results_P(:,9))/size(Results_P,1);
        end

        matrix7(sub,loc,condition) = confidence;

        clear Results_P
        clear condition_mean 
        clear confidence
    end
    end

    end

    for condition = 1:2 

        grandmatrix7(1,:,condition) = nanmean(matrix7(:,:,condition));
        grandmatrix7(2,:,condition)= within_se(matrix7(:,:,condition),num_sub,3);
        
    end
    


    %incongruent --------------------------------------------------------------------------------
     matrix2 = zeros(num_sub,3,2);

    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP|Find_Incongruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5) == 3 & R_indv(:,8) == 1;
            Results_P = R_indv(Find_patch,:);
            confidence = sum(Results_P(:,9))/size(Results_P,1);
        else
            Find_patch = R_indv(:,5) == 2 & R_indv(:,8)==-1;
            Results_P = R_indv(Find_patch,:);
            confidence = sum(Results_P(:,9))/size(Results_P,1);
        end

        matrix2(sub,loc,condition) = confidence;

        clear Results_P
        clear condition_mean 
        clear confidence
    end
    end

    end

    for condition = 1:2 

        grandmatrix2(1,:,condition) = nanmean(matrix2(:,:,condition));
        grandmatrix2(2,:,condition)= within_se(matrix2(:,:,condition),num_sub,3);

    end

    % calculate hit and CR collapsed across eccentricity levels
    all_cong(1,:) = [mean(grandmatrix7(1,:,1)) mean(grandmatrix7(1,:,2))];
    all_incong(1,:) = [mean(grandmatrix2(1,:,1)) mean(grandmatrix2(1,:,2))];
    
        for c = 1:2
        all_cong(2,c) = std(mean(matrix1(:,:,c),2))/sqrt(num_sub);
        all_incong(2,c) = std(mean(matrix2(:,:,c),2))/sqrt(num_sub);
        end
    
    % plot figures
    % congruent
    
figure(out);
subplot(2,4,7),errorbar(-5,all_cong(1,1),all_cong(2,1),'d','MarkerSize',8,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
hold on
errorbar(-5, all_cong(1,2), all_cong(2,2),'o','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(1,:),'LineWidth',1);
errorbar(-5,all_incong(1,1),all_incong(2,1),'d','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
    colours(2,:),'Color',colours(2,:),'LineWidth',1);
errorbar(-5,all_incong(1,2),all_incong(2,2),'o','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
    colours(1,:),'Color',colours(1,:),'LineWidth',1);
hold off
ylabel('Confidence level');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([1 4]);
legend('off');


hold on
errorbar([0 6.48 9.16],grandmatrix7(1,:,1),grandmatrix7(2,:,1),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix7(1,:,2),grandmatrix7(2,:,2),'o-','MarkerSize',6,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',...
    colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix2(1,:,1),grandmatrix2(2,:,1),'d--','MarkerSize',7,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix2(1,:,2),grandmatrix2(2,:,2),'o--','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
    colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
legend({'Congruent Hit','Congruent CR','Incongruent Hit','Incongruent CR'});
title('g)','FontName','Arial');

%% hypothesis 2 type 2 AUC 
matrix7 = zeros(num_sub,2);


for condition = 1:2
    
    for sub = 1:num_sub
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

matrix9 = zeros(num_sub, 3);
location = [0 6.48 9.16];
for condition = 1:3
    
    for sub = 1:num_sub
  
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


matrix0 = zeros(num_sub, 3);

for condition = 1:3
    
    for sub = 1:num_sub
  
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



figure(out);
subplot(2,4,8),d = errorbar(-5,mean(matrix7(:,1)),std(matrix7(:,1))/sqrt(num_sub),'d','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
ylabel('AUC');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.5','9.2'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar(-5,mean(matrix7(:,2)),std(matrix7(:,2))/sqrt(num_sub),'d','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.48 9.16],nanmean(matrix9),within_se(matrix9,num_sub,3),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.48 9.16],nanmean(matrix0),within_se(matrix0,num_sub,3),'d--','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
plot([-7 11],[0.5 0.5],'k--');
legend({'Congruent','Incongruent'});
title('h)','FontName','Arial');
hold off

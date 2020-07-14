 %% final output...
clear all;
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

colours = cbrewer('qual', 'Set1', 8); 
colour = cbrewer('qual','Set2',8);

%% percentage hypo 2

  
for i = 1:9

    for sub = 1:15
        percentage_CP(sub,i) = size(Results(Results(:,1)== sub & Find_Congruent_CP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Congruent_CP,:));
        percentage_CI(sub,i) = size(Results(Results(:,1)== sub &Find_Congruent_IP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Congruent_IP,:));
        percentage_IP(sub,i) =  size(Results(Results(:,1)== sub &Find_Incongruent_IP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Incongruent_IP,:));
        percentage_II(sub,i) =  size(Results(Results(:,1)== sub & Find_Incongruent_CP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Incongruent_CP,:));
    end
end

percentage_CP = [percentage_CP(:,1:4) percentage_CP(:,6:9)];
percentage_CI = [percentage_CI(:,1:4) percentage_CI(:,6:9)];
percentage_IP = [percentage_IP(:,1:4) percentage_IP(:,6:9)];
percentage_II = [percentage_II(:,1:4) percentage_II(:,6:9)];


for a = 1:8
    ax(2.*a - 1)= a - 0.14;
    ax(2.*a) = a + 0.14;
end
    
for a = 1:8
    std_CPI(2.*a-1) = std(percentage_CP(:,a))/sqrt(15);
    std_CPI(2.*a) = std(percentage_CI(:,a))/sqrt(15);
    std_IPI(2.*a-1) = std(percentage_IP(:,a))/sqrt(15);
    std_IPI(2.*a) = std(percentage_II(:,a))/sqrt(15);
    pcg_CPI(2.*a-1) = mean(percentage_CP(:,a));
    pcg_CPI(2.*a) = mean(percentage_CI(:,a));
    pcg_IPI(2.*a -1) = mean(percentage_IP(:,a));
    pcg_IPI(2.*a) = mean(percentage_II(:,a));
end

percentage_CP = mean(percentage_CP,1);
percentage_CI = mean(percentage_CI,1);
percentage_IP = mean(percentage_IP,1);
percentage_II = mean(percentage_II,1);

percentage_CPI = [reshape(percentage_CP,[8,1]) reshape(percentage_CI,[8,1])];
percentage_IPI = [reshape(percentage_IP,[8,1]) reshape(percentage_II,[8,1])];


subplot(1,3,1),e = bar(percentage_CPI,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
e(1).FaceColor = colours(2,:);
e(2).FaceColor = colours(1,:);
xlabel('Decision x confidence'),
ylabel('Percentage');
hold on
errorbar(ax,pcg_CPI,std_CPI,'k.','MarkerSize',1,'LineWidth',0.8);
hold off
title('a)','FontWeight','normal');
legend({'Original','Different'},'Location','northwest');

subplot(1,3,2),f = bar(percentage_IPI,'grouped','BarWidth',1,'LineWidth',1.2);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
f(1).FaceColor = 'white';
f(2).FaceColor = 'white';
f(1).EdgeColor = colours(2,:);
f(2).EdgeColor = colours(1,:);
xlabel('Decision x confidence'),
ylim([0 0.5])
hold on
errorbar(ax,pcg_IPI,std_IPI,'k.','MarkerSize',1,'LineWidth',0.8);
hold off
title('b)','FontWeight','normal');
legend({'Original','Different'},'Location','northwest');



%hit and correct rejection for hypo 2

% grandmatrix = zeros(2,2,2);
% matrix1 = zeros(15,2,2);
% 
% for congruency = 1:2
%     if congruency == 1
%         Find_Original = Find_Congruent_CP;
%         Find_Different = Find_Congruent_IP;
%     else
%         Find_Original = Find_Incongruent_IP;
%         Find_Different = Find_Incongruent_CP;
%     end
%     
% for sub = 1:15
% 
% 
%     hit = sum(Results(:,1)==sub & Find_Original & Results(:,8)== 1)/size(Results(Results(:,1)==sub & Find_Original,:),1);
%     cr = sum(Results(:,1)==sub & Find_Different & Results(:,8)== -1)/size(Results(Results(:,1)==sub & Find_Different,:),1);
%     
%     matrix1(sub,1,congruency) = hit;
%     matrix1(sub,2,congruency) = cr;
%     
% clear hit
% clear cr
% 
% end
% end
% 
% 
% % calculate standard errors
% for condition = 1:2
% grandmean = repmat(mean(mean(matrix1(:,:,condition))),15,1);
% individualmean = mean(matrix1(:,:,condition),2);
% indi_diff = individualmean - grandmean;
% indi_diff = [indi_diff indi_diff];
% matrix1_corr(:,:,condition) = matrix1(:,:,condition) - indi_diff;
% end
% 
% for cond = 1:2
% for i = 1:2
%     grandmatrix(1,i,cond) = mean(matrix1_corr(:,i,cond));
%     grandmatrix(2,i,cond) = std(matrix1_corr(:,i,cond))/sqrt(15);
% end
% end
% 
% 
% subplot(2,2,3),d = bar(0.75,grandmatrix(1,1,1),'BarWidth',0.45);
% ylabel('Percentage of accurate judgment');
% ylim([0 1]);
% title('c)','FontWeight','normal');
% xlim([0 4]),xticks([1 3]),xticklabels({'Congruent','Incongruent'})
% set(gca,'FontSize',12);
% d.FaceColor = colours(2,:);
% hold on
% bar(1.25,grandmatrix(1,2,1),'BarWidth',0.45,'FaceColor',colours(1,:));
% bar(2.75,grandmatrix(1,1,2),'BarWidth',0.41,'FaceColor','white','EdgeColor',colours(2,:),'LineWidth',1.2);
% bar(3.25,grandmatrix(1,2,2),'BarWidth',0.41,'FaceColor','white','EdgeColor',colours(1,:),'LineWidth',1.2);
% hold off
% hold on
% errorbar([0.75 1.25 2.75 3.25],[grandmatrix(1,:,1) grandmatrix(1,:,2)],[grandmatrix(2,:,1) grandmatrix(2,:,2)],'k.',...
%     'MarkerSize',2,'LineWidth',1.2);
% hold off
% 
% [h,p,ci,stats]= ttest(matrix1(:,1,1),matrix1(:,2,1))
% [h,p,ci,stats]= ttest(matrix1(:,1,2),matrix1(:,2,2))
% [h,p,ci,stats]= ttest(matrix1(:,1,1),matrix1(:,1,2)) 
% [h,p,ci,stats]= ttest(matrix1(:,2,1),matrix1(:,2,2),0.0125); 


%AUC for hypo 2
matrix1 = zeros(15,2);

for condition = 1:2
    if condition == 1
    Results_NC = Results(Find_Congruent_IP,:);
    Results_APC = Results(Find_Congruent_CP,:);
    else
    Results_NC = Results(Find_Incongruent_CP,:);
    Results_APC = Results(Find_Incongruent_IP,:);
    end
    
    
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
    matrix1(sub,condition)= AUC;



end

end

population_m = mean(mean(mean(matrix1))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix1,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix1_corr = matrix1 - repmat(indv_dff, [1,2]);

%calculate within-subject sem on each location level
se = [0 0];
for i = 1:2
    se(i)= std(matrix1_corr(:,i))/sqrt(15);
end

subplot(1,3,3),c = bar(mean(matrix1,1));
c.FaceColor = colours(2,:);
hold on
errorbar(mean(matrix1,1),se,'k.', 'MarkerSize',2,'LineWidth',1);
ylim([0 1]);
hold off
set(gca,'XTickLabel',{'C','I'},'FontSize',12);
H = sigstar({{'C','I'}},0.01);
ylabel('AUC');
title('c)','FontWeight','normal');

%% for hypothesis 4 hits and correct rejection
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;

Results = [Results eccentricity];
location = [0 6.48 9.16];

grandmatrix = zeros(2,4,3);
matrix1 = zeros(15,4,3);
matrix2 = zeros(180,5);

for loc = 1:3
for congruency = 1:2
    if congruency == 1
        Find_Original = Find_Congruent_CP;
        Find_Different = Find_Congruent_IP;
    else
        Find_Original = Find_Incongruent_IP;
        Find_Different = Find_Incongruent_CP;
    end
    
for sub = 1:15


    hit = sum(Results(:,1)==sub & Find_Original & Results(:,8)== 1 & Results(:,end)==location(loc))/size(Results(Results(:,1)==sub & Find_Original...
        & Results(:,end)==location(loc),:),1);
    cr = sum(Results(:,1)==sub & Find_Different & Results(:,8)== -1&Results(:,end)==location(loc))/size(Results(Results(:,1)==sub & Find_Different...
        & Results(:,end)==location(loc),:),1);
    
    if congruency == 1
    matrix1(sub,1,loc) = hit;
    matrix1(sub,2,loc) = cr;
    else
    matrix1(sub,3,loc) = hit;
    matrix1(sub,4,loc) = cr;
    end
    if loc == 1 && congruency ==1 && sub == 1
        count = 1;
    else
        count = count + 2;
    end
    matrix2(count,1)= hit;
    matrix2(count,2)= sub;
    matrix2(count,3)= location(loc);
    matrix2(count,4)= 0;
    matrix2(count,5)= congruency;
    matrix2(count + 1, 1)= cr;
    matrix2(count+1,2)= sub;
    matrix2(count +1,3)= location(loc);
    matrix2(count+1,4)= 1;
    matrix2(count+1,5) = congruency;
    
clear hit
clear cr

end
end

end
matrix3 = matrix2(matrix2(:,4)==0,:);
matrix4 = matrix2(matrix2(:,4)==1,:);

data = table(matrix2(:,1),matrix2(:,2),matrix2(:,3),categorical(matrix2(:,4)),categorical(matrix2(:,5)),'VariableName',{'Accuracy','Subject','Eccentricity','PatchType','Congruency'});
lm1 = fitlme(data,'Accuracy ~ Eccentricity*PatchType*Congruency + (1|Subject)');

data2 = table(matrix3(:,1),matrix3(:,2),matrix3(:,3),categorical(matrix3(:,5)),'VariableName',{'Accuracy','Subject','Eccentricity','Congruency'});
lm2 = fitlme(data2,'Accuracy ~ Congruency*Eccentricity + (1|Subject)')
lm3 = fitlme(data2,'Accuracy ~ Congruency + Eccentricity + (1|Subject)');
lm4 = fitlme(data2, 'Accuracy~ Eccentricity + Congruency:Eccentricity + (1|Subject)');
lm5 = fitlme(data2,'Accuracy~ Congruency + Congruency:Eccentricity + (1|Subject)');
ctest_interaction = compare(lm3,lm2)
ctest_congruent = compare(lm4,lm2)
ctest_eccentricity = compare(lm5,lm2)

data3 = table(matrix4(:,1),matrix4(:,2),matrix4(:,3),categorical(matrix4(:,5)),'VariableName',{'Accuracy','Subject','Eccentricity','Congruency'});
lm6 = fitlme(data3,'Accuracy ~ Congruency*Eccentricity + (1|Subject)')
lm7 = fitlme(data3,'Accuracy ~ Congruency + Eccentricity + (1|Subject)');
lm8 = fitlme(data3, 'Accuracy~ Eccentricity + Congruency:Eccentricity + (1|Subject)');
lm9 = fitlme(data3,'Accuracy~ Congruency + Congruency:Eccentricity + (1|Subject)');
itest_interaction = compare(lm7,lm6)
itest_congruent = compare(lm8,lm6)
itest_eccentricity = compare(lm9,lm6)


% calculate standard errors
for loc = 1:3
for condition = 1:2
grandmean = repmat(mean(mean(matrix1(:,(2.*condition -1):(2.*condition),loc))),15,1);
individualmean = mean(matrix1(:,(2.*condition -1):(2.*condition),loc),2);
indi_diff = individualmean - grandmean;
indi_diff = [indi_diff indi_diff];
matrix1_corr(:,(2.*condition -1):(2.*condition),loc) = matrix1(:,(2.*condition -1):(2.*condition),loc) - indi_diff;
end
end

for loc = 1:3
for cond = 1:2

    grandmatrix(1,2.*cond-1,loc) = mean(matrix1_corr(:,2.*cond-1,loc));
    grandmatrix(1,2.*cond,loc) = mean(matrix1_corr(:,2.*cond,loc));
    grandmatrix(2,2.*cond-1,loc) = std(matrix1_corr(:,2.*cond-1,loc))/sqrt(15);
    grandmatrix(2,2.*cond,loc) = std(matrix1_corr(:,2.*cond,loc))/sqrt(15);

end
end
%'^-'(blue line, red marker),'^-'(orange line, red marker),'^--'(blue line, red marker),'o--',(orange line, red marker)

subplot(1,2,1),errorbar([0 6.48 9.16], reshape(grandmatrix(1,1,:),[1,3]),reshape(grandmatrix(2,1,:),[1,3]),'^-', 'Color',colours(2,:),...
    'MarkerSize',3,'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',8,'LineWidth',1);
xlim([0 10]), xticks([0 6.48 9.16]),ylim([0 1]);
set(gca,'FontSize',12);
ylabel('Percentage of accurate judgments'),xlabel('Eccentricity (in dva)');
hold on
errorbar([0 6.48 9.16], reshape(grandmatrix(1,2,:),[1,3]),reshape(grandmatrix(2,2,:),[1,3]),'^--','Color',colours(2,:), 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
errorbar([0 6.48 9.16], reshape(grandmatrix(1,3,:),[1,3]),reshape(grandmatrix(2,3,:),[1,3]),'^-','Color',colours(5,:), 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
errorbar([0 6.48 9.16], reshape(grandmatrix(1,4,:),[1,3]),reshape(grandmatrix(2,4,:),[1,3]),'^--','Color',colours(5,:), 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
hold off
title('a)','FontWeight','normal');

%% for hypothesis 4

%congruent
clear indvN
clear indvP
clear indvN_loc
clear indvP_loc

matrix4 = zeros(15, 3);
location = [0 6.48 9.16];
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Congruent_IP,:);
    indvP = Results(Results(:,1)==sub & Find_Congruent_CP,:); % trial classification
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
    matrix4(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

% sort AUC, location and individual ID into a table, then fit lme model
% AUC = reshape(matrix3,[45,1]);
% subjects(1:15,1)= 1:1:15;
% subjects = repmat(subjects,[3,1]);
% location = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];
% 
% 
% data = table(AUC, subjects, location,'VariableName',{'AUC','subject','location'});
% lm1 = fitlme(data, 'AUC ~ location + (location|subject) + (1|subject)');
% lm2 = fitlme(data, 'AUC ~ 1 +(location|subject) + (1|subject) ');
% 
% compare(lm2,lm1);


% visualising the effect by plotting within-subjects errorbar

population_m = mean(mean(mean(matrix4))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix4,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix4_corr = matrix4 - repmat(indv_dff, [1,3]);

%calculate within-subject sem on each location level

for i = 1:3
    se(i)= std(matrix4_corr(:,i))/sqrt(15);
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

% plotting within-subjects errorbar

population_m = mean(mean(mean(matrix5))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix5,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix5_corr = matrix5 - repmat(indv_dff, [1,3]);

%calculate within-subject sem on each location level

for i = 1:3
    se1(i)= std(matrix5_corr(:,i))/sqrt(15);
end

% visualisation for congruent and incongruent condition


subplot(1,2,2),errorbar([0 6.48 9.16], mean(matrix4,1),se,'d-', 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
xlim([0 10]), xticks([0 6.48 9.16]),ylim([0.4 0.8]);
set(gca,'FontSize',12);
ylabel('AUC'),xlabel('Eccentricity (in dva)');
hold on
errorbar([0 6.48 9.16], mean(matrix5,1),se1,'d-', 'MarkerSize',3,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
plot([0 10],[0.5 0.5],'k--');
hold off
title('b)','FontWeight','normal');

%% for hypothesis 3
% percentage distribution
colours = cbrewer('qual', 'Set1', 8); 
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;

Results = [Results eccentricity];
location = [0 6.48 9.16];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N(sub,i) = size(Results(Results(:,end)==location(loc)& Find_N & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            length(Results(Results(:,end)==location(loc)&Find_N & Results(:,1)==sub,:));
        pcg_AP(sub,i) = size(Results(Results(:,end)==location(loc)&(Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            length(Results(Results(:,end)==location(loc)&(Find_CAP|Find_IAP) & Results(:,1)== sub,:));
        end
    end

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = mean(pcg_N,1);
percentage_AP = mean(pcg_AP,1);


    for a = 1:8
        std_hypo3(2.*a -1)= std(pcg_AP(:,a))/sqrt(15);
        std_hypo3(2.*a) = std(pcg_N(:,a))/sqrt(15);
        m_hypo3(2.*a-1) = percentage_AP(:,a);
        m_hypo3(2.*a) = percentage_N(:,a);
        ax(2.*a - 1)= a - 0.14;
        ax(2.*a) = a + 0.14;
    end
    
percentage_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
subplot(2,3,loc),b = bar(percentage_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
b(1).FaceColor = colours(2,:);
b(2).FaceColor = colours(1,:);
xlabel('Decision x confidence'),
ylim([0 0.6]);

hold on
errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',1);
hold off
    if loc ==1
    title('Eccentricity = 0 dva','FontWeight','normal');
    ylabel('Percentage %')
    elseif loc == 2
    title('Eccentricity = 6.48 dva','FontWeight','normal');
    else
    legend({'Patch Present','Patch Absent'},'Location','north');
    title('Eccentricity = 9.16 dva','FontWeight','normal')
    end
    
% clear pcg_AP
% clear pcg_N
% clear percentage_N
% clear percentage_AP
% clear std_hypo3
% clear m_hypo3
% clear percentage_All

end

grandmatrix = zeros(2,2);
grandmean1 = sum(Results(Find_CAP|Find_IAP|Find_N,8)==1)/length(Results(Find_CAP|Find_IAP|Find_N,:));

matrix1 = zeros(90,4);
location = [0 6.48 9.16];

for loc = 1:3
for sub = 1:15
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
    rowindex = size(matrix1,1)-sum(matrix1(:,1)==0)+1;
    
    matrix1(rowindex,1) = sub;
    matrix1(rowindex,2) = condition - 1;
    matrix1(rowindex,3) = location(loc);
    matrix1(rowindex,4) = accuracy;
    
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

end


% calculate standard errors for hit rate
matrix2 = [matrix1(matrix1(:,2)==0 & matrix1(:,3)==0,4) matrix1(matrix1(:,2)==0 & matrix1(:,3)==6.48,4)...
          matrix1(matrix1(:,2)==0 & matrix1(:,3)==9.16,4)];

grandmean1= repmat(mean(mean(matrix2)),15,1);
individualmean1 = mean(matrix2,2);
indi_diff1 = individualmean1 - grandmean1;
indi_diff1 = [indi_diff1 indi_diff1 indi_diff1];
matrix2_corr = matrix2 - indi_diff1;

for i = 1:3
    grandmatrix1(1,i) = mean(matrix2_corr(:,i));
    grandmatrix1(2,i) = std(matrix2_corr(:,i))/sqrt(15);
end

% calculate standard errors for CR
matrix3 = [matrix1(matrix1(:,2)==1 & matrix1(:,3)==0,4) matrix1(matrix1(:,2)==1 & matrix1(:,3)==6.48,4)...
          matrix1(matrix1(:,2)==1 & matrix1(:,3)==9.16,4)];

grandmean2= repmat(mean(mean(matrix3)),15,1);
individualmean2 = mean(matrix3,2);
indi_diff2 = individualmean2 - grandmean2;
indi_diff2 = [indi_diff2 indi_diff2 indi_diff2];
matrix3_corr = matrix3 - indi_diff2;

for i = 1:3
    grandmatrix2(1,i) = mean(matrix3_corr(:,i));
    grandmatrix2(2,i) = std(matrix3_corr(:,i))/sqrt(15);
end

subplot(2,2,3),b = errorbar([0 6.48 9.16],grandmatrix1(1,:),grandmatrix1(2,:),'o-',...
    'MarkerSize',0.8,'LineWidth',1.2);
b.Color = colours(2,:);
hold on
r = errorbar([0 6.48 9.16],grandmatrix2(1,:),grandmatrix2(2,:),'o-',...
    'MarkerSize',0.8,'LineWidth',1.2);
r.Color = colours(1,:);
hold off

set(gca,'FontSize',12);
ylim([0.4 1]), xticks([0 6.48 9.16]),ylabel('Percentage of accurate judgment'),xlabel('Eccentricity (in dva)'),...
title('b)','FontSize',14);
    

% calculate AUC on each location level

location = [0 6.48 9.16];


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

% sort AUC, location and individual ID into a table, then fit lme model
AUC = reshape(matrix3,[45,1]);
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[3,1]);
location = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];


data = table(AUC, subjects, location,'VariableName',{'AUC','subject','location'});
lm1 = fitlme(data, 'AUC ~ location + (location|subject) + (1|subject)');
lm2 = fitlme(data, 'AUC ~ 1 +(location|subject) + (1|subject) ');

compare(lm2,lm1);

% visualising the effect by plotting within-subjects errorbar

population_m = mean(mean(mean(matrix3))); % obtain grand mean across participant & condition
population_m = population_m + zeros(15,1);
indv_mean = mean(matrix3,2); % obtain individual mean
indv_dff = indv_mean - population_m; 
matrix3_corr = matrix3 - repmat(indv_dff, [1,3]);

%calculate within-subject sem on each location level

for i = 1:3
    se(i)= std(matrix3_corr(:,i))/sqrt(15);
end

subplot(2,2,4),errorbar([0 6.48 9.16], mean(matrix3,1),se,'d-', 'MarkerSize',2,'MarkerEdgeColor','red',...
    'MarkerFaceColor','red','LineWidth',1,'CapSize',8);
xlim([0 10]), xticks([0 6.48 9.16]),ylim([0.4 1]);
set(gca,'FontSize',12);
ylabel('AUC'),xlabel('Eccentricity (in dva)');
title('c)','FontSize',14);
hold on
plot([0 10], [0.5 0.5],'k--');
hold off

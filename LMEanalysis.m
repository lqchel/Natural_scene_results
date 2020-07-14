% this code builds LME models for the relationship between RGB differences
% and discrimination performance

clear all
Results = importdata('AUC and RGB matrix.mat');

data = table(categorical(Results(:,1)),Results(:,2),Results(:,3),categorical(Results(:,4)),'VariableNames',{'Image','RGB','AUC','PatchType'});




lm1 = fitlme(data,'AUC~ RGB*PatchType + (1|Image)');
lm2 = fitlme(data,'AUC ~ RGB:PatchType + PatchType + (1|Image)');
compare(lm2,lm1)
lm3 = fitlme(data,'AUC ~ RGB + RGB:PatchType + (1|Image)');
compare(lm3,lm1)

%% difference patch types, analyse the effect of RGB separately

% present v.s. randomly-selected only
data1 = table(categorical(Results(Results(:,4)==0,1)),Results(Results(:,4)==0,2),Results(Results(:,4)==0,3),...
    'VariableNames',{'Image','RGB','AUC'});
lm3 = fitlme(data1,'AUC~ RGB + (1|Image)');
lm4 = fitlme(data1,'AUC ~ (1|Image)');
compare(lm4,lm3)

% present, critical object absent v.s. modified absent 
data2 = table(categorical(Results(Results(:,4)==1,1)),Results(Results(:,4)==1,2),Results(Results(:,4)==1,3),...
    'VariableNames',{'Image','RGB','AUC'});
lm5 = fitlme(data2,'AUC~ RGB + (1|Image)');
lm6 = fitlme(data2,'AUC ~ (1|Image)');
compare(lm6,lm5)

% original vs. modified
data3 = table(categorical(Results(Results(:,4)==2,1)),Results(Results(:,4)==2,2),Results(Results(:,4)==2,3),...
    'VariableNames',{'Image','RGB','AUC'});
lm7 = fitlme(data3,'AUC~ RGB + (1|Image)');
lm8 = fitlme(data3,'AUC ~ (1|Image)');
compare(lm8,lm7)

subplot(2,1,1),gscatter(Results(:,2),Results(:,3),Results(:,4),'grb','.',[7 7]);
legend('off');
xlabel('RGB difference'),ylabel('AUC'),title('Raw patch-pair AUC vs. RGB difference distribution');
ylim([0 1]),yticks([0:0.2:1]);
set(gca,'FontSize',12);
subplot(2,1,2),gscatter(data.RGB,fitted(lm1),data.PatchType,'grb','.',[7 7]);
legend('present vs. randomly-selected absent','present, critical object absent vs. modified absent',...
    'present critical object vs. modified absent');
% legend('off');
xlabel('RGB difference'),ylabel('AUC'),title('LME fitted values of AUC vs. RGB difference');
ylim([0 1]),yticks([0:0.2:1]);
set(gca,'FontSize',12);

%% build model for patch pairs with equivalent RGB differences only

Results_filtered = Results(Results(:,2)<= 9906386,:);
data4 = table(categorical(Results_filtered(:,1)),Results_filtered(:,2),Results_filtered(:,3),categorical(Results_filtered(:,4)),...
    'VariableNames',{'Image','RGB','AUC','PatchType'});

lm9 = fitlme(data4,'AUC~ RGB*PatchType + (1|Image)');
lm10 = fitlme(data4,'AUC ~ PatchType + RGB:PatchType + (1|Image)');
PD = compare(lm10,lm9)
lm11 = fitlme(data4,'AUC ~ RGB + RGB:PatchType + (1|Image)');
PairType = compare(lm11,lm9)
lm12 = fitlme(data4,'AUC ~ RGB + PatchType + (1|Image)');
interaction = compare(lm12,lm9)

%plot the effect
subplot(1,2,1),gscatter(Results_filtered(:,2),Results_filtered(:,3),Results_filtered(:,4),'grb','.',[7 7]);
legend('off');
xlabel('RGB difference'),ylabel('AUC'),title('Raw patch-pair AUC vs. RGB difference distribution');
ylim([0 1]),yticks([0:0.2:1]);
set(gca,'FontSize',12);
subplot(1,2,2),gscatter(data4.RGB,fitted(lm9),data4.PatchType,'grb','.',[7 7]);
legend('present vs. randomly-selected absent','present, critical object absent vs. modified absent',...
    'present critical object vs. modified absent');
% legend('off');
xlabel('RGB difference'),ylabel('AUC'),title('LME fitted values of AUC vs. RGB difference');
ylim([0 1]),yticks([0:0.2:1]);
set(gca,'FontSize',12);

%% build model for patch pairs only
data = table(categorical(Results(:,1)),Results(:,2),Results(:,3),categorical(Results(:,4)),'VariableNames',{'Image','RGB','AUC','PatchType'});

lm4 = fitlme(data,'AUC~ PatchType + (1|Image)');
lm2 = fitlme(data,'AUC ~ RGB:PatchType + PatchType + (1|Image)');
compare(lm2,lm1)
lm3 = fitlme(data,'AUC ~ RGB + RGB:PatchType + (1|Image)');
compare(lm3,lm1)

%% anlysing the effect of RGB difference on each eccentricity level


Results = data1(data1.PresentLoc==data1.AbsentLoc,:);

location1 = (Results.AbsentLoc==2 | Results.AbsentLoc==4| Results.AbsentLoc==6|Results.AbsentLoc== 8) .* 6.50;
location2 = (Results.AbsentLoc==1 | Results.AbsentLoc==3| Results.AbsentLoc==7|Results.AbsentLoc== 9) .* 9.20;
eccentricity = zeros(height(Results),1)+ location1 + location2;
Results = addvars(Results, eccentricity,'NewVariableNames',{'Eccentricity'});
location = [0 6.50 9.20];


level_results1 = Results(Results.Eccentricity == 0,:);
data0 = table(categorical(level_results1.Image),level_results1.RGB,level_results1.AUC,'VariableNames',{'Image','RGB','AUC'});
lm1 = fitlme(data0,'AUC~RGB + (1|Image)');
lm_test1 = fitlme(data0,'AUC~ 1+(1|Image)');
compare(lm_test1,lm1)

level_results2 = Results(Results.Eccentricity == 6.50,:);
data2 = table(categorical(level_results2(:,1)),level_results2(:,2),level_results2(:,3),'VariableNames',{'Image','RGB','AUC'});
lm2 = fitlme(data2,'AUC~RGB + (1|Image)');
lm_test2 = fitlme(data2,'AUC~ 1+(1|Image)');
compare(lm_test2,lm2)

level_results3 = Results(Results.Eccentricity == 9.20,:);
data3 = table(categorical(level_results3(:,1)),level_results3(:,2),level_results3(:,3),'VariableNames',{'Image','RGB','AUC'});
lm3 = fitlme(data3,'AUC~RGB + (1|Image)');
lm_test3 = fitlme(data3,'AUC~1+(1|Image)');
compare(lm_test3,lm3)


% plotting the results
colours = cbrewer('qual','Set2',8);

subplot(2,2,1),scatter(data0.RGB,data0.AUC,7,colours(1,:),'filled');
xlabel('RGB difference'),ylabel('AUC');
set(gca,'FontName','Arial','FontSize',12);
title('Raw data, Eccentricity = 0');

subplot(2,2,2),scatter(data2.RGB,data2.AUC,7,colours(2,:),'filled');
xlabel('RGB difference'),ylabel('AUC');
set(gca,'FontName','Arial','FontSize',12);
title('Raw data, Eccentricity = 6.50');

subplot(2,2,3),scatter(data3.RGB,data3.AUC,7,colours(3,:),'filled');
xlabel('RGB difference'),ylabel('AUC');
set(gca,'FontName','Arial','FontSize',12);
title('Raw data, Eccentricity = 9.20');

subplot(2,2,4),scatter(data1.RGB,fitted(lm1),7,colours(1,:),'filled');
hold on
scatter(data2.RGB, fitted(lm2),7,colours(2,:),'filled');
scatter(data3.RGB,fitted(lm3),7,colours(3,:),'filled');
hold off
xlabel('RGB difference'),ylabel('AUC');
set(gca,'FontName','Arial','FontSize',12);
title('Fitted data');
legend({'0','6.50','9.20'});



%% fitting lme model to confidence x decision to absent patches only, rather than AUC
% Results = data_sub1; %% uncomment for one subject data

Results = data1(data1.PresentLoc == data1.AbsentLoc & data1.PatchType == 2,:);
location1 = (Results.AbsentLoc==2 | Results.AbsentLoc==4| Results.AbsentLoc==6|Results.AbsentLoc== 8) .* 6.50;
location2 = (Results.AbsentLoc==1 | Results.AbsentLoc==3| Results.AbsentLoc==7|Results.AbsentLoc== 9) .* 9.20;
eccentricity = zeros(height(Results),1)+ location1 + location2;
Results = addvars(Results, eccentricity,'NewVariableNames',{'Eccentricity'});
location = [0 6.50 9.20];

level_results1 = Results(Results.Eccentricity == 0,:);
data4 = table(categorical(level_results1.Image),level_results1.RGB,level_results1.Decision,'VariableNames',{'Image','RGB','Decision'});
lm4 = fitlme(data4,'Decision~RGB + (1|Image)');
lm_test4 = fitlme(data4,'Decision~ 1+(1|Image)');
compare(lm_test4,lm4)

level_results2 = Results(Results.Eccentricity == 6.50,:);
data5 = table(categorical(level_results2.Image),level_results2.RGB,level_results2.Decision,'VariableNames',{'Image','RGB','Decision'});
lm5 = fitlme(data5,'Decision~RGB + (1|Image)');
lm_test5 = fitlme(data5,'Decision~ 1+(1|Image)');
compare(lm_test5,lm5)

level_results3 = Results(Results.Eccentricity == 9.20,:);
data6 = table(categorical(level_results3.Image),level_results3.RGB,level_results3.Decision,'VariableNames',{'Image','RGB','Decision'});
lm6 = fitlme(data6,'Decision~RGB + (1|Image)');
lm_test6 = fitlme(data6,'Decision~1+(1|Image)');
compare(lm_test6,lm6)

% plot raw data and fitted values
colours = cbrewer('qual','Set2',8);

subplot(2,2,1),scatter(data4.RGB,data4.Decision,7,colours(1,:),'filled');
xlabel('RGB difference'),ylabel('Decision x Confidence');
set(gca,'FontName','Arial','FontSize',12);
title('Raw data, Eccentricity = 0');

subplot(2,2,2),scatter(data5.RGB,data5.Decision,7,colours(2,:),'filled');
xlabel('RGB difference'),ylabel('Decision x Confidence');
set(gca,'FontName','Arial','FontSize',12);
title('Raw data, Eccentricity = 6.50');

subplot(2,2,3),scatter(data6.RGB,data6.Decision,7,colours(3,:),'filled');
xlabel('RGB difference'),ylabel('Decision x Confidence');
set(gca,'FontName','Arial','FontSize',12);
title('Raw data, Eccentricity = 9.20');

subplot(2,2,4),scatter(data4.RGB,fitted(lm4),7,colours(1,:),'filled');
hold on
scatter(data5.RGB, fitted(lm5),7,colours(2,:),'filled');
scatter(data6.RGB,fitted(lm6),7,colours(3,:),'filled');
hold off
xlabel('RGB difference'),ylabel('Decision x Confidence');
set(gca,'FontName','Arial','FontSize',12);
title('Fitted data');
legend({'0','6.50','9.20'});
ylim([-4 4])

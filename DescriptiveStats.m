clear all
Results = importdata('Pooled Results.mat');
Results(:,9)= Results(:,9)-0.5;
Results(:,9) = Results(:,8).*Results(:,9);

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

%% patch difference
PDN = Results(Find_N,11);
PDC = Results(Find_Incongruent_CP,11);
PD = zeros(16,2);
t = 1000000;
for a = 1:16
    range = [(a-1).*t a.*t];
    if a == 1
    PD(a,1)= size(PDC(PDC<range(2)&PDC>range(1)),1)/size(PDC,1);
    PD(a,2) = size(PDN(PDN<range(2)&PDN>range(1)),1)/size(PDN,1);
    else
    PD(a,1)= PD(a-1,1)+size(PDC(PDC<range(2)&PDC>range(1)),1)/size(PDC,1);
    PD(a,2)= PD(a-1,2)+size(PDN(PDN<range(2)&PDN>range(1)),1)/size(PDN,1);
    end
end

% bar([1:1:16],PD,'stacked'),title('Normalized cumulative histogram for RGB difference'),xlabel('Patch dffierence(x 1000000)'),...
%     ylabel('Cumulative probability');

PD = [0 0;PD];
plot([0:1:16],PD(:,1),'bo-','LineWidth',1.5),title('Normalized cumulative plot for RGB difference'),xlabel('Patch dffierence(x 1000000)'),...
    ylabel('Cumulative probability'),xlim([0 16]),xticks([0:1:16]);
hold on
plot([0:1:16],PD(:,2),'ro-','LineWidth',1.5);
hold off
legend('Critical patches','N patches');

histogram(Results(Find_Incongruent_CP,11)),title('sum RGB difference of critical object patches'),...
    xlabel('RGB difference'), ylabel('Count');
histogram(Results(Find_Incongruent_CP|Find_N,11)),title('sum RGB difference of critical object and N patches'),...
    xlabel('RGB difference'), ylabel('Count');

%% mean and std of confidence ratings
rc = Results(:,9).*(Results(:,8));
rc_m = round(mean(rc),2);
rc_std = round(std(rc),2);
histogram(rc),...
    xlabel('Response x confidence'), ylabel('Count'),title(['Response x confidence distribution, mean = ',...
    num2str(rc_m),' std = ',num2str(rc_std)]);

conf_m = round(mean(Results(:,9)),2);
std_m = round(std(Results(:,9)),2);
histogram(Results(:,9)),xticks([1:1:4]),title(['Confidence distribution, mean = ',num2str(conf_m),' std = ',num2str(std_m)]);

%% percentage correct over question number

PooledR = importdata('Pooled & Cleaned Results.mat');


for i = 1:21
    Results = PooledR(PooledR(:,3)==i,:);
    
    %Find trials that presented N patches -- signal absent for hypo 1
    %Find trials that presented N patches -- signal absent for hypo 1
    Find_N = Results(:,5) ==1;

    % present patch trials -- signal present for hypo 1
    Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
    Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

    %Incongruent trial with congruent object, congruent trial with incongruent
    %object -- signal absent for hypo 2

    Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
    Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
    
    present_correct = (Find_CAP|Find_IAP) & Results(:,8)==1;
    absent_correct = (Find_N|Find_Congruent_IP|Find_Incongruent_CP) & Results(:,8)==-1;
    
    pcg(i) = (length(Results(present_correct|absent_correct,:))/length(Results)).* 100;
    conf(i) = round(mean(Results(:,9)),2);
    se(i) =std(Results(:,9)) / sqrt(length(Results(:,9)));
end

bar([1:1:21],pcg),xticks([1:1:21]),ylim([0 100]),title('Percentage correct for 1st-21st patch'),xlabel('Patch number'),...
    ylabel('% Percentage correct');
bar([1:1:21],conf),xticks([1:1:21]),ylim([0 4]),title('Average confidence for 1st-21st patch'),xlabel('Patch number'),...
    ylabel('Confidence level');

% errorbar([1:1:21],conf,se,'b.-','MarkerSize',5,'LineWidth',1),xlim([0 21]),ylim([0 4]),title('Average confidence for 1st-21st patch'),xlabel('Patch number'),...
%     ylabel('Confidence level');

%% percentage correct over question number, separating conditions

PooledR = importdata('Pooled & Cleaned Results.mat');


for i = 1:21
    Results = PooledR(PooledR(:,3)==i,:);
    
    %Find trials that presented N patches -- signal absent for hypo 1
    %Find trials that presented N patches -- signal absent for hypo 1
    Find_N = Results(:,5) ==1;

    % present patch trials -- signal present for hypo 1
    Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
    Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

    %Incongruent trial with congruent object, congruent trial with incongruent
    %object -- signal absent for hypo 2

    Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
    Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
    
    %percentage correct and confidence level for N patch, present patch,
    %and critical lures
    pcg(i,1) = length(Results(Find_N&Results(:,8)== -1,:))/length(Results);
    pcg(i,2) = length(Results((Find_IAP|Find_CAP)&Results(:,8)==1,:))/length(Results);
    pcg(i,3) = length(Results((Find_Congruent_IP|Find_Incongruent_CP)&Results(:,8)== -1,:))/length(Results);
    conf(i,1) = round(mean(Results(Find_N,9)),2);
    conf(i,2) = round(mean(Results(Find_IAP|Find_CAP,9)),2);
    conf(i,3) = round(mean(Results(Find_Congruent_IP|Find_Incongruent_CP,9)),2);
    se(i) =std(Results(:,9)) / sqrt(length(Results(:,9)));

end
pcg = pcg .*100;
bar([1:1:21],pcg,'stacked'),xticks([1:1:21]),ylim([0 100]),title('Percentage correct for 1st-21st patch'),xlabel('Patch number'),...
    ylabel('% Percentage correct'),legend('N','A and P','Critical object lures');
bar([1:1:21],conf,'stacked'),xticks([1:1:21]),ylim([0 12]),title('Average confidence for 1st-21st patch'),xlabel('Patch number'),...
    ylabel('Confidence level'),legend('N','A and P','Critical object lures');


%% average confidence for levels of periphery
central = Results(:,7)==5;
periphery_1= 2.*(Results(:,7)==2|Results(:,7)==4|Results(:,7)==6|Results(:,7)==8);
periphery_2 = 3.* (Results(:,7)==1|Results(:,7)==3|Results(:,7)==7|Results(:,7)==9);
location = central + periphery_1 + periphery_2;
Results(:,9) = Results(:,8).*Results(:,9);
m1 = mean(Results(location ==1,9));
se1 = std(Results(location ==1,9))/sqrt(length(Results(location ==1,9)));
m2 = mean(Results(location ==2,9));
se2 = std(Results(location ==2,9))/sqrt(length(Results(location ==2,9)));
m3 = mean(Results(location ==3,9));
se3 = std(Results(location ==3,9))/sqrt(length(Results(location ==3,9)));

bar([1 2 3],[m1 m2 m3]);
hold on
errorbar([m1 m2 m3],[se1 se2 se3],'k.');
hold off

%% raw distribution of response values

colours = cbrewer('qual', 'Set1', 8); 
%distribution of all response values
for i = 1:9
    frequency_all(i) = size(Results(Results(:,9)== i-5,:),1); 
end

frequency_all = [frequency_all(1:4) frequency_all(6:9)];
bar([1:8], frequency_all),xticklabels({'-4','-3','-2','-1','1','2','3','4'});

%response value distribution for N patches

for i = 1:9
    frequency_N(i) = size(Results(Find_N & Results(:,9)== i-5,:),1);
end

frequency_N = [frequency_N(1:4) frequency_N(6:9)];
%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});

%response value distribution for AP

for i = 1:9
    frequency_AP(i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5,:),1);
end
%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});

frequency_AP = [frequency_AP(1:4) frequency_AP(6:9)];
frequency_All = [reshape(frequency_AP,[8,1]) reshape(frequency_N,[8,1])];
subplot(1,3,1),b = bar(frequency_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
b(1).FaceColor = colours(2,:);
b(2).FaceColor = colours(1,:);
xlabel('Decision x confidence'),
ylabel('Count');
ylim([0 9000])
title('a)','FontName','Arial','FontWeight','normal');
legend({'Patch present','Patch absent'},'Location','northwest');

%proportion of responses
pcgN_matrix = zeros(15,9);
for sub = 1:15
for i = 1:9
    pcgN_matrix(sub,i) = size(Results(Find_N & Results(:,9)== i-5 & Results(:,1)==sub,:),1)/length(Results(Find_N & Results(:,1)==sub,:));
end
end
percentage_N = mean(pcgN_matrix,1);
percentage_N = [percentage_N(1:4) percentage_N(6:9)];
std_N = std(pcgN_matrix)/sqrt(15);
std_N = [std_N(1:4) std_N(6:9)];
%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});

%accuracy check
for i = 1:9
    percentage_NO(i) = size(Results(Find_N & Results(:,9)== i-5,:),1)/length(Results(Find_N,:));
end

%response value distribution for AP
pcgAP_matrix = zeros(15,9);

for sub = 1:15
for i = 1:9
    pcgAP_matrix(sub,i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
        length(Results((Find_CAP|Find_IAP) & Results(:,1)== sub,:));
end
end

%accuracy check
for i = 1:9
percentage_yes(i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5,:),1)/...
        length(Results((Find_CAP|Find_IAP),:));
end

%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});
percentage_AP = mean(pcgAP_matrix,1);
percentage_AP = [percentage_AP(1:4) percentage_AP(6:9)];
std_AP = std(pcgAP_matrix)/sqrt(15);
std_AP = [std_AP(1:4) std_AP(6:9)];



percentage_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
%manually setting the x, y, and error for error bar
for a = 1:8
    ax(2.*a - 1)= a - 0.14;
    ax(2.*a) = a + 0.14;
    pcg(2.*a-1) = percentage_AP(a);
    pcg(2.*a) = percentage_N(a);
    std_hypo1(2.*a - 1) = std_AP(a);
    std_hypo1(2.*a) = std_N(a);
end

subplot(1,3,2),b = bar(percentage_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
b(1).FaceColor = colours(2,:);
b(2).FaceColor = colours(1,:);
xlabel('Decision x confidence'),
ylabel('Percentage of responses');
title('b)','FontName','Arial','FontWeight','normal');
hold on
errorbar(ax,pcg,std_hypo1,'k.','MarkerSize',1,'LineWidth',1);
hold off
%legend({'Patch present','Patch absent'});


%% critical object distributions, in raw frequency and percentages

for i = 1:9
    frequency_CP(i) = size(Results(Find_Congruent_CP & Results(:,9)== i-5,:),1);
    frequency_CI(i) = size(Results(Find_Congruent_IP & Results(:,9)== i-5,:),1);
    frequency_IP(i) =  size(Results(Find_Incongruent_IP & Results(:,9)== i-5,:),1);
    frequency_II(i) =  size(Results(Find_Incongruent_CP & Results(:,9)== i-5,:),1);
    
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

frequency_CP = [frequency_CP(1:4) frequency_CP(6:9)];
frequency_CI = [frequency_CI(1:4) frequency_CI(6:9)];
frequency_CPI = [reshape(frequency_CP,[8,1]) reshape(frequency_CI,[8,1])];

frequency_IP = [frequency_IP(1:4) frequency_IP(6:9)];
frequency_II = [frequency_II(1:4) frequency_II(6:9)];
frequency_IPI = [reshape(frequency_IP,[8,1]) reshape(frequency_II,[8,1])];

percentage_CP = [percentage_CP(:,1:4) percentage_CP(:,6:9)];
percentage_CI = [percentage_CI(:,1:4) percentage_CI(:,6:9)];
percentage_IP = [percentage_IP(:,1:4) percentage_IP(:,6:9)];
percentage_II = [percentage_II(:,1:4) percentage_II(:,6:9)];

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


subplot(2,2,1), a = bar(frequency_CPI,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
a(1).FaceColor = colours(2,:);
a(2).FaceColor = colours(1,:);
ylabel('Count');
title('Congruent condition');
legend({'Original P1','Different P1'},'Location','northwest');

subplot(2,2,2), c = bar(frequency_IPI,'grouped','BarWidth',1,'LineWidth',1.2);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
c(1).EdgeColor = colours(2,:);
c(2).EdgeColor = colours(1,:);
c(1).FaceColor = 'white';
c(2).FaceColor = 'white';
title('Incongruent condition');
legend({'Original P1','Different P1'},'Location','northwest');

subplot(2,2,3),e = bar(percentage_CPI,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',12);
e(1).FaceColor = colours(2,:);
e(2).FaceColor = colours(1,:);
xlabel('Decision x confidence'),
ylabel('Percentage');
hold on
errorbar(ax,pcg_CPI,std_CPI,'k.','MarkerSize',1,'LineWidth',0.8);
hold off

subplot(2,2,4),f = bar(percentage_IPI,'grouped','BarWidth',1,'LineWidth',1.2);
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


%% error bars for hypothesis 1
colours = cbrewer('qual', 'Pastel1', 8); 
grandmatrix = zeros(2,2);
grandmean1 = mean(Results(Find_IAP|Find_CAP|Find_N,9));
matrix1 = zeros(15,2);

for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
    indv_mean = mean(R_indv(:,9));
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
    else
        Find_patch = R_indv(:,5)==1;
    end
    
    Results_P = R_indv(Find_patch,:);
    condition_mean = mean(Results_P(:,9))-indv_mean + grandmean1;
    matrix1(sub,condition) = condition_mean;
    
    clear Results_P
    clear condition_mean 

end
end

for i = 1:2
    grandmatrix(1,i) = mean(matrix1(:,i));
    grandmatrix(2,i) = std(matrix1(:,i))/sqrt(15);
end



%% connecting individual data points
for sub = 1:15
    plot([1 2],[matrix1(sub,:) matrix2(sub,:)],'bo');
    hold on
    if sub == 15
        hold off
    end
end
xlim([0 3]), xticks([1:2]),ylim([-4 4]);

%% how many sigstar??

Results_h1 = Results(Find_IAP|Find_CAP|Find_N,:);
Response = Results_h1(:,9);
patch = Results_h1(:,5)~= 1;
data = table(Response, patch,Results_h1(:,1),'VariableName',{'Response','patch','subjects'});
lm1 = fitlme(data,'Response ~ patch + (1|subjects)');

%% draw graph

subplot(1,3,3),d = bar(grandmatrix(1,:),'BarWidth',0.7);
ylabel('Decision x confidence'),...
    ylim([-4 4]);
set(gca,'XTickLabel',{'Present','Absent'},'FontSize',12);
%H = sigstar({{'Patch present','Patch absent'}},0.001);
d.FaceColor = colours(5,:);
title('c)');
hold on
errorbar(grandmatrix(1,:),grandmatrix(2,:),'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off



%% error bars for hypo 2


grandmatrix2 = zeros(2,2);
grandmean2 = mean(Results(Results(:,6)==3,9));
matrix2= zeros(15,2);

for sub = 1:15
    R_indv = Results(Results(:,1)==sub & Results(:,6)==3,:);
    indv_mean = mean(R_indv(:,9));
for condition = 1:2
    if condition == 1
        Find_patch = (R_indv(:,4)==0 & R_indv(:,5)==2) | (R_indv(:,4)==1 & R_indv(:,5)==3);
    else
        Find_patch = (R_indv(:,4)==0 & R_indv(:,5)==3) | (R_indv(:,4)==1 & R_indv(:,5)==2);
    end
    
    Results_P = R_indv(Find_patch,:);
    condition_mean = mean(Results_P(:,9))-indv_mean + grandmean2;
    matrix2(sub,condition) = condition_mean;
    
    clear Results_P
    clear condition_mean 

end
end

for i = 1:2
    grandmatrix2(1,i) = mean(matrix2(:,i));
    grandmatrix2(2,i) = std(matrix2(:,i))/sqrt(15);
end

histogram(matrix2(:,1));
hold on
histogram(matrix2(:,2));
hold off
xlim([-4 4]), xticks([-4:1:4])
%% how many sigstar??

Results_h2 = Results(Results(:,6)==3,:);
Response = Results_h2(:,9);
patch = (Results_h2(:,4)==0 & Results_h2(:,5)==2) | (Results_h2(:,4)==1 & Results_h2(:,5)==3);
data2 = table(Response, patch,Results_h2(:,1),'VariableName',{'Response','patch','subjects'});
lm2 = fitlme(data2,'Response ~ patch + (1|subjects)');
lm3 = fitlme(data2, 'Response ~ 1+(1|subjects)')
compare(lm3,lm2)

[h,p] = ttest(Results_h2(patch,9),Results_h2(~patch,9));

sum(patch & Results_h2(:,8)==1)/length(Results_h2(patch,:));
sum(~patch & Results_h2(:,8)== -1)/length(Results_h2(~patch,:))

%% plot the graph
d = bar(grandmatrix2(1,:),'BarWidth',0.7);
ylabel('Decision x confidence'),...
    ylim([-4 4]);
set(gca,'XTickLabel',{'Same','Different'},'FontSize',14);
H = sigstar({{'Same','Different'}},0.001);
d.FaceColor = colours(2,:);

hold on
errorbar(grandmatrix2(1,:),grandmatrix2(2,:),'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off

%AUC distribution

histogram(matrix1,'BinNumber',5)


%% raincloud plot for descriptives

colours = cbrewer('qual', 'Set1', 6); 

 fig_position = [200 200 600 400];
 f7 = figure('Position', fig_position);
h1 = raincloud_plot(Results(Find_N,9), 'box_on', 1, 'color', colours(2,:), 'alpha', 0.6, 'bandwidth', .2,'density_type', 'ks',...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(Results(Find_CAP|Find_IAP,9), 'box_on', 1, 'color', colours(6,:), 'alpha', 0.6,'bandwidth', .2,'density_type', 'ks',...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
legend([h1{1} h2{1}], {'Patch present', 'Patch absent'});
title(['Individual mean responses for present and absent patches']),...
    xlabel('Response value (decision x confidence)'),ylabel('Density'),...
    xticklabels({'-4','-3','-2','-1','1','2','3','4'});
box off

b(1).FaceColor = colours(6,:);
b(2).FaceColor = colours(2,:);

%% hit and false alarm

% hypothesis 1

%sum(Results(Find_CAP|Find_IAP,8)== 1)/length(Results(Find_CAP|Find_IAP,:))

colours = cbrewer('qual', 'Pastel1', 8); 
grandmatrix = zeros(2,2);
grandmean1 = sum(Results(Find_CAP|Find_IAP|Find_N,8)==1)/length(Results(Find_CAP|Find_IAP|Find_N,:));
matrix1 = zeros(15,2);

for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
    indv_mean = sum(R_indv(:,8)==1)/length(R_indv);
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
    else
        Find_patch = R_indv(:,5)==1;
    end
    
    Results_P = R_indv(Find_patch,:);
    condition_mean = (sum(Results_P(:,8)==1)/length(Results_P))-indv_mean + grandmean1;
    matrix1(sub,condition) = condition_mean;
    
    clear Results_P
    clear condition_mean 

end
end

for i = 1:2
    grandmatrix(1,i) = mean(matrix1(:,i));
    grandmatrix(2,i) = std(matrix1(:,i))/sqrt(15);
end

d = bar(grandmatrix(1,:),'BarWidth',0.7);
ylabel('Percentage of present judgment');
ylim([0 1]);
set(gca,'XTickLabel',{'Patch present','Patch absent'},'FontSize',14);
d.FaceColor = colours(2,:);
H = sigstar({{'Patch present','Patch absent'}},0.001);
hold on
errorbar(grandmatrix(1,:),grandmatrix(2,:),'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off

Results_h1 = Results(Find_CAP|Find_IAP|Find_N,:);
data = table(Results_h1(:,8)==1,Results_h1(:,5)==1,Results_h1(:,1),'VariableName',{'Response','Presence','Subjects'});
[h,p,ci] = ttest(matrix1(:,1),matrix1(:,2))
 
% for hypothesis 2
%% error bars for hypo 2


grandmatrix2 = zeros(2,2);
grandmean2 = sum(Results(Results(:,6)==3,8)==1)/length(Results(Results(:,6)==3,:));
matrix2= zeros(15,2);

for sub = 1:15
    R_indv = Results(Results(:,1)==sub & Results(:,6)==3,:);
    indv_mean = sum(R_indv(:,8)==1)/length(R_indv);
for condition = 1:2
    if condition == 1
        Find_patch = (R_indv(:,4)==0 & R_indv(:,5)==2) | (R_indv(:,4)==1 & R_indv(:,5)==3);
    else
        Find_patch = (R_indv(:,4)==0 & R_indv(:,5)==3) | (R_indv(:,4)==1 & R_indv(:,5)==2);
    end
    
    Results_P = R_indv(Find_patch,:);
    condition_mean = (sum(Results_P(:,8)==1)/length(Results_P))-indv_mean + grandmean2;
    matrix2(sub,condition) = condition_mean;
    
    clear Results_P
    clear condition_mean 

end
end

for i = 1:2
    grandmatrix2(1,i) = mean(matrix2(:,i));
    grandmatrix2(2,i) = std(matrix2(:,i))/sqrt(15);
end

d = bar(grandmatrix2(1,:),'BarWidth',0.7);
ylabel('Percentage of present judgment');
ylim([0 1]);
set(gca,'XTickLabel',{'Original object','Replaced object'},'FontSize',14);
d.FaceColor = colours(2,:);
H = sigstar({{'Original object','Replaced object'}},0.001);
hold on
errorbar(grandmatrix2(1,:),grandmatrix2(2,:),'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off
[h,p]= ttest(matrix2(:,1),matrix2(:,2));

%% combining the two hypothesis

d = bar([grandmatrix(1,1) grandmatrix2(1,2) grandmatrix(1,2)],'BarWidth',0.7);
ylabel('Percentage of present judgment');
ylim([0 1]);
set(gca,'XTickLabel',{'Present','Modified','Random'},'FontSize',12);
d.FaceColor = colours(2,:);
hold on
errorbar([grandmatrix(1,1) grandmatrix2(1,2) grandmatrix(1,2)],[grandmatrix(2,1) grandmatrix2(2,2) grandmatrix(2,2)],'k.',...
    'MarkerSize',2,'LineWidth',1.2);
hold off
H = sigstar({{'Present','Modified'}},0.001);
H = sigstar({{'Present','Random'}},0.001);

%% percentage of yes judgement along question number
decay_matrix = zeros(3,21);
for q = 1:21
    decay_matrix(1,q) = size(Results((Find_CAP|Find_IAP) & Results(:,3)== q & Results(:,8)==1,:),1)/size(Results((Find_CAP|Find_IAP)...
        & Results(:,3)== q,:),1);
    decay_matrix(2,q) = size(Results((Find_Congruent_IP|Find_Incongruent_CP) & Results(:,3)== q & Results(:,8)==1,:),1)/size(Results((Find_Congruent_IP|Find_Incongruent_CP)...
        & Results(:,3)== q,:),1);
    decay_matrix(3,q)= size(Results(Find_N & Results(:,3)== q & Results(:,8)==1,:),1)/size(Results(Find_N...
        & Results(:,3)== q,:),1);
    
    
end

for condition = 1:3
    subplot(1,2,1),plot([1:1:21],decay_matrix(condition,:),'.-','LineWidth',1,'MarkerSize',8);
    hold on
    if condition == 3
        hold off
    end
end
legend('off');
ylim([0 1]),xlim([1 21]),xticks([1:4:21]),ylabel('Percentage of present judgment'),xlabel('Question number');
set(gca,'FontSize',12);

%% mean confidence ratings along question number

decay_matrix = zeros(3,21,15);

for q = 1:21
    for sub = 1:15
    decay_matrix(1,q,sub) = nanmean(Results((Find_CAP|Find_IAP) & Results(:,3)== q & Results(:,1)==sub,9));
    decay_matrix(2,q,sub) = nanmean(Results((Find_Congruent_IP|Find_Incongruent_CP) & Results(:,3)== q & Results(:,1)==sub,9));
    decay_matrix(3,q,sub)= nanmean(Results(Find_N & Results(:,3)== q & Results(:,1)==sub,9));
    end
end

condition_mean = nanmean(decay_matrix,3);
se_matrix = zeros(3,21);

for condition = 1:3
    indv_mean = nanmean(decay_matrix(condition,:,:),2);
    grand_mean = nanmean(nanmean(nanmean(decay_matrix(condition,:,:))))+zeros(1,1,15);
    indv_diff = indv_mean - grand_mean;
    decay_matrix(condition,:,:) = decay_matrix(condition,:,:)- repmat(indv_diff, [1,21,1]);
    for num = 1:21
    se_matrix(condition,num) = nanstd(decay_matrix(condition,num,:))/sqrt(15);
    end
    clear indv_mean
    clear grand_mean
    clear indv_diff
end
    
for condition = 1:3
    %plot([1:1:21],decay_matrix(condition,:),'.-','LineWidth',1,'MarkerSize',8);
    errorbar([1:1:21],condition_mean(condition,:),se_matrix(condition,:),'.-','LineWidth',1,'MarkerSize',8);
    hold on
    if condition == 3
        hold off
    end
end
legend({'Present','Similar','Absent'});
ylim([-3.5 3.5]),yticks([-3.5:1:3.5]),xlim([1 21]),xticks([1:4:21]),ylabel('Decision x onfidence'),xlabel('Question number');
set(gca,'FontSize',12);

%% percentage of yes judgement along trial number

decay_matrix = zeros(3,80);
for q = 1:80
    decay_matrix(1,q) = size(Results((Find_CAP|Find_IAP) & Results(:,2)== q & Results(:,8)==1,:),1)/size(Results((Find_CAP|Find_IAP)...
        & Results(:,2)== q,:),1);
    decay_matrix(2,q) = size(Results((Find_Congruent_IP|Find_Incongruent_CP) & Results(:,2)== q & Results(:,8)==1,:),1)/size(Results((Find_Congruent_IP|Find_Incongruent_CP)...
        & Results(:,2)== q,:),1);
    decay_matrix(3,q)= size(Results(Find_N & Results(:,2)== q & Results(:,8)==1,:),1)/size(Results(Find_N...
        & Results(:,2)== q,:),1);
    
    
end

for condition = 1:3
    plot([1:1:80],decay_matrix(condition,:),'d-');
    hold on
    if condition == 3
        hold off
    end
end
legend({'Present','Similar','Absent'});
ylim([0 1]),ylabel('Percentage of present judgment'),xlabel('Trial number');
set(gca,'FontSize',12);


%% hypothesis 3
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
subplot(1,3,loc),b = bar(percentage_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',14);
b(1).FaceColor = colours(2,:);
b(2).FaceColor = colours(1,:);
xlabel('Decision x confidence'),
ylim([0 0.6]);

hold on
errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',1);
hold off
    if loc ==1
    legend({'Patch Present','Patch Absent'},'Location','north');
    title('Eccentricity = 0 dva');
    ylabel('Percentage %')
    elseif loc == 2
    title('Eccentricity = 6.48 dva');
    else
    title('Eccentricity = 9.16 dva')
    end
    
% clear pcg_AP
% clear pcg_N
% clear percentage_N
% clear percentage_AP
% clear std_hypo3
% clear m_hypo3
% clear percentage_All

end

%% hypothesis 4
% congruent
location = [0 6.48 9.16];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N(sub,i) = size(Results(Results(:,end)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)& Find_Congruent_IP & Results(:,1)==sub,:),1);
        pcg_AP(sub,i) = size(Results(Results(:,end)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)&Find_Congruent_CP & Results(:,1)== sub,:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = mean(pcg_N,1);
percentage_AP = mean(pcg_AP,1);

if loc == 1
    disp(pcg_N)
    disp(pcg_AP)
    disp(percentage_N)
    disp(percentage_AP)
end
    for a = 1:8
        std_hypo3(2.*a -1)= std(pcg_AP(:,a))/sqrt(15);
        std_hypo3(2.*a) = std(pcg_N(:,a))/sqrt(15);
        m_hypo3(2.*a-1) = percentage_AP(:,a);
        m_hypo3(2.*a) = percentage_N(:,a);
        ax(2.*a - 1)= a - 0.14;
        ax(2.*a) = a + 0.14;
    end
    
percentage_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
subplot(1,3,loc),b = bar(percentage_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',14);
b(1).FaceColor = colours(2,:);
b(2).FaceColor = colours(1,:);
ylim([0 0.7]);

hold on
errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',1);
hold off
    if loc ==1
    legend({'Original P1','Different P1'},'Location','northwest');
    title('Eccentricity = 0 dva');
    elseif loc == 2
    title('Eccentricity = 6.48 dva');
    else
    title('Eccentricity = 9.16 dva')
    end
    
clear pcg_AP
clear pcg_N
clear percentage_N
clear percentage_AP
clear std_hypo3
clear m_hypo3
clear percentage_All

end


%incongruent

for loc = 1:3
    
    for i = 1:9
        percentage_N(i) = size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5,:),1)/length(Results(Results(:,end)==location(loc)&Find_Incongruent_CP,:));
        percentage_AP(i) = size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5,:),1)/length(Results(Results(:,end)==location(loc)& Find_Incongruent_IP,:));
    end

percentage_N = [percentage_N(1:4) percentage_N(6:9)];

percentage_AP = [percentage_AP(1:4) percentage_AP(6:9)];
percentage_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
subplot(1,3,loc),b = bar(percentage_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',14);
b(1).EdgeColor = colours(2,:);
b(2).EdgeColor = colours(1,:);
xlabel('Decision x confidence'),

ylim([0 0.7]),yticks([0:0.1:0.8]);
   
end

for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N(sub,i) = size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,1)==sub,:),1);
        pcg_AP(sub,i) = size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,1)== sub,:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = mean(pcg_N,1);
percentage_AP = mean(pcg_AP,1);

if loc == 1
    disp(pcg_N)
    disp(pcg_AP)
    disp(percentage_N)
    disp(percentage_AP)
end
    for a = 1:8
        std_hypo3(2.*a -1)= std(pcg_AP(:,a))/sqrt(15);
        std_hypo3(2.*a) = std(pcg_N(:,a))/sqrt(15);
        m_hypo3(2.*a-1) = percentage_AP(:,a);
        m_hypo3(2.*a) = percentage_N(:,a);
        ax(2.*a - 1)= a - 0.14;
        ax(2.*a) = a + 0.14;
    end
    
percentage_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
subplot(2,3,loc),b = bar(percentage_All,'grouped','BarWidth',1,'LineWidth',1.2);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',14);
b(1).EdgeColor = colours(2,:);
b(2).EdgeColor = colours(1,:);
b(1).FaceColor = 'white';
b(2).FaceColor = 'white';
ylim([0 0.7]);
xlabel('Decision x confidence');
hold on
errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',0.8);
hold off
 
% clear pcg_AP
% clear pcg_N
% clear percentage_N
% clear percentage_AP
% clear std_hypo3
% clear m_hypo3
% clear percentage_All

end

%% confidence distribution

Confidence = Results(:,9).*Results(:,8);
confidence_h1P = Confidence(Find_CAP|Find_IAP,:);
confidence_h1A = Confidence(Find_N,:);

for c = 1:4
    frequency_h1P(c,1)=sum(confidence_h1P == c);
    frequency_h1A(c,1) = sum(confidence_h1A == c);
end
frequency_h1 = [frequency_h1P frequency_h1A];
t = bar([1 2 3 4],frequency_h1,'grouped','BarWidth',1);
t(1).FaceColor= colours(:,2);
t(2).FaceColor = colours(:,1);
set(gca,'FontSize',12);
ylabel('Count'),xlabel('Confidence');
legend({'Patch present','Patch absent'},'Location','northwest');

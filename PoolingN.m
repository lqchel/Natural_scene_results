%% check for confidence ratings
clear all;
Results = importdata('Pooled Results.mat');
Results(:,9) = Results(:,8).*Results(:,9);

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

% % hypo 1, not by congruency
% Results_N = Results(Find_N,:);
% Results_AP = Results(Find_IAP|Find_CAP,:);
% 
% for i = -4:4
%     Counts_N(i+5) = 100.*length(Results_N(Results_N(:,9)== i,:))/length(Results_N);
%     Counts_AP(i+5) = 100.*length(Results_AP(Results_AP(:,9) == i, :))/length(Results_AP);
% end
% 
% subplot(1,2,1), N_plot1 = plot([-4 -3 -2 -1], Counts_N(1:4),'LineWidth',2); xlabel('Response x confidence'), ylabel('%Percentage of responses'),...
%     title('Percentages of "No" to N and AP across confidence levels'),...
%     xticks([-4:1:1]),ylim([0 50]);
%     set(gca,'FontSize',16);
%     axis square;
% hold on
% AP_plot1 = plot([-4 -3 -2 -1], Counts_AP(1:4),'LineWidth',1.5);
% hold off
% legend([N_plot1,AP_plot1], 'No to N (CR)','No to AP (Miss)');
% 
% subplot(1,2,2),N_plot2 = plot([1 2 3 4], Counts_N(6:9),'LineWidth',1.5); xlabel('Response x confidence'), ylabel('%Percentage of responses'),...
%     title('Percentages of "Yes" to N and AP across confidence levels'),...
%     xticks([1:1:4]),ylim([0 50]);
%     set(gca,'FontSize',16);
%     axis square;
% hold on
% AP_plot2 = plot([1 2 3 4], Counts_AP(6:9),'LineWidth',2);
% hold off
% legend([N_plot2,AP_plot2], 'Yes to N (FA)','Yes to AP (Hit)');

%% pooling N patches, and get bits per second
Results_pooled = importdata('Pooled Results.mat');
Results_pooled(:,9) = Results_pooled(:,8).*Results_pooled(:,9);
% signal present and absent trial classification
Idv_bits = zeros(1, 80);
for q = 1:80
    qnum = q + 2;
    
Results = Results_pooled(Results_pooled(:,2) == qnum,:);
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

%% hypo 1, not by congruency
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);
disp(length(Results_AP));
% type 1 ROC curve analysis

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
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

Idv_AUC(q) = AUC;

% if q == 12
%     subplot(1,2,1), plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC curve for Image 12 A/P discrimination AUC = ', num2str(AUC)]);
%     axis square;
% elseif q == 24
%     subplot(1,2,2), plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC curve for Image 24 A/P discrimination AUC = ', num2str(AUC)]);
%     axis square;
% end

Critical_patch = Results_AP(Results_AP(:,6)==3,:);
disp(Critical_patch)

if isinteger(Critical_patch(1,7)/2)
    bits = log2(2.*AUC).*(4 + length(Results_N));
    Idv_bits(1,q) = bits/0.133;
else
    bits = log2(2.*AUC).*(5 + length(Results_N));
    Idv_bits(1,q) = bits/0.133;
end



Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];
bits = [];
end

average_bits = mean(Idv_bits);

%histogram(Idv_AUC), title('Image AUC distribution'), xticks([0 0.69 0.83 0.96 1]),xlabel('AUC'), ylabel('Frequency');
%histogram(Idv_AUC), title('Image AUC distribution'), ylim([0 1]), xlabel('AUC'), ylabel('Frequency');
% 
% %histogram(Results_pooled(:,9))
% % plot([1:1:80],Idv_bits,'LineWidth',1.2),ylim([0 1800])
% colours = cbrewer('qual', 'Set1', 8); 
% h = histogram(Idv_bits);
% h.FaceColor = colours(2,:);
% xlim([0 1900]),ylim([0 35]),xlabel('Bits per second'), ylabel('Number of images'),set(gca,'FontSize',12);
% hold on
% plot([40 40],[0 40],'r--','LineWidth',1.8);
% hold off

plot(Idv_bits)

subplot(2,1,1), plot(1:80,Idv_AUC,'bo-');
ylim([0.5 1]), xlabel('Image number'), ylabel('AUC');
title('AUC of image 1-80, old calculation method');


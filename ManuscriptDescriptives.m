clear all
Results = importdata('Pooled Results.mat');
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

%% hypothesis 1
%proportion of responses
pcgN_matrix = zeros(3,9);
for sub = 1:3
for i = 1:9
    pcgN_matrix(sub,i) = size(Results(Find_N & Results(:,9)== i-5 & Results(:,1)==sub,:),1)/length(Results(Find_N & Results(:,1)==sub,:));
end
end

pcgN_matrix = [pcgN_matrix(:,1:4) pcgN_matrix(:,6:9)];
se_N = within_se(pcgN_matrix,size(pcgN_matrix,1),size(pcgN_matrix,2));
m_N = mean(pcgN_matrix);


%response value distribution for Present patches
pcgAP_matrix = zeros(3,9);

for sub = 1:3
for i = 1:9
    pcgAP_matrix(sub,i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
        length(Results((Find_CAP|Find_IAP) & Results(:,1)== sub,:));
end
end
pcgAP_matrix = [pcgAP_matrix(:,1:4) pcgAP_matrix(:,6:9)];

%accuracy check
% % for i = 1:9
% % percentage_yes(i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5,:),1)/...
% %         length(Results((Find_CAP|Find_IAP),:));
% % end

%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});
m_AP = mean(pcgAP_matrix);
se_AP = within_se(pcgAP_matrix,3,8);


m_All = [reshape(m_AP,[8,1]) reshape(m_N,[8,1])];
%manually setting the x, y, and error for error bar
for a = 1:8
    ax(2.*a - 1)= a - 0.14;
    ax(2.*a) = a + 0.14;
    m_hypo1(2.*a-1) = m_AP(:,a);
    m_hypo1(2.*a) = m_N(a);
    se_hypo1(2.*a - 1) = se_AP(a);
    se_hypo1(2.*a) = se_N(a);
end

%% 
colours = cbrewer('qual', 'Set1', 8); 
subplot(3,4,1),b = bar(m_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'[]'},'FontSize', 10);
b(1).FaceColor = colours(6,:);
b(2).FaceColor = colours(5,:);
xlabel('Decision x confidence'),
ylabel('Percentage of responses'),ylim([0 0.6]);
title('Across eccentricities','FontName','Arial');

hold on
errorbar(ax,m_hypo1,se_hypo1,'k.','MarkerSize',1,'LineWidth',1);
hold off
legend({'Present patch','Null patch'});

%% hypothesis 1 on eccentricity 
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;

Results = [Results eccentricity];
location = [0 6.48 9.16];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:3
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

se_AP = within_se(pcg_AP,3,8);
se_N =  within_se(pcg_N,3,8); 

    for a = 1:8
        std_hypo3(2.*a-1) = se_AP(:,a);
        std_hypo3(2.*a) = se_N(:,a);
        m_hypo3(2.*a-1) = percentage_AP(:,a);
        m_hypo3(2.*a) = percentage_N(:,a);
        ax(2.*a - 1)= a - 0.14;
        ax(2.*a) = a + 0.14;
    end

m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];

subplot(3,4,loc+1),b = bar(m_All,'grouped','BarWidth',1);
xticks([]);
b(1).FaceColor = colours(6,:);
b(2).FaceColor = colours(5,:);
xlabel('Decision x confidence'),
ylim([0 0.6]);

if loc ~= 1
yticks([]);
end

hold on
errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',1);
hold off
    if loc ==1
    title('Eccentricity = 0 dva','FontWeight','normal');
    ylabel(' ');
    elseif loc == 2
    title('Eccentricity = 6.50 dva','FontWeight','normal');
    else
    title('Eccentricity = 9.20 dva','FontWeight','normal');
    end
    legend('off');
 
% clear pcg_AP
% clear pcg_N
% clear percentage_N
% clear percentage_AP
% clear std_hypo3
% clear m_hypo3
% clear percentage_All

end


%% hypohtesis 2

for i = 1:9
   
    for sub = 1:3
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

se_CP = within_se(percentage_CP,3,8);
se_CI = within_se(percentage_CI,3,8);
se_IP = within_se(percentage_IP,3,8);
se_II = within_se(percentage_II,3,8);

for a = 1:8
    std_CPI(2.*a-1) = se_CP(:,a);
    std_CPI(2.*a) = se_CI(:,a);
    std_IPI(2.*a-1) = se_IP(:,a);
    std_IPI(2.*a) = se_II(:,a);
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



subplot(3,4,5),e = bar(percentage_CPI,'grouped','BarWidth',1);
e(1).FaceColor = colours(2,:);
e(2).FaceColor = colours(1,:);
ylim([0 0.7]);
hold on
errorbar(ax,pcg_CPI,std_CPI,'k.','MarkerSize',1,'LineWidth',0.8);
hold off
legend({'Congruent original', 'Congruent modified'});

subplot(3,4,9),f = bar(percentage_IPI,'grouped','BarWidth',1,'LineWidth',1.2);
xticks([]);
f(1).FaceColor = 'white';
f(2).FaceColor = 'white';
f(1).EdgeColor = colours(2,:);
f(2).EdgeColor = colours(1,:);
ylim([0 0.7])
hold on
errorbar(ax,pcg_IPI,std_IPI,'k.','MarkerSize',1,'LineWidth',0.8);
hold off
legend({'Incongruent original', 'Incongruent modified'});


%% hypothesis 2 on eccentricity
% congruent
location = [0 6.48 9.16];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:3
        pcg_N(sub,i) = size(Results(Results(:,end)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)& Find_Congruent_IP & Results(:,1)==sub,:),1);
        pcg_AP(sub,i) = size(Results(Results(:,end)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)&Find_Congruent_CP & Results(:,1)== sub,:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = nanmean(pcg_N,1);
percentage_AP = nanmean(pcg_AP,1);
se_AP = within_se(pcg_AP,3,8);
se_N = within_se(pcg_N,3,8);

if loc == 1
    disp(pcg_N)
    disp(pcg_AP)
    disp(percentage_N)
    disp(percentage_AP)
end
    for a = 1:8
        std_hypo3(2.*a -1)= se_AP(:,a);
        std_hypo3(2.*a) = se_N(:,a);
        m_hypo3(2.*a-1) = percentage_AP(:,a);
        m_hypo3(2.*a) = percentage_N(:,a);
        ax(2.*a - 1)= a - 0.14;
        ax(2.*a) = a + 0.14;
    end
    
m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
subplot(3,4,loc+5),b = bar(m_All,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 10);
b(1).FaceColor = colours(2,:);
b(2).FaceColor = colours(1,:);
ylim([0 0.7]);

if loc ~= 1
    yticks([]);
end

hold on
errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',1);
hold off

if loc ==1
     legend({'Original P1','Different P1'},'Location','northwest');
 end
% %     title('Eccentricity = 0 dva');
%     elseif loc == 2
%     title('Eccentricity = 6.48 dva');
%     else
%     title('Eccentricity = 9.16 dva')
%     end
    
clear pcg_AP
clear pcg_N
clear percentage_N
clear percentage_AP
clear std_hypo3
clear m_hypo3
clear percentage_All

end


%incongruent

% for loc = 1:3
%     
%     for i = 1:9
%         percentage_N(i) = size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5,:),1)/length(Results(Results(:,end)==location(loc)&Find_Incongruent_CP,:));
%         percentage_AP(i) = size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5,:),1)/length(Results(Results(:,end)==location(loc)& Find_Incongruent_IP,:));
%     end
% 
% percentage_N = [percentage_N(1:4) percentage_N(6:9)];
% 
% percentage_AP = [percentage_AP(1:4) percentage_AP(6:9)];
% m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
% subplot(3,4,loc+9),b = bar(m_All,'grouped','BarWidth',1);
% set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',14);
% b(1).EdgeColor = colours(2,:);
% b(2).EdgeColor = colours(1,:);
% xlabel('Decision x confidence'),
% 
% ylim([0 0.7]),yticks([0:0.1:0.8]);
%    
% end

for loc = 1:3
    
    for i = 1:9
        for sub = 1:3
        pcg_N(sub,i) = size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,1)==sub,:),1);
        pcg_AP(sub,i) = size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,1)== sub,:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = nanmean(pcg_N,1);
percentage_AP = nanmean(pcg_AP,1);
se_AP = within_se(pcg_AP,3,8);
se_N = within_se(pcg_N,3,8);

% if loc == 1
%     disp(pcg_N)
%     disp(pcg_AP)
%     disp(percentage_N)
%     disp(percentage_AP)
% end
    for a = 1:8
        std_hypo3(2.*a -1)= se_AP(:,a);
        std_hypo3(2.*a) = se_N(:,a);
        m_hypo3(2.*a-1) = percentage_AP(:,a);
        m_hypo3(2.*a) = percentage_N(:,a);
        ax(2.*a - 1)= a - 0.14;
        ax(2.*a) = a + 0.14;
    end
    
m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
subplot(3,4,loc+9),b = bar(m_All,'grouped','BarWidth',1,'LineWidth',1.2);
set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize',10);
b(1).EdgeColor = colours(2,:);
b(2).EdgeColor = colours(1,:);
b(1).FaceColor = 'white';
b(2).FaceColor = 'white';
ylim([0 0.7]);
    
if loc ~= 1
    yticks([]);
end
% xlabel('Decision x confidence');
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


%% trying out new function

data = importdata('OnlinePilotData.mat');
massive_descriptives(data,10);






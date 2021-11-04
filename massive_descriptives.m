%%% this function produces the descriptive statistics of Qianchen's Natural
%%% scene patch experiment. for the figure produced, x axis is Decision x
%%% confidence category (-4 to 4), and y axis is percentage of responses.
%%% data = the data matrix, num_sub = number of participants.


function [out] = massive_descriptives(data)
Results = data(data(:,11)~=0,:);
num_sub = length(unique(Results(:,1)));
Results(:,9) = Results(:,8).*Results(:,9);

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1; % Results(:,6) == 1 for experiment 2
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;

% uncomment for exp 1 results
% ecc_level_1 = Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8;
% ecc_level_2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 2;
% ecc_levels = zeros(length(Results),1)+ ecc_level_1 + ecc_level_2;
% Results = [Results ecc_levels];


%% hypothesis 1
%proportion of responses
pcgN_matrix = zeros(num_sub,9);
subject_code = unique(Results(:,1));
for sub = 1:num_sub
for i = 1:9
    pcgN_matrix(sub,i) = size(Results(Find_N & Results(:,9)== i-5 & Results(:,1)==subject_code(sub),:),1)/length(Results(Find_N & Results(:,1)==subject_code(sub),:));
end
end

pcgN_matrix = [pcgN_matrix(:,1:4) pcgN_matrix(:,6:9)];

if num_sub>=3
se_N = std(pcgN_matrix)/sqrt(num_sub);
end

m_N = mean(pcgN_matrix,1);


%response value distribution for AP
pcgAP_matrix = zeros(num_sub,9);

for sub = 1:num_sub
for i = 1:9
    pcgAP_matrix(sub,i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== subject_code(sub),:),1)/...
        length(Results((Find_CAP|Find_IAP) & Results(:,1)== subject_code(sub),:));
end
end
pcgAP_matrix = [pcgAP_matrix(:,1:4) pcgAP_matrix(:,6:9)];


%bar([1:8], frequency_N),xticklabels({'-4','-3','-2','-1','1','2','3','4'});
m_AP = mean(pcgAP_matrix,1);

if num_sub>=3
se_AP = std(pcgAP_matrix)/sqrt(num_sub);
end

% m_All = [reshape(m_AP,[8,1]) reshape(m_N,[8,1])];
% %manually setting the x, y, and error for error bar
% for a = 1:8
%     ax(2.*a - 1)= a - 0.16;
%     ax(2.*a) = a + 0.16;
%     m_hypo1(2.*a-1) = m_AP(:,a);
%     m_hypo1(2.*a) = m_N(a);
%     se_hypo1(2.*a - 1) = se_AP(a);
%     se_hypo1(2.*a) = se_N(a);
% end


colours = cbrewer('qual', 'Set1', 8); 
out = figure;

if num_sub<3
    
subplot(3,4,1),plot([1:8],m_AP,'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),'MarkerFaceColor',colours(5,:),'LineWidth',1);

hold on
plot(1:8,m_N,'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta','LineWidth',1);
hold off

else

subplot(3,4,1),errorbar([1:8],m_AP,se_AP,'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),'MarkerFaceColor',colours(5,:),'LineWidth',1);

hold on
errorbar(1:8,m_N,se_N,'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta','LineWidth',1);
hold off

end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]);
xlim([0.5 8.5]);
title('Across eccentricities','FontName','Arial','FontSize',12);
legend({'Present + original','Null'},'Location','northwest','FontSize',12);
legend('boxoff');

% colours = cbrewer('qual', 'Set1', 8); 
% out = figure;
% subplot(3,4,1),b = bar(m_All,'grouped','BarWidth',1);
% set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
% b(1).FaceColor = colours(6,:);
% b(2).FaceColor = colours(5,:);
% ylabel('Percentage of responses'),ylim([0 0.6]);
% title('Across eccentricities','FontName','Arial');
% 
% hold on
% errorbar(ax,m_hypo1,se_hypo1,'k.','MarkerSize',1,'LineWidth',1);
% hold off
% legend({'Present patch','Null patch'});

%% hypothesis 1 on eccentricity 

location = [0 1 2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
        pcg_N(sub,i) = size(Results(Results(:,13)==location(loc)& Find_N & Results(:,9)== i-5 & Results(:,1)== subject_code(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_N & Results(:,1)==subject_code(sub),:),1);
        pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)&(Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== subject_code(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&(Find_CAP|Find_IAP) & Results(:,1)== subject_code(sub),:),1);
        end
    end

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];
 
percentage_N = mean(pcg_N,1);
percentage_AP = mean(pcg_AP,1);

if loc == 1
    sum(percentage_AP)
    sum(pcg_AP(1,:))
end

if num_sub>=3
se_AP = std(pcg_AP)/sqrt(num_sub);
se_N =  std(pcg_N)/sqrt(num_sub); 
end

%     for a = 1:8
%         std_hypo3(2.*a-1) = se_AP(:,a);
%         std_hypo3(2.*a) = se_N(:,a);
%         m_hypo3(2.*a-1) = percentage_AP(:,a);
%         m_hypo3(2.*a) = percentage_N(:,a);
%         ax(2.*a - 1)= a - 0.16;
%         ax(2.*a) = a + 0.16;
%     end

m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];

figure(out);

if num_sub<3
    subplot(3,4,loc+1), plot([1:8],percentage_N,'.-', 'MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta',...
    'LineWidth',1);
    % xticks([]);
    hold on
    plot(1:8,percentage_AP, '.-','MarkerSize',12,'MarkerEdgeColor', colours(5,:),'MarkerFaceColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    hold off
else
    subplot(3,4,loc+1), errorbar([1:8],percentage_N,se_N,'.-', 'MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta',...
        'LineWidth',1);

    % xticks([]);
    hold on
    errorbar(1:8,percentage_AP,se_AP, '.-','MarkerSize',12,'MarkerEdgeColor', colours(5,:),'MarkerFaceColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    hold off
end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]);
xlim([0.5 8.5]);

if loc ~= 1
yticks([]);
end

% hold on
% errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',1);
% hold off
    if loc ==1
    title('Foveal','FontWeight','normal','FontSize',12);
    ylabel(' ');
    elseif loc == 2
    title('Peripheral','FontWeight','normal','FontSize',12);
    else
    title('Para-peripheral','FontWeight','normal','FontSize',12);
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
   
    for sub = 1:num_sub
        percentage_CP(sub,i) = size(Results(Results(:,1)== subject_code(sub) & Find_Congruent_CP & Results(:,9)== i-5,:),1)/...
            size(Results(Results(:,1)== subject_code(sub) & Find_Congruent_CP,:),1);
        percentage_CI(sub,i) = size(Results(Results(:,1)== subject_code(sub) &Find_Congruent_IP & Results(:,9)== i-5,:),1)/...
            size(Results(Results(:,1)== subject_code(sub) & Find_Congruent_IP,:),1);
        percentage_IP(sub,i) =  size(Results(Results(:,1)== subject_code(sub) &Find_Incongruent_IP & Results(:,9)== i-5,:),1)/...
            size(Results(Results(:,1)== subject_code(sub) & Find_Incongruent_IP,:),1);
        percentage_II(sub,i) =  size(Results(Results(:,1)== subject_code(sub) & Find_Incongruent_CP & Results(:,9)== i-5,:),1)/...
            size(Results(Results(:,1)== subject_code(sub) & Find_Incongruent_CP,:),1);
    end
end


percentage_CP = [percentage_CP(:,1:4) percentage_CP(:,6:9)];
percentage_CI = [percentage_CI(:,1:4) percentage_CI(:,6:9)];
percentage_IP = [percentage_IP(:,1:4) percentage_IP(:,6:9)];
percentage_II = [percentage_II(:,1:4) percentage_II(:,6:9)];

if num_sub >=3
se_CP = std(percentage_CP)/sqrt(num_sub);
se_CI = std(percentage_CI)/sqrt(num_sub);
se_IP = std(percentage_IP)/sqrt(num_sub);
se_II = std(percentage_II)/sqrt(num_sub);
end

% for a = 1:8
%     std_CPI(2.*a-1) = se_CP(:,a);
%     std_CPI(2.*a) = se_CI(:,a);
%     std_IPI(2.*a-1) = se_IP(:,a);
%     std_IPI(2.*a) = se_II(:,a);
%     pcg_CPI(2.*a-1) = mean(percentage_CP(:,a));
%     pcg_CPI(2.*a) = mean(percentage_CI(:,a));
%     pcg_IPI(2.*a -1) = mean(percentage_IP(:,a));
%     pcg_IPI(2.*a) = mean(percentage_II(:,a));
% end

percentage_CP = mean(percentage_CP,1);
percentage_CI = mean(percentage_CI,1);
percentage_IP = mean(percentage_IP,1);
percentage_II = mean(percentage_II,1);

percentage_CPI = [reshape(percentage_CP,[8,1]) reshape(percentage_CI,[8,1])];
percentage_IPI = [reshape(percentage_IP,[8,1]) reshape(percentage_II,[8,1])];


figure(out);

if num_sub <3
    
    subplot(3,4,5),plot(1:8,percentage_CP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    hold on
    plot(1:8,percentage_CI,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest');
    legend('boxoff')

    subplot(3,4,9),plot(1:8,percentage_IP,'.--','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),'LineWidth',1);
    hold on
    plot(1:8,percentage_II,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerEdgeColor',colours(1,:),'LineWidth',1.2);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest','FontSize',12);
    legend('boxoff')
else
    
    subplot(3,4,5),errorbar(1:8,percentage_CP,se_CP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
        colours(2,:),'LineWidth',1);
    hold on
    errorbar(1:8,percentage_CI,se_CI,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest','FontSize',12);
    legend('boxoff')
    
    subplot(3,4,9),errorbar(1:8,percentage_IP,se_IP,'.--','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),'LineWidth',1);
    hold on
    errorbar(1:8,percentage_II,se_II,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerEdgeColor',colours(1,:),'LineWidth',1.2);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest','FontSize',12);
    legend('boxoff')

end





% subplot(3,4,5),e = bar(percentage_CPI,'grouped','BarWidth',1);
% set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
% e(1).FaceColor = colours(2,:);
% e(2).FaceColor = colours(1,:);
% ylim([0 0.6]);
% hold on
% errorbar(ax,pcg_CPI,std_CPI,'k.','MarkerSize',1,'LineWidth',0.8);
% hold off
% legend({'Congruent original', 'Congruent modified'});
% 
% subplot(3,4,9),f = bar(percentage_IPI,'grouped','BarWidth',1,'LineWidth',1.2);
% set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
% f(1).FaceColor = 'white';
% f(2).FaceColor = 'white';
% f(1).EdgeColor = colours(2,:);
% f(2).EdgeColor = colours(1,:);
% ylim([0 0.6])
% hold on
% errorbar(ax,pcg_IPI,std_IPI,'k.','MarkerSize',1,'LineWidth',0.8);
% hold off
% legend({'Incongruent original', 'Incongruent modified'});


%% hypothesis 2 on eccentricity
% congruent
location = [0 1 2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
        pcg_N(sub,i) = size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== subject_code(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,1)==subject_code(sub),:),1);
        pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== subject_code(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,1)== subject_code(sub),:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = nanmean(pcg_N,1);
percentage_AP = nanmean(pcg_AP,1);

if num_sub >=3
se_AP = nanstd(pcg_AP)/sqrt(num_sub);
se_N = nanstd(pcg_N)/sqrt(num_sub);
end

% 
%     for a = 1:8
%         std_hypo3(2.*a -1)= se_AP(:,a);
%         std_hypo3(2.*a) = se_N(:,a);
%         m_hypo3(2.*a-1) = percentage_AP(:,a);
%         m_hypo3(2.*a) = percentage_N(:,a);
%         ax(2.*a - 1)= a - 0.16;
%         ax(2.*a) = a + 0.16;
%     end
%     
% m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
figure(out)

if num_sub<3
    
    subplot(3,4,loc+5),plot(1:8,percentage_AP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    hold on
    plot(1:8,percentage_N,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
else
    
    subplot(3,4,loc+5),errorbar(1:8,percentage_AP,se_AP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
        colours(2,:),'LineWidth',1);
    hold on
    errorbar(1:8,percentage_N,se_N,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);

end


set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]),xlim([0.5 8.5]);

if loc > 1
    yticks([]);
end
hold off
% 
% subplot(3,4,loc+5),b = bar(m_All,'grouped','BarWidth',1);
% set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
% b(1).FaceColor = colours(2,:);
% b(2).FaceColor = colours(1,:);
% ylim([0 1]);
% 
% if loc ~= 1
%     yticks([]);
% end


    
clear pcg_AP
clear pcg_N
clear percentage_N
clear percentage_AP
clear std_hypo3
clear m_hypo3
clear percentage_All

end



for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
        pcg_N(sub,i) = size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== subject_code(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,1)==subject_code(sub),:),1);
        pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== subject_code(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,1)== subject_code(sub),:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = nanmean(pcg_N,1);
percentage_AP = nanmean(pcg_AP,1);

if num_sub >=3
se_AP = nanstd(pcg_AP)/sqrt(num_sub);
se_N = nanstd(pcg_N)/sqrt(num_sub);
end

if num_sub <3
    subplot(3,4,loc+9),plot(1:8,percentage_AP,'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    hold on
    plot(1:8,percentage_N,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    hold off

else
    subplot(3,4,loc+9),errorbar(1:8,percentage_AP,se_AP,'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
        colours(2,:),'LineWidth',1);
    hold on
    errorbar(1:8,percentage_N,se_N,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    hold off

end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]),xlim([0.5 8.5]);

if loc ~= 1
    yticks([]);
end


% 
%     for a = 1:8
%         std_hypo3(2.*a -1)= se_AP(:,a);
%         std_hypo3(2.*a) = se_N(:,a);
%         m_hypo3(2.*a-1) = percentage_AP(:,a);
%         m_hypo3(2.*a) = percentage_N(:,a);
%         ax(2.*a - 1)= a - 0.16;
%         ax(2.*a) = a + 0.16;
%     end
%     
% m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];
% figure(out);
% subplot(3,4,loc+9),b = bar(m_All,'grouped','BarWidth',1,'LineWidth',1.2);
% set(gca,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial');
% b(1).EdgeColor = colours(2,:);
% b(2).EdgeColor = colours(1,:);
% b(1).FaceColor = 'white';
% b(2).FaceColor = 'white';
% ylim([0 1]);
%     
% if loc ~= 1
%     yticks([]);
% end
% 
% if loc == 2
%     xlabel('Decision x Confidencce');
% end

% hold on
% errorbar(ax,m_hypo3,std_hypo3,'k.','LineWidth',0.8);
% hold off
%  

% clear pcg_AP
% clear pcg_N
% clear percentage_N
% clear percentage_AP
% clear std_hypo3
% clear m_hypo3
% clear percentage_All

end
end

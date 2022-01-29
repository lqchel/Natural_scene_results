%%% this function produces the descriptive statistics of Qianchen's Natural
%%% scene patch experiment. for the figure produced, x axis is Decision x
%%% confidence category (-4 to 4), and y axis is percentage of responses.
%%% data = the data matrix, num_sub = number of participants.


function [out] = massive_descriptives(data)
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results'));
Results = data;
num_sub = length(unique(Results(:,1)));

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



%%% set fixed axes for each small figure
for i = 1:3
ax_1st_col{i} = [0.05 0.1+(i-1).*0.27 0.2 0.21];
end
ax_1st_col = flip(ax_1st_col);

for i = 1:3
    for s = 1:3
        ax_rest_cols{i,s}= [0.32+(s-1).*0.23 0.1+(i-1).*0.27 0.2 0.21];
    end
end
ax_rest_cols = flip(ax_rest_cols,1);

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

colours = cbrewer('qual', 'Set1', 8); 
out = figure;
ax = axes('Position',ax_1st_col{1}); % choose axes position

if num_sub<3
    
plot(ax,[1:8],m_AP,'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),'MarkerFaceColor',colours(5,:),'LineWidth',1);

hold on
plot(ax,1:8,m_N,'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta','LineWidth',1);
hold off

else

errorbar(ax,[1:8],m_AP,se_AP,'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),'MarkerFaceColor',colours(5,:),'LineWidth',1);

hold on
errorbar(ax,1:8,m_N,se_N,'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta','LineWidth',1);
hold off

end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]);
xlim([0.5 8.5]);
title('Across eccentricities','FontName','Arial','FontSize',12);
legend({'Present + original','Null'},'Location','northwest','FontSize',12);
legend('boxoff');
ax = [];
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
 
percentage_N = mean(pcg_N,1)
percentage_AP = mean(pcg_AP,1);

if loc == 1
    sum(percentage_AP)
    sum(pcg_AP(1,:))
end

if num_sub>=3
se_AP = std(pcg_AP)/sqrt(num_sub);
se_N =  std(pcg_N)/sqrt(num_sub); 
end

m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];

figure(out);
ax = axes('Position',ax_rest_cols{1,loc});

if num_sub<3
   plot(ax,[1:8],percentage_N,'.-', 'MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta',...
    'LineWidth',1);
    % xticks([]);
    hold on
    plot(ax,1:8,percentage_AP, '.-','MarkerSize',12,'MarkerEdgeColor', colours(5,:),'MarkerFaceColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
    hold off
else
    errorbar(ax,[1:8],percentage_N,se_N,'.-', 'MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta',...
        'LineWidth',1);

    % xticks([]);
    hold on
    errorbar(ax,1:8,percentage_AP,se_AP, '.-','MarkerSize',12,'MarkerEdgeColor', colours(5,:),'MarkerFaceColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
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

ax = [];
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

percentage_CP = mean(percentage_CP,1);
percentage_CI = mean(percentage_CI,1);
percentage_IP = mean(percentage_IP,1);
percentage_II = mean(percentage_II,1);

percentage_CPI = [reshape(percentage_CP,[8,1]) reshape(percentage_CI,[8,1])];
percentage_IPI = [reshape(percentage_IP,[8,1]) reshape(percentage_II,[8,1])];


figure(out);

if num_sub <3
    
    ax = axes('Position',ax_1st_col{2});
    plot(ax,1:8,percentage_CP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    hold on
    plot(ax,1:8,percentage_CI,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest');
    legend('boxoff')

    ax = axes('Position',ax_1st_col{3});
    plot(ax,1:8,percentage_IP,'.--','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),'LineWidth',1);
    hold on
    plot(ax,1:8,percentage_II,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerEdgeColor',colours(1,:),'LineWidth',1.2);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest','FontSize',12);
    legend('boxoff')
else
    ax = axes('Position',ax_1st_col{2});
    errorbar(ax,1:8,percentage_CP,se_CP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
        colours(2,:),'LineWidth',1);
    hold on
    errorbar(ax,1:8,percentage_CI,se_CI,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest','FontSize',12);
    legend('boxoff')
    
    ax = axes('Position',ax_1st_col{3});
    errorbar(ax,1:8,percentage_IP,se_IP,'.--','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),'LineWidth',1);
    hold on
    errorbar(ax,1:8,percentage_II,se_II,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerEdgeColor',colours(1,:),'LineWidth',1.2);
    set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
    xlim([0.5 8.5]),ylim([0 1]);
    hold off
    legend({'Original','Modified'},'Location','northwest','FontSize',12);
    legend('boxoff')

end
ax = [];

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

figure(out)

if num_sub<3
    ax = axes('Position',ax_rest_cols{2,loc});
    plot(ax,1:8,percentage_AP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    hold on
    plot(ax,1:8,percentage_N,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
else
    ax = axes('Position',ax_rest_cols{2,loc});
   errorbar(ax,1:8,percentage_AP,se_AP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
        colours(2,:),'LineWidth',1);
    hold on
    errorbar(ax,1:8,percentage_N,se_N,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);

end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]),xlim([0.5 8.5]);

if loc > 1
    yticks([]);
end
hold off

    
clear pcg_AP
clear pcg_N
clear percentage_N
clear percentage_AP
clear std_hypo3
clear m_hypo3
clear percentage_All

end

ax = [];

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
    ax = axes('Position',ax_rest_cols{3,loc});
    plot(ax,1:8,percentage_AP,'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    hold on
    plot(ax,1:8,percentage_N,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    hold off

else
    ax = axes('Position',ax_rest_cols{3,loc});
    errorbar(ax,1:8,percentage_AP,se_AP,'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
        colours(2,:),'LineWidth',1);
    hold on
    errorbar(ax,1:8,percentage_N,se_N,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
        'LineWidth',1);
    hold off

end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]),xlim([0.5 8.5]);

if loc ~= 1
    yticks([]);
end


end
end

%%% this function produces the descriptive statistics of Qianchen's Natural
%%% scene patch experiment. for the figure produced, x axis is Decision x
%%% confidence category (-4 to 4), and y axis is percentage of responses.
%%% data = the data matrix, num_sub = number of participants.


function [out] = massive_descriptives(data,select_trials)

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

% Conditional selection defined by select_trials
%%%% select_trials = 0, plot all trials only; = 1, plot first 6 patches of
%%%% in each trial; = 2, plot first 24 trials
if select_trials == 1
    Con_sel = Results(:,3) <= 6;
elseif select_trials == 2
    Con_sel = Results(:,2) <= 27;    
end

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
subject_code = unique(Results(:,1));

% response selection indices: present, null, present patches under selected
% conditions, null under selected conditions

index{1} = Find_IAP | Find_CAP;
index{2} = Find_N;
if select_trials > 0
    index{3} = index{1} & Con_sel;
    index{4} = index{2} & Con_sel;
end

for sub = 1:num_sub
    for i = 1:9
        for n = 1:size(index,2) 
            percentage_matrix{n}(sub,i) = size(Results(index{n} & Results(:,9)== i-5 & Results(:,1)==subject_code(sub),:),1)/...
                length(Results(index{n} & Results(:,1)==subject_code(sub),:));
        end
    end
end

for n = 1:size(index,2)
    percentage_matrix{n} = [percentage_matrix{n}(:,1:4) percentage_matrix{n}(:,6:9)]; % select DxC = -4 to -1 and 1 to 4, omit 0
    se{n} = std(percentage_matrix{n})/sqrt(num_sub); % caculate mean and error bars for plotting the figures
    m{n} = mean(percentage_matrix{n},1);
end

%%%%%%%% plot figure
colours = cbrewer('qual', 'Set1', 8); 
c_a = cbrewer('qual','Set3',12); % customize a more transparent colour pallate for all trial results
c_b = cbrewer('qual','Set2',8);
colours_2 = [c_a(6,:); c_b(4,:); c_a(5,:); c_a(4,:)]; % colour palatte for all trial results

out = figure;
ax = axes('Position',ax_1st_col{1}); % choose axes position

h1 = errorbar(ax,1:8,m{1},se{1},'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
hold on
h2 = errorbar(ax,1:8,m{2},se{2},'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'Color','magenta','LineWidth',1);

if select_trials > 0 % plot additional lines under selected conditions
    
h1.Color = colours_2(1,:); % reduce transparency of all trial results
h1.MarkerEdgeColor = colours_2(1,:);
h1.MarkerFaceColor = colours_2(1,:);

h2.Color = colours_2(2,:);
h2.MarkerEdgeColor = colours_2(2,:);
h2.MarkerFaceColor = colours_2(2,:);

% plot lines of selected trials
errorbar(ax,1:8,m{3},se{3},'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
errorbar(ax,1:8,m{4},se{4},'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'Color','magenta','LineWidth',1);
end

hold off
%%%%%%%%%%%%%%%%%%

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]);
xlim([0.5 8.5]);
title('Across eccentricities','FontName','Arial','FontSize',12);

if select_trials == 0
legend({'Present + original','Null'},'Location','northwest','FontSize',12);
legend('boxoff');
end

ax = [];

clear percentage_matrix
clear m
clear se

%% hypothesis 1 on eccentricity 

location = [0 1 2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
            for n = 1:size(index,2)
                percentage_matrix{n}(sub,i) = size(Results(Results(:,13)==location(loc)& index{n} & Results(:,9)== i-5 & Results(:,1)== subject_code(sub),:),1)/...
                    size(Results(Results(:,13)==location(loc)& index{n}& Results(:,1)==subject_code(sub),:),1);
            end
        end
    end

    for n = 1:size(index,2)
        percentage_matrix{n} = [percentage_matrix{n}(:,1:4) percentage_matrix{n}(:,6:9)];
        m{n} = mean(percentage_matrix{n},1);
        se{n} =  std(percentage_matrix{n})/sqrt(num_sub); 
    end


figure(out);
ax = axes('Position',ax_rest_cols{1,loc});
h3{loc} = errorbar(ax,1:8,m{1},se{1}, '.-','MarkerSize',12,'MarkerEdgeColor', colours(5,:),'MarkerFaceColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
hold on
h4{loc} = errorbar(ax,1:8,m{2},se{2},'.-', 'MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta','LineWidth',1);

if select_trials > 0
h3{loc}.Color = colours_2(1,:);
h3{loc}.MarkerEdgeColor = colours_2(1,:);
h3{loc}.MarkerFaceColor = colours_2(1,:);

h4{loc}.Color = colours_2(2,:);
h4{loc}.MarkerEdgeColor = colours_2(2,:);
h4{loc}.MarkerFaceColor = colours_2(2,:);

errorbar(ax,1:8,m{3},se{3},'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),...
    'MarkerFaceColor',colours(5,:),'LineWidth',1);
hold on
errorbar(ax,1:8,m{4},se{4},'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'Color','magenta','LineWidth',1);
end
hold off

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]);
xlim([0.5 8.5]);

if loc ~= 1
yticks([]);
end

if loc ==1
title('Foveal','FontWeight','normal','FontSize',12);
ylabel(' ');
elseif loc == 2
title('Para-fovea','FontWeight','normal','FontSize',12);
else
title('Peripheral','FontWeight','normal','FontSize',12);
end
legend('off');
 
clear percentage_matrix
clear m
clear se

end

ax = [];
%% hypohtesis 2
% response selection indices: congruent original, congruent modified, incongruent original and incongruent modified. + each category under selected
% conditions
index_2 = {Find_Congruent_CP; Find_Congruent_IP; Find_Incongruent_IP; Find_Incongruent_CP};
if select_trials > 0 
    for n = 1:4
        index_2{n+4} = index_2{n} & Con_sel;
    end
end

for i = 1:9   
    for sub = 1:num_sub        
        for n = 1:size(index_2,1)
            percentage_matrix{n}(sub,i) = size(Results(Results(:,1)== subject_code(sub) & index_2{n} & Results(:,9)== i-5,:),1)/...
                size(Results(Results(:,1)== subject_code(sub) & index_2{n},:),1);
        end
    end
end

for n = 1:size(index_2,1)
    percentage_matrix{n} = [percentage_matrix{n}(:,1:4) percentage_matrix{n}(:,6:9)];
    m{n} = mean(percentage_matrix{n},1);
    se{n} = std(percentage_matrix{n})/sqrt(num_sub);
end


figure(out);

% plot congruent results
ax = axes('Position',ax_1st_col{2});
h5 = errorbar(ax,1:8,m{1},se{1},'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
hold on
h6 = errorbar(ax,1:8,m{2},se{2},'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
    'LineWidth',1);

if select_trials > 0
    h5.Color = colours_2(3,:);
    h5.MarkerEdgeColor = colours_2(3,:);
    h5.MarkerFaceColor = colours_2(3,:);

    h6.Color = colours_2(4,:);
    h6.MarkerEdgeColor = colours_2(4,:);
    h6.MarkerFaceColor = colours_2(4,:);
    
    errorbar(ax,1:8,m{5},se{5},'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    errorbar(ax,1:8,m{6},se{6},'.-','MarkerSize',12,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',...
    colours(1,:),'LineWidth',1);
end
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
xlim([0.5 8.5]),ylim([0 1]);
hold off

if select_trials == 0
legend({'Original','Modified'},'Location','northwest','FontSize',12);
legend('boxoff')
end

% plot incongruent results
ax = axes('Position',ax_1st_col{3});
h7 = errorbar(ax,1:8,m{3},se{3},'.--','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),'LineWidth',1);
hold on
h8 = errorbar(ax,1:8,m{4},se{4},'.--','MarkerSize',12,'Color',colours(1,:),'MarkerEdgeColor',colours(1,:),'LineWidth',1.2);

if select_trials > 0
    h7.Color = colours_2(3,:);
    h7.MarkerEdgeColor = colours_2(3,:);
    h7.MarkerFaceColor = colours_2(3,:);

    h8.Color = colours_2(3,:);
    h8.MarkerEdgeColor = colours_2(4,:);
    h8.MarkerFaceColor = colours_2(4,:);
    
    errorbar(ax,1:8,m{7},se{7},'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    errorbar(ax,1:8,m{8},se{8},'.--','MarkerSize',12,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',...
    colours(1,:),'LineWidth',1);
end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
xlim([0.5 8.5]),ylim([0 1]);
hold off

if select_trials == 0
legend({'Original','Modified'},'Location','northwest','FontSize',12);
legend('boxoff')
end

ax = [];

clear percentage_matrix
clear m
clear se
%% hypothesis 2 on eccentricity
location = [0 1 2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
            for n = 1:size(index_2,1)
                percentage_matrix{n}(sub,i) = size(Results(Results(:,13)==location(loc)& index_2{n} & Results(:,9)== i-5 &...
                    Results(:,1)== subject_code(sub),:),1)/size(Results(Results(:,13)==location(loc)& index_2{n} & Results(:,1)==subject_code(sub),:),1);
            end
        end
    end
    
    for n = 1:size(index_2,1)
        percentage_matrix{n} = [percentage_matrix{n}(:,1:4) percentage_matrix{n}(:,6:9)];
        m{n} = nanmean(percentage_matrix{n},1);
        se{n} = nanstd(percentage_matrix{n})/sqrt(num_sub);
    end



figure(out)

ax = axes('Position',ax_rest_cols{2,loc});
h9{loc} = errorbar(ax,1:8,m{1},se{1},'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
hold on
h10{loc} = errorbar(ax,1:8,m{2},se{2},'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
    'LineWidth',1);

if select_trials > 0
    h9{loc}.Color = colours_2(3,:);
    h9{loc}.MarkerEdgeColor = colours_2(3,:);
    h9{loc}.MarkerFaceColor = colours_2(3,:);

    h10{loc}.Color = colours_2(4,:);
    h10{loc}.MarkerEdgeColor = colours_2(4,:);
    h10{loc}.MarkerFaceColor = colours_2(4,:);
    
    errorbar(ax,1:8,m{5},se{5},'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    errorbar(ax,1:8,m{6},se{6},'.-','MarkerSize',12,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',...
    colours(1,:),'LineWidth',1);
end

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]),xlim([0.5 8.5]);

if loc > 1
    yticks([]);
end
hold off

ax = axes('Position',ax_rest_cols{3,loc});
h11{loc} = errorbar(ax,1:8,m{3},se{3},'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
hold on
h12{loc} = errorbar(ax,1:8,m{4},se{4},'.--','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
    'LineWidth',1);

if select_trials > 0
    h11{loc}.Color = colours_2(3,:);
    h11{loc}.MarkerEdgeColor = colours_2(3,:);
    h11{loc}.MarkerFaceColor = colours_2(3,:);

    h12{loc}.Color = colours_2(4,:);
    h12{loc}.MarkerEdgeColor = colours_2(4,:);
    h12{loc}.MarkerFaceColor = colours_2(4,:);
    
    errorbar(ax,1:8,m{7},se{7},'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
    errorbar(ax,1:8,m{8},se{8},'.--','MarkerSize',12,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',...
    colours(1,:),'LineWidth',1);
end

hold off
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 1]),xlim([0.5 8.5]);

if loc ~= 1
    yticks([]);
end
clear m
clear se
clear percentage_matrix
ax = [];

end



end

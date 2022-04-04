%%%%%%%%%%% this is the master script for all figures and analysis in stage
%%%%%%%%%%% 2 report 


%% load experiment 1 and 2 data
data{1} = get_data2(1);
data{2} = get_data2(2);

%% DxC distribution figure (figure 2)

massive_descriptives(data{1});
massive_descriptives(data{2});

%% Type 1 and Type 2 AUC analysis

% produce figure (Figure 3)
AUC_figures(data);

% AUC analysis 
for exp = 1:2
    sub_AUC = AUC_matrix(data{exp});
    AUC_lme(sub_AUC.hypo1_Type1,1,exp) % present + original vs. null, Type 1 AUC - effect size of eccentricity
    AUC_lme(sub_AUC.hypo1_Type2,1,exp) % present + original vs. null, Type 2 AUC - effect size of eccentricity
    AUC_lme(sub_AUC.hypo2_Type1,2, exp) % original vs. modified, Type 1 - effect sizes of ecc, congruence & interaction
    AUC_lme(sub_AUC.hypo2_Type2,2, exp) % original vs. modified, Type 2 - effect sizes of ecc, congruence & interaction
end

%% follow-up analyses after AUC-LME analysis 

% hit and CR along eccentricity for present + original vs. null

subject_id = unique(Results(:,1));
% null patches
Find_N = Results(:,5) ==1;
% present patch trials
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 
% collate into two columns
index{1} = Find_CAP | Find_IAP;  % present
index{2} = Find_N; % null
location = [0 1 2];
accurate_judgement = [1 -1]; % yes to present, no to null
for in = 1:2
    condition_data = Results(index{in},:);
    b = 1;
    for loc = 1:3
        for sub = 1:length(subject_id)
                current_sub_data = condition_data(condition_data(:,1) == subject_id(sub) & condition_data(:,13) == location(loc),:);
                mt(b,:)= [sub location(loc)  sum(current_sub_data(:,8) == accurate_judgement(in))/length(current_sub_data(:,8))]; 
                b = b+1;
        end
    end
    AUC_lme(mt,1)
    clear mt
end

% comparing between patch response on each DxC category (* on Figure 2) -
% run 'AdditionalTtests.m'

%% Image-based analysis

% tDxC figures
run('Plot_tDxC_distributions.m');

% plot tDxC histogram
colours = cbrewer('qual', 'Set2', 8); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg
num_bins = 20;
h2 = histogram(incong_delta.delta,num_bins,'EdgeAlpha',0,'FaceAlpha',0.8,'FaceColor',colours(2,:));
h2.NumBins = num_bins;
hold on
h1 = histogram(cong_delta.delta,num_bins,'EdgeAlpha',0,'FaceAlpha',0.7,'FaceColor',colours(3,:));
h1.NumBins = num_bins;
hold off
legend('off');
xlabel(['\Delta' 'tDxC']),ylabel('Count');
set(gca,'FontName','Arial','FontSize',14,'Box','off');

% delta tDxC against eccentricity
addpath(genpath('C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\RainCloudPlots-master'));
location = categorical([0 1 2]);
cong = categorical([0 1]);
%[0.1300 0.1100 0.7750 0.8150]
%set(gca,'Color','none','XColor','none','YColor','none');
figure;
set(gcf, 'Color', 'white');
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
% set(gcf, 'RendererMode', 'manual');
% set(gcf, 'Renderer', 'painters');
labels = {'F','P-F','P'};
colours = cbrewer('qual', 'Set2', 8); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg
 rain_spread = 0.1; % jitter amount of scatter plot

for ecc_level = 1 : size(location,2)
   
    subplot(1,3,ecc_level);
    set(gca, 'XAxisLocation', 'top');

 
    % Create raincloud
      rain_scatter = (rand(sum(data.eccentricity == location(ecc_level) & data.congruence == cong(1)), 1) - 0.5) * rain_spread; % jitter for the scatter plot
      h1 = raincloud_plot(data.delta(data.eccentricity == location(ecc_level) & data.congruence == cong(1),:), 'box_on', 1, 'color', colours(3, :), 'alpha', 0.8,...
        'box_dodge', 1, 'box_dodge_amount', 0.7, 'dot_dodge_amount', 0.7,...
        'box_col_match', 0,'density_type','rash');
      h2 = raincloud_plot(data.delta(data.eccentricity == location(ecc_level) & data.congruence == cong(2),:), 'box_on', 1, 'color', colours(2, :), 'alpha', 0.8,...
        'box_dodge', 1, 'box_dodge_amount', 1.5, 'dot_dodge_amount',1.5,...
        'box_col_match', 0,'density_type','rash');

     % adjust raincloud
     h1{1}.ShowBaseLine = 'off'; % hide baseline
     h1{1}.EdgeAlpha = 0; % hide cloud outline
     h2{1}.ShowBaseLine = 'off'; % hide baseline
     h2{1}.EdgeAlpha = 0; % hide cloud outline
     
     % adjust rain
     rain_scatter = reshape(rain_scatter, [1 length(rain_scatter)]); % reshape jitter array dimensions to match with YData of the scatter
     h1{2}.YData = h1{2}.YData + rain_scatter; % set jitter to the data scatter
     h2{2}.YData = h2{2}.YData + rain_scatter;
     
     % adjust boxplots
      h1{3}.Position([2 4]) = [h1{3}.Position(2)-rain_spread/2 h1{3}.Position(4)+ rain_spread]; % set box width to match rain spread
      h1{4}.YData = [h1{4}.YData(1) - rain_spread/2 h1{4}.YData(2) + rain_spread/2]; % set median line to match rain spread
      h2{3}.Position([2 4]) = [h2{3}.Position(2)-rain_spread/2 h2{3}.Position(4)+ rain_spread]; 
      h2{4}.YData = [h2{4}.YData(1) - rain_spread/2 h2{4}.YData(2) + rain_spread/2];

    xlim([-4 6.2]);
    ylim([-1.2 1.2]);
    yticks(0);
    yticklabels(labels{ecc_level})
    view([90 -90]);
    set(gca,'box','off','FontName','Arial','FontSize',14);
end

% set another empty axes for collating the rainclouds
figure('Color','white');
axes;
xlim([-4 6.2]);
ylim([-1.2 1.2]);
yticks([-0.8 0 0.8]);
yticklabels(labels);
xlabel('\DeltatDxC');
set(gca,'box','off','FontName','Arial','FontSize',14);
view([90 -90]);
[h,p] = corrcoef(cong_delta.eccentricity, object_size(:,2));


% image stat estimation and data analysis - refer to 'Exploratory_analysis.m'





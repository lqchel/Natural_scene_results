% raincloud plots for the data

% clear all;
% Results = importdata('Pooled Results.mat');
% 
% %% patch index
% Results(:,9) = Results(:,8).* Results(:,9);
% % signal present and absent trial classification
% 
% %Find trials that presented N patches -- signal absent for hypo 1
% %Find trials that presented N patches -- signal absent for hypo 1
% Find_N = Results(:,5) ==1;
% 
% % present patch trials -- signal present for hypo 1
% Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
% Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 
% 
% %Congruent trial with Congruent object, and incongruent trial with
% %incongruent object -- signal present for hypo 2
% Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
% Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;
% 
% %Incongruent trial with congruent object, congruent trial with incongruent
% %object -- signal absent for hypo 2
% Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
% Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
% 
% %% location type
% central = Results(:,7)==5;
% periphery_1= 2.*(Results(:,7)==2|Results(:,7)==4|Results(:,7)==6|Results(:,7)==8);
% periphery_2 = 3.* (Results(:,7)==1|Results(:,7)==3|Results(:,7)==7|Results(:,7)==9);
% location = central + periphery_1 + periphery_2;
% 
% Results = [Results location];
% 
% %% AP v.s. N
% 
% matrix1 = zeros(15,3,2);
% matrix2 = zeros(15,3,2);
% for cong = 1:2
%    if cong == 1
%        Find_P = Find_CAP;
%        Find_CN = Find_N & Results(:,4) ==0;
%    else
%        Find_P = Find_IAP;
%        Find_CN = Find_N & Results(:,4) ==1;
%    end
%    
%     Results_P = Results(Find_P,:);
%     for sub = 1:15
%         R_indv = Results_P(Results_P(:,1)==sub,:);
%         for condition = 1:3
%             matrix1(sub,condition,cong)= mean(R_indv(R_indv(:,end)== condition,9));
%         end
%     end
%     
%     clear R_indv
%   
% 
%     Results_N = Results(Find_CN, :);
%     for sub = 1:15
%         R_indv = Results_N(Results_N(:,1)==sub,:);
%         for condition = 1:3
%             matrix2(sub,condition,cong)= mean(R_indv(R_indv(:,end)== condition,9));
%         end
%     end
%     
%     clear Results_P
%     clear Results_N
% end
% 
% 
% for i = 1:3
%     for j = 1:2
%         if j == 1
%         dataC{i,j}= matrix1(:,i,1);
%         dataI{i,j}= matrix1(:,i,2);
%         else
%         dataC{i,j}=matrix2(:,i,1);
%         dataI{i,j}= matrix2(:,i,2);
%         end
%     end
% end
% 
% 
% 
% 
% for cong = 1:2
%     if cong == 1
%         congruency = 'congruent';
%     else 
%         congruency = 'incongruent';
%     end
% subplot(1,2,cong),errorbar([0 7.74 10.95],grandmatrix(1,:,cong),grandmatrix(2,:,cong),'bo-',...
%     'MarkerSize',2, 'LineWidth',1), ylim([-4 4]), xticks([0 7.74 10.95]),xlabel('Eccentricity by degree of visual angle'),...
%     ylabel('Confidence x response');
% hold on
% errorbar([0 7.74 10.95],grandmatrix(3,:,cong),grandmatrix(4,:,cong),'ro-',...
%     'MarkerSize',2, 'LineWidth',1);
% hold off
% title(['AP v.s N for ',congruency,' images']);
% legend('AP','N');
% end
% 
% 
% 
% %% congruent v.s. incongruent
% 
% matrix3 = zeros(15,3,2);
% matrix4= zeros(15,3,2);
% 
% for cong = 1:2
%    if cong == 1
%        Find_P = Find_Congruent_CP;
%        Find_CN = Find_Congruent_IP;
%    else
%        Find_P = Find_Incongruent_IP;
%        Find_CN = Find_Incongruent_CP;
%    end
%    
%     Results_P = Results(Find_P,:);
%     
%     for sub = 1:15
%         R_indv = Results_P(Results_P(:,1)==sub,:);
%         for condition = 1:3
%             condition_m = mean(R_indv(R_indv(:,end)== condition,9));
%             matrix1(sub,condition)= condition_m- indv_m + grandmean1;
%         end
%     end
% 
%     for condition = 1:3
%         conditionmean1(condition) = mean(matrix1(:,condition));
%         se1(condition) = std(matrix1(:,condition))/sqrt(15);
%     end
%     
%     clear R_indv
%     clear condition_m
%     clear indv_m
% 
%     Results_N = Results(Find_CN, :);
%     matrix2 = zeros(15,3);
%     grandmean2 = mean(Results_N(:,9));
% 
%     for sub = 1:15
%         R_indv = Results_N(Results_N(:,1)==sub,:);
%         indv_m = mean(R_indv(:,9));
%         for condition = 1:3
%             condition_m = mean(R_indv(R_indv(:,end)== condition,9));
%             matrix2(sub,condition)= condition_m- indv_m + grandmean2;
%         end
%     end
% 
%     for condition = 1:3
%         conditionmean2(condition) = mean(matrix2(:,condition));
%         se2(condition) = std(matrix2(:,condition))/sqrt(15);
%     end
%     grandmatrix(1:4,:,cong)= [conditionmean1; se1; conditionmean2; se2];
%     
%     clear Results_P
%     clear Results_N
%     clear matrix1
%     clear matrix2
% end
% 
% for cong = 1:2
%     if cong == 1
%         congruency = 'congruent';
%     else 
%         congruency = 'incongruent';
%     end
% subplot(1,2,cong),errorbar([0 7.74 10.95],grandmatrix(1,:,cong),grandmatrix(2,:,cong),'bo-',...
%     'MarkerSize',2, 'LineWidth',1), ylim([-4 4]),xticks([0 7.74 10.95]),xlabel('Eccentricity by degree of visual angle'),...
%     ylabel('Confidence x response');
% hold on
% errorbar([0 7.74 10.95],grandmatrix(3,:,cong),grandmatrix(4,:,cong),'ro-',...
%     'MarkerSize',2, 'LineWidth',1);
% hold off
% title(['Critical objects for ',congruency,' images']);
% legend('present','absent');
% end
% 
% %% not by congruency
% %hypothesis 1
% matrix1 = mean(matrix1,3);
% matrix1 = mean(matrix1,2);
% matrix2 = mean(matrix2,3);
% matrix2 = mean(matrix2,2);
% %hypothesis 2
% 
% matrix3 = zeros(15,1);
% matrix4= zeros(15,1);
% 
%    Find_P = Find_Congruent_CP|Find_Incongruent_IP;
%    Find_CN = Find_Congruent_IP|Find_Incongruent_CP;
% 
%     Results_P = Results(Find_P,:);
%     
%     for sub = 1:15
%         R_indv = Results_P(Results_P(:,1)==sub,:);
%         matrix3(sub,1)= mean(R_indv(:,9));
%         clear R_indv
%     end
% 
%     
%     Results_N = Results(Find_CN, :);
%     for sub = 1:15
%         R_indv = Results_N(Results_N(:,1)==sub,:);
%         matrix4(sub,1)= mean(R_indv(:,9));
%         clear R_indv
%     end
% 
%  %% raincloud plots
%  fig_position = [200 200 600 400];
%  f7 = figure('Position', fig_position);
%  cb = cbrewer('qual', 'Pastel1', 8); 
%  ksdensity(Results(Find_CAP|Find_IAP,9));
%  hold on
%  ksdensity(Results(Find_N,9),'Support',[-4 4]);
%  hold off
%  xlabel('Response value (decision x confidence)'),ylabel('Density');
% set(gca,'XLim', [-4 4],'XTick',[-4:1:4],'YLim',[0 0.7]);
%  
% h1 = raincloud_plot(matrix1, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.6,'support',[-4 4],'density_type', 'rash',...
%      'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%      'box_col_match', 0);
% h2 = raincloud_plot(matrix2, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.6,'bandwidth',8,'support',[-4 4],'density_type', 'rash',...
%      'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
% legend([h1{1} h2{1}], {'Patch present', 'Patch absent'});
%     xlabel('Response value (decision x confidence)'),ylabel('Density');
% set(gca,'XLim', [-5 4],'XTick',[-4:1:4],'YLim',[-0.2 0.7]);
% box off
% 
%  f8 = figure('Position', fig_position);
% h1 = raincloud_plot(matrix3, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'bandwidth', .2,'density_type', 'ks',...
%      'box_dodge', 1, 'box_dodge_amount', .20, 'dot_dodge_amount', .20,...
%      'box_col_match', 0);
% h2 = raincloud_plot(matrix4, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,'bandwidth', .2,'density_type', 'ks',...
%      'box_dodge', 1, 'box_dodge_amount', .45, 'dot_dodge_amount', .45, 'box_col_match', 0);
% legend([h1{1} h2{1}], {'Patch present', 'Patch absent'});
% title(['Individual mean responses for critical object patches']),...
%     xlabel('Response value (decision x confidence)'),ylabel('Density');
% set(gca,'XLim', [-4 4.5],'XTick',[-4:1:4],'YLim',[-0.3 1]);
% box off

%% new plots, plotting RGB on each DxC level 
% 
% location = [0 6.50 9.20];
% 
% % How to set multiple baselines - https://stackoverflow.com/questions/44195924/mutiple-area-with-different-baseline-in-matlab
% figure;
% 
% set(gcf, 'Color', [1 1 1]);
% set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
% set(gcf, 'RendererMode', 'manual');
% set(gcf, 'Renderer', 'painters');
% 
% for e = 1:3
% 
% Results = data_sub1(data_sub1.Eccentricity == location(e),:);    
% decision_option = unique(Results.Decision);
% 
% colours = cbrewer('qual', 'Dark2', 8); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg
% 
% base_offset = 37; % 17;
% baselines = (0:-base_offset:-(size(decision_option,1)-1)*base_offset); % baselines of clouds
% rain_spread = 10;
% cloud_rain_dist = 1;
% rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain
% rain_scatter = (rand(size(Results.Decision, 1), 1) - 0.5) * rain_spread;
% 
% handles = cell(size(unique(Results.Decision), 1), 1);
% 
% for decision_cat = 1 : size(decision_option, 1)
%   
%       % this locks the axes limits
%        subplot(1,3,e),axis([0 4 min(data_sub1.RGB)-10^6 max(data_sub1.RGB)+ 10^6]);
%         set(gca, 'XAxisLocation', 'top');
%        hold on;  
%      
%     % Create raincloud
%     handles{decision_cat} = raincloud_plot(Results.RGB(Results.Decision == decision_option(decision_cat)), 'box_on', 1, 'color', colours(decision_cat, :), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
%         'box_col_match', 0);
%     
%     % Shift baseline (cloud)
%     handles{decision_cat}{1}.BaseValue = baselines(decision_cat); % move baseline
%     handles{decision_cat}{1}.YData = handles{decision_cat}{1}.YData + baselines(decision_cat); % move cloud
%     handles{decision_cat}{1}.ShowBaseLine = 'off'; % hide baseline
%     handles{decision_cat}{1}.EdgeAlpha = 0; % hide cloud outline
%     
%     % Shift baseline (rain)
%     handles{decision_cat}{2}.YData = rain_scatter;
%     handles{decision_cat}{2}.YData = handles{decision_cat}{2}.YData + baselines(decision_cat) - rain_offset; % move rain
%     handles{decision_cat}{2}.SizeData = 2; % raindrop size
%     
%     % Shift boxplots
%     handles{decision_cat}{3}.Position([2 4]) = [-rain_spread/2 rain_spread]; % set box width to match rain spread
%     handles{decision_cat}{4}.YData = [-rain_spread/2 rain_spread/2]; % set median line to match rain spread
%     handles{decision_cat}{3}.Position(2) = handles{decision_cat}{3}.Position(2) + baselines(decision_cat) - rain_offset; % move box
%     handles{decision_cat}{4}.YData = handles{decision_cat}{4}.YData + baselines(decision_cat) - rain_offset; % move median line
%     handles{decision_cat}{5}.YData = handles{decision_cat}{5}.YData + baselines(decision_cat) - rain_offset; % move top whisker
%     handles{decision_cat}{6}.YData = handles{decision_cat}{6}.YData + baselines(decision_cat) - rain_offset; % move bot whisker
%     
%     % Hide all axes except the first
%     if decision_cat > 1
%         axis off; % hide all axes except the first
%     end
%     
%     view([-90 90]);
%     
% end
% 
% % Set same axis limits for all axes
% % limits = cell2mat(get(ax,'YLim')); % get both axes limits
% % set(ax,'YLim',[min(limits(:)) max(limits(:))]); % set the same values for both axes
% % linkaxes(ax);
% 
% % Link view angle, axes limits
% %linkprop(ax, {'CameraPosition','CameraUpVector','YLim','XLim'});
% 
% %view([-90 90]); % view([90 -90]); % Swap x and y axes
% 
% 
% % Formatting
% % ylim(ax(1), [min(handles{end}{2}.YData)-2 max(handles{1}{1}.YData)+1]);
% % xlim(ax(1), [min(Results.RGB(:))-0.01 max(Results.RGB(:))+0.01]);
% 
% title(['Eccentricity = ' num2str(location(e))]);
% 
% xlabel('Decision x confidence');
% ylabel('RGB difference');
% set(gca,'YTick', fliplr(baselines), 'YTickLabel', fliplr({'-3.5', '-2.5', '-1.5', '-0.5', '0.5', '1.5', '2.5'}));
% 
% % Draw separating lines
% % for base = max_order+1 : length(baselines)
% %     ypos = baselines(base)+(base_offset/2);
% %     plot([0 1], [ypos ypos], 'k');
% % end
% 
% 
% end

figure;

Results = data_sub1(data_sub1.Eccentricity == 6.50,:);    
decision_option = unique(Results.Decision);
all_pos_decision = [-3.5:1:3.5];

colours = cbrewer('div', 'Spectral', 9); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg

base_offset = 6.2 .* 10^-7; % 17;
baselines = 0:-base_offset:-7*base_offset; % baselines of clouds
rain_spread = 0.5 .* 10^-7;
cloud_rain_dist = 0.2.* 10^-7;
rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain

for decision_cat = 1 : 2 %size(decision_option, 1)
    
    ax(decision_cat) = axes; % Not possible to set multiple baselines within same axis;
    set(gca, 'XAxisLocation', 'top');
    
    if decision_cat > 1
        ax(decision_cat).YLim = ax(1).YLim;
        ax(decision_cat).XLim = ax(1).XLim;
        hold on; % this locks the axes limits
    end
    
    % this locks the axes limits
    
    
    % Create raincloud
    handles{decision_cat} = raincloud_plot(Results.RGB(Results.Decision == decision_option(decision_cat)), 'box_on', 1, 'color', colours(decision_cat+1, :), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    
    % check if need to skip a baseline
    if decision_option(decision_cat) == all_pos_decision(decision_cat)
        index = decision_cat;
    else
        index = decision_cat + 1;
    end
    
    % Shift baseline (cloud)
    handles{decision_cat}{1}.BaseValue = baselines(index); % move baseline
    handles{decision_cat}{1}.YData = handles{decision_cat}{1}.YData + baselines(index); % move cloud
%     handles{decision_cat}{1}.ShowBaseLine = 'off'; % hide baseline
%     handles{decision_cat}{1}.EdgeAlpha = 0; % hide cloud outline
%     
%     % Shift baseline (rain)
%     rain_scatter = (rand(sum(Results.Decision == decision_option(decision_cat)), 1) - 0.5) * rain_spread;
%     handles{decision_cat}{2}.YData = rain_scatter;
%     handles{decision_cat}{2}.YData = reshape(rain_scatter,[1,size(rain_scatter,1)]);
%     handles{decision_cat}{2}.YData = handles{decision_cat}{2}.YData + baselines(index) - rain_offset; % move rain
%     handles{decision_cat}{2}.SizeData = 2; % raindrop size
%     
%     % Shift boxplots
%     handles{decision_cat}{3}.Position([2 4]) = [-rain_spread/2 rain_spread]; % set box width to match rain spread
%     handles{decision_cat}{4}.YData = [-rain_spread/2 rain_spread/2]; % set median line to match rain spread
%     handles{decision_cat}{3}.Position(2) = handles{decision_cat}{3}.Position(2) + baselines(index) - rain_offset; % move box
%     handles{decision_cat}{4}.YData = handles{decision_cat}{4}.YData + baselines(index) - rain_offset; % move median line
%     handles{decision_cat}{5}.YData = handles{decision_cat}{5}.YData + baselines(index) - rain_offset; % move top whisker
%     handles{decision_cat}{6}.YData = handles{decision_cat}{6}.YData + baselines(index) - rain_offset; % move bot whisker
%     
%       axis([min(data_sub1.RGB)-8.*10^6 max(data_sub1.RGB)+ 8.*10^6 min(baselines)-(4.*10^-7) max(baselines)+4.*10^-7]);
%         set(gca, 'XAxisLocation', 'top');
%         
%         if decision_cat < size(decision_option, 1)
%         hold on; 
%         else 
%         hold off
%         end
    
    % Hide all axes except the first
%     if decision_cat > 1
%         axis off; % hide all axes except the first
%     end
    
%     view([-90 90]);
    
end

%% Plot raincloud (for all; single concept, full constellation, big phi)

%uncomment for all-subject results
% load('Second order data matrix.mat');
% data9 = data1(data1.PresentLoc == data1.AbsentLoc & data1.PatchType == 0,:);
% location1 = (data9.AbsentLoc==2 | data9.AbsentLoc==4| data9.AbsentLoc==6|data9.AbsentLoc== 8) .* 6.50;
% location2 = (data9.AbsentLoc==1 | data9.AbsentLoc==3| data9.AbsentLoc==7|data9.AbsentLoc== 9) .* 9.20;
% eccentricity = zeros(height(data9),1)+ location1 + location2;
% data9 = addvars(data9, eccentricity,'NewVariableNames',{'Eccentricity'});

Results = data_allab(data_allab.PatchType == 0& data_allab.Eccentricity == 9.20,:);  

%uncomment for subject 1 results
% load('null_data_sub1.mat');
% Results = data_sub1(data_sub1.Eccentricity == 9.20,:);


decision_option = unique(Results.Decision);
acc_all = cell(1, length(decision_option));

for decision_cat = 1 : length(decision_option)

acc_all{decision_cat} = Results.Contrast(Results.Decision == decision_option(decision_cat));

end

% uncomment when there're less than 3 data points in a category
% exception = acc_all{4};
% acc_all = {acc_all{1:3} acc_all{5:8}};

% How to set multiple baselines - https://stackoverflow.com/questions/44195924/mutiple-area-with-different-baseline-in-matlab

figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
set(gcf, 'RendererMode', 'manual');
set(gcf, 'Renderer', 'painters');

colours = cbrewer('qual', 'Dark2', 8); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg
% % colours when there're less than 3 data points in a category
% 
% colours_1 = cbrewer('qual', 'Dark2', 8);
% colours(1:6,:) = colours_1(1:6,:);
% colours(7,:) = colours_1(8,:);

% uncomment for RGB
% base_offset = 6.2.*10^-7; % 17
% baselines = [0:-base_offset:-1.*(size(acc_all,2)-1).*base_offset]; %uncomment when all response category has >2 data points are there
% % baselines = [0:-base_offset:-2*base_offset -4.*base_offset:-base_offset:-7.*base_offset]; % baselines of clouds
% % baselines = [0:-base_offset:-5*base_offset -7.*base_offset];
% rain_spread = 1.5.*10^-7;
% cloud_rain_dist = 0.5.*10^-7;
% rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain

% % for luminance
% base_offset =  0.025; % 17;
% baselines = 0:-base_offset:-7*base_offset; % baselines of clouds
% rain_spread = 0.005;
% cloud_rain_dist = 0.0015;
% rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain
% location = [0 6.50 9.20];
% all_pos_decision = -3.5:1:3.5;

% for contrast
base_offset =  2.* 10^-4; % 17;
baselines = 0:-base_offset:-7*base_offset; % baselines of clouds
rain_spread = base_offset/6;
cloud_rain_dist = rain_spread/3;
rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain
location = [0 6.50 9.20];
all_pos_decision = -3.5:1:3.5;

handles = cell(size(acc_all, 2), 1);

for feature_type = 1 : size(acc_all, 2)
   
    ax(feature_type) = axes; % Not possible to set multiple baselines within same axis;
    set(gca, 'XAxisLocation', 'top');
    
    if feature_type > 1
        ax(feature_type).YLim = ax(1).YLim;
        ax(feature_type).XLim = ax(1).XLim;
        hold on; % this locks the axes limits
    end
    
    % Create raincloud
    rain_scatter = (rand(size(acc_all{feature_type}, 1), 1) - 0.5) * rain_spread;
    handles{feature_type} = raincloud_plot(acc_all{feature_type}, 'box_on', 1, 'color', colours(feature_type, :), 'alpha', 0.8,...
        'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0,'density_type','rash');
    
    % Shift baseline (cloud)
    handles{feature_type}{1}.BaseValue = baselines(feature_type); % move baseline
    handles{feature_type}{1}.YData = handles{feature_type}{1}.YData + baselines(feature_type); % move cloud
    handles{feature_type}{1}.ShowBaseLine = 'off'; % hide baseline
    handles{feature_type}{1}.EdgeAlpha = 0; % hide cloud outline
    
    % Shift baseline (rain)
    handles{feature_type}{2}.YData = rain_scatter;
    handles{feature_type}{2}.YData = handles{feature_type}{2}.YData + baselines(feature_type) - rain_offset; % move rain
    handles{feature_type}{2}.SizeData = 3; % raindrop size
    
    % Shift boxplots
    handles{feature_type}{3}.Position([2 4]) = [-rain_spread/2 rain_spread]; % set box width to match rain spread
    handles{feature_type}{4}.YData = [-rain_spread/2 rain_spread/2]; % set median line to match rain spread
    handles{feature_type}{3}.Position(2) = handles{feature_type}{3}.Position(2) + baselines(feature_type) - rain_offset; % move box
    handles{feature_type}{4}.YData = handles{feature_type}{4}.YData + baselines(feature_type) - rain_offset; % move median line
    handles{feature_type}{5}.YData = handles{feature_type}{5}.YData + baselines(feature_type) - rain_offset; % move top whisker
    handles{feature_type}{6}.YData = handles{feature_type}{6}.YData + baselines(feature_type) - rain_offset; % move bot whisker
    
    % Hide all axes except the first
    if feature_type > 1
        axis off; % hide all axes except the first
    end
    
    view([-90 90]);
    
end

% Set same axis limits for all axes

% for l = 1:size(ax,2)
%     limits(l,:) = ax(l).YLim;
% end
% 
% for i = 1:size(ax,2)
%     ax(i).YLim = [min(limits(:)) max(limits(:))];
% end
% 
% linkaxes(ax);

limits = cell2mat(get(ax,'YLim')); % get both axes limits
set(ax,'YLim',[min(limits(:)) max(limits(:))]); % set the same values for both axes
linkaxes(ax);

% Link view angle, axes limits
%linkprop(ax, {'CameraPosition','CameraUpVector','YLim','XLim'});

%view([-90 90]); % view([90 -90]); % Swap x and y axes

% % Plot chance line
% plot([0.5 0.5], [-200 200], 'k--');
% 
% Formatting

%uncomment for RGB difference
% ylim(ax(1),[min(handles{end}{2}.YData)-3.*10^-7 max(handles{1}{1}.YData+3.*10^-7)]);

%for Luminance
ylim(ax(1),[min(handles{end}{2}.YData)-0.5.*10^-4 max(handles{1}{1}.YData)+0.5.*10^-4]);
title(ax(1), ['Eccentricity = ' num2str(unique(Results.Eccentricity)) ', n = ' num2str(size(Results,1))]);
xlabel(ax(1), 'Contrast difference');
ylabel(ax(1), 'Decision x confidence');
set(ax(1),'FontName','Arial','FontSize',18);
set(ax(1), 'YTick', fliplr(0:-base_offset:-7.*base_offset), 'YTickLabel', fliplr({'-4', '-3', '-2', '-1', '1', '2', '3','4'}));

%uncomment for plotting individual data points
% hold on
% plot(exception,[-3.*base_offset -3.*base_offset + 2.5.*10^-7],'o','MarkerFaceColor',colours(8,:),'MarkerSize',5,'MarkerEdgeColor','black');
% hold off
clear ax





%% for baselines
% decision = -3.5, spread = 1.5 x 10^-7, rainwidth = 0.5 x 10-7
% decision = -2.5, spread = 1.25 x 10^-7
% decision = -1.5, spread = 2.2 x 10^-7
% -0.5, 1.4
% 0.5, 1.4
% 1.5, 2.2
% 3.5,2.5

figure, raincloud_plot(Results.Contrast(Results.Decision == -3.5), 'box_on', 1, 'color', colours(2, :), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    
    
    
    
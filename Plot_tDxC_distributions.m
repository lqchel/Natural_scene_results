addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results'));
for exp = 1:2
Results = get_data2(exp);
Results(:,9) = Results(:,8).*Results(:,9);
Results(:,9) = Results(:,9) - 0.5;
Results(:,9) = Results(:,8).*Results(:,9);
img_id = unique(Results(:,11));

% signal present and absent trial classification
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1; %% Results(:,6) == 1 for exp 2
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;

%% original vs. modified: tDxC
% s = zeros(length(img_id),1);
% for p = 1:max(img_id)
%     s(p,:) = sum(Results(:,11) == p & Results(:,6) == 1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% congruent tDxC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_1 = zeros(length(img_id),6);
%size = zeros(80,1);
for img = 1:length(img_id)
%     current_cong= imread(fullfile(folder1,theFiles1(img-2).name));
%     current_incong= imread(fullfile(folder2,theFiles2(img-2).name));
    col_length = min([length(Results(Results(:,11)==img_id(img) & Find_Congruent_CP,9)) length(Results(Results(:,11)==img_id(img) & Find_Congruent_IP,9))]);
    dxc_c_c = Results(Results(:,11)==img_id(img) & Find_Congruent_CP,9);
    dxc_c_i = Results(Results(:,11)==img_id(img) & Find_Congruent_IP,9);
    temp = dxc_c_c(1:col_length,:) - dxc_c_i(1:col_length,:);
    delta_1(img,1) = mean(temp);
    delta_1(img,2) = img;
    delta_1(img,3) = col_length;
    [h,delta_1(img,4)] = ttest(temp,0);
    delta_1(img,5) = unique(Results(Results(:,11)==img_id(img) & Find_Congruent_CP,13));
    delta_1(img,6) = unique(Results(Results(:,11)==img_id(img) & Find_Congruent_CP,7));
    clear temp
    clear col_length
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% incongruent tDxC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
delta_2 = zeros(length(img_id),6);
%condition_diff = zeros(80,1);
for img = 1:length(img_id)
    col_length = min([length(Results(Results(:,11)==img_id(img) & Find_Incongruent_IP,9)) length(Results(Results(:,11)==img_id(img) & Find_Incongruent_CP,9))]);
    dxc_i_i = Results(Results(:,11)==img_id(img) & Find_Incongruent_IP,9);
    dxc_i_c = Results(Results(:,11)==img_id(img) & Find_Incongruent_CP,9);
    temp = dxc_i_i(1:col_length,:)-dxc_i_c(1:col_length,:);
    delta_2(img,1) = mean(temp);
    delta_2(img,2) = img_id(img);
    delta_2(img,3) = col_length;
    %condition_diff(img-2,1) = delta_1(img-2,1)-delta_2(img-2,1);
    [h,delta_2(img,4)] = ttest(temp,0);
    delta_2(img,5) = unique(Results(Results(:,11)==img_id(img) & Find_Incongruent_IP,13));
    delta_2(img,6) = unique(Results(Results(:,11)==img_id(img) & Find_Incongruent_IP,7));
    clear temp
end

%%%%%%%%%%%%%%%%%% tDxC tables for each condition %%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

cong_delta = table(delta_1(:,1),delta_1(:,2),delta_1(:,3),delta_1(:,4),delta_1(:,5),delta_1(:,6),'VariableNames',{'delta','img','num_sub','significance','eccentricity','location'});
incong_delta = table(delta_2(:,1),delta_2(:,2),delta_2(:,3),delta_2(:,4),delta_1(:,5),delta_1(:,6),'VariableNames',{'delta','img','num_sub','significance','eccentricity','location'});

%% compare delta between conditions

[r1,p1] = corrcoef(cong_delta.delta, incong_delta.delta)

condition_significance = zeros(80,1);
for img = 1:length(img_id)
    col_length = min([length(Results(Results(:,11)==img_id(img) & Find_Congruent_CP,9)) length(Results(Results(:,11)==img_id(img) & Find_Congruent_IP,9))]);
    dxc_c_c = Results(Results(:,11)==img_id(img) & Find_Congruent_CP,9);
    dxc_c_i = Results(Results(:,11)==img_id(img) & Find_Congruent_IP,9);
    temp1 = dxc_c_c(1:col_length,:) - dxc_c_i(1:col_length,:);
    clear col_length
    
    col_length = min([length(Results(Results(:,11)==img_id(img) & Find_Incongruent_IP,9)) length(Results(Results(:,11)==img_id(img) & Find_Incongruent_CP,9))]);
    dxc_i_i = Results(Results(:,11)==img_id(img) & Find_Incongruent_IP,9);
    dxc_i_c = Results(Results(:,11)==img_id(img) & Find_Incongruent_CP,9);
    temp2 = dxc_i_i(1:col_length,:)-dxc_i_c(1:col_length,:);
    
    [h,condition_significance(img,1)]= ttest2(temp1,temp2);
    clear temp1
    clear temp2
end

c_n_i_n = cong_delta.significance > 0.05 & incong_delta.significance > 0.05;
c_n_i_y = cong_delta.significance > 0.05 & incong_delta.significance < 0.05;
c_y_i_n = cong_delta.significance < 0.05 & incong_delta.significance > 0.05;
c_y_i_y = cong_delta.significance < 0.05 & incong_delta.significance < 0.05;


%% create scatterplot, by cong and incong conditions

c_n_i_n = cong_delta.significance > 0.05 & incong_delta.significance > 0.05;
c_n_i_y = cong_delta.significance > 0.05 & incong_delta.significance < 0.05;
c_y_i_n = cong_delta.significance < 0.05 & incong_delta.significance > 0.05;
c_y_i_y = cong_delta.significance < 0.05 & incong_delta.significance < 0.05;
x_lim = [-2 6];
y_lim = [-4 6];

temp_colour = cbrewer('qual', 'Set1', 9);
line_colour = cbrewer('qual','Set2',8);
sz = 50;
subplot(2,2,exp.*2-1),plot([0 0],y_lim,'--','Color',line_colour(3,:),'LineWidth',1.5); % add x = 0 line
hold on
plot(x_lim,[0 0],'--','Color',line_colour(3,:),'LineWidth',1.5); % y = 0
plot(x_lim,x_lim,'--','Color',line_colour(3,:),'LineWidth',1.5); % x = y
scatter(cong_delta.delta(c_n_i_n,1),incong_delta.delta(c_n_i_n,1),[],temp_colour(9,:),'filled','LineWidth',1.2)
scatter(cong_delta.delta(c_n_i_y,1),incong_delta.delta(c_n_i_y,1),[],'red','filled','LineWidth',1.2);
scatter(cong_delta.delta(c_y_i_n,1),incong_delta.delta(c_y_i_n,1),[],'blue','filled','LineWidth',1.2);
scatter(cong_delta.delta(c_y_i_y,1),incong_delta.delta(c_y_i_y,1),[],'black','filled','LineWidth',1.2);
scatter(cong_delta.delta(condition_significance < 0.05,1),incong_delta.delta(condition_significance < 0.05,1),95,'ks','LineWidth',0.8);
hold off

xlabel(['\Delta' 'tDxC in congruent images']), ylabel(['\Delta' 'tDxC in incongruent images']);
set(gca,'FontName','Arial','FontSize',14,'Box','off');
xlim(x_lim), ylim(y_lim);
[r,p] = corrcoef(cong_delta.delta(~c_n_i_n,1), incong_delta.delta(~c_n_i_n,1))
[r,p] = corrcoef(cong_delta.delta, incong_delta.delta);

%% plot results plotting the lines
colours = cbrewer('qual', 'Set1', 8);
[Y1,edges] = histcounts(delta_1(:,1),max(img_id));
Y1_cumulative  = cumsum(Y1);
subplot(2,2,exp.*2),plot(edges,[Y1_cumulative max(img_id)],'LineWidth',1.2,'Color',colours(2,:));

hold on
[Y2,edges] = histcounts(delta_2(:,1),max(img_id));
Y2_cumulative  = cumsum(Y2);
plot(edges,[Y2_cumulative max(img_id)],'LineWidth',1.2,'Color',colours(1,:));
hold off
box off

xlabel(['\Delta','tD×C (Original - Modified)']), ylabel('Cumulative counts');
set(gca,'FontName','Arial','FontSize',14);
xlim([-7 7]);
ylim([0 length(img_id)])
legend({'Congruent','Incongruent'},'Box','off','Location','northwest');

clear delta_1
clear delta_2
clear cong_delta
clear incong_delta
clear Results

% %% create scatterplot, by stimulus locations
% 
% plot([0 0],y_lim,'--','Color',line_colour(3,:),'LineWidth',1.5);
% hold on
% plot(x_lim,[0 0],'--','Color',line_colour(3,:),'LineWidth',1.5);
% plot(x_lim,x_lim,'--','Color',line_colour(3,:),'LineWidth',1.5);
% scatter(cong_delta.delta(cong_delta.eccentricity == 0),incong_delta.delta(incong_delta.eccentricity == 0),[],'red','filled','LineWidth',1.2);
% scatter(cong_delta.delta(cong_delta.eccentricity == 1),incong_delta.delta(incong_delta.eccentricity == 1),[],'blue','filled','LineWidth',1.2);
% scatter(cong_delta.delta(cong_delta.eccentricity == 2),incong_delta.delta(incong_delta.eccentricity == 2),[],'black','filled','LineWidth',1.2);
% scatter(cong_delta.delta(condition_significance < 0.05,1),incong_delta.delta(condition_significance < 0.05,1),95,'ks','LineWidth',0.8);
% hold off
% xlabel(['\Delta' 'tDxC in congruent images']), ylabel(['\Delta' 'tDxC in incongruent images']);
% set(gca,'FontName','Arial','FontSize',16,'Box','off');
% xlim(x_lim), ylim(y_lim)
% 
% [r,p] = corrcoef([cong_delta.eccentricity; incong_delta.eccentricity],[cong_delta.delta; incong_delta.delta])
end
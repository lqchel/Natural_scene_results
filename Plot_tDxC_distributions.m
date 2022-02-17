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
index = {Find_Congruent_CP; Find_Congruent_IP; Find_Incongruent_IP; Find_Incongruent_CP};
for img = 1:length(img_id)
    for cong = 1:2
        dxc_org = Results(Results(:,11)==img_id(img) & index{2.*cong-1},9);
        dxc_mod = Results(Results(:,11)==img_id(img) & index{2.*cong},9);
        img_num_sub = length(unique(Results(Results(:,11)==img_id(img) & index{2.*cong-1},1)));
        % delta tDxC, image id, number of participants, difference between
        % original and modified patch tDxC, eccentricity level, and location
        delta_mt{cong}(img,:) = [mean(dxc_org) - mean(dxc_mod) img_id(img) img_num_sub ttest2(dxc_org,dxc_mod)...
            unique(Results(Results(:,11)==img_id(img) & index{2.*cong-1},13))...
            unique(Results(Results(:,11)==img_id(img) & index{2.*cong-1},7))]; 
    end
end

%%%%%%%%%%%%%%%%%% tDxC tables for each condition %%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

cong_delta = array2table(delta_mt{1},'VariableNames',{'delta','img','num_sub','significance','eccentricity','location'});
incong_delta = array2table(delta_mt{2},'VariableNames',{'delta','img','num_sub','significance','eccentricity','location'});

%% determine color scale for scatter plot
c_n_i_n = cong_delta.significance == 0 & incong_delta.significance == 0;
c_n_i_y = cong_delta.significance == 0 & incong_delta.significance == 1;
c_y_i_n = cong_delta.significance == 1 & incong_delta.significance == 0;
c_y_i_y = cong_delta.significance == 1 & incong_delta.significance == 1;

%% compare delta between conditions

[r1,p1] = corrcoef(cong_delta.delta, incong_delta.delta)
org_or_mod = [0 1 0 1];
cong_incong = [0 0 1 1];
response = [];
patch_type = [];
congruence = [];
condition_significance = [];

for img = 1:length(img_id)
    for in = 1:size(index,1)
        if in == 1
            response = Results(Results(:,11)==img_id(img) & index{in},9);
            patch_type = zeros(length(Results(Results(:,11)==img_id(img) & index{in},9)),1) + org_or_mod(in);
            congruence = zeros(length(Results(Results(:,11)==img_id(img) & index{in},9)),1) + cong_incong(in);
        else
            response = [response; Results(Results(:,11)==img_id(img) & index{in},9)];
            patch_type = [patch_type; zeros(length(Results(Results(:,11)==img_id(img) & index{in},9)),1) + org_or_mod(in)];
            congruence = [congruence; zeros(length(Results(Results(:,11)==img_id(img) & index{in},9)),1) + cong_incong(in)];
        end        
    end
    anova_test = anovan(response,{categorical(patch_type),categorical(congruence)},'model',2,'varnames',{'patch_type','congruence'},'display','off');
    condition_significance(img,1) = anova_test(3,:); % take the interaction term only
end


%% create scatterplot, by cong and incong conditions

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
[Y1,edges] = histcounts(cong_delta.delta,max(img_id));
Y1_cumulative  = cumsum(Y1);
subplot(2,2,exp.*2),plot(edges,[Y1_cumulative max(img_id)],'LineWidth',1.2,'Color',colours(2,:));

hold on
[Y2,edges] = histcounts(incong_delta.delta,max(img_id));
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

%% display images based on tDxC
exp = 2;
% load image list
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\congruent cropped'; 
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\incongruent cropped';
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);

% load patch list
file_path_Congruent_Patch = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\congruent patch\';
filePattern3 = fullfile(file_path_Congruent_Patch,'*.jpg');
theFiles3 = dir(filePattern3);

% record critical object location in each image
for img = 1:max(img_id)
    Presentation_image_patch_name = sprintf('cong_%d_*.jpg',img);
    img_path_list = dir(strcat(file_path_Congruent_Patch,Presentation_image_patch_name));
    num_pch = length(img_path_list);  
    temp_array = ones(1,num_pch);
    for kk = 1:num_pch
        temp_array(kk) = img_path_list(kk).name(length(img_path_list(kk).name)-4);
    end
    c(img,1) = find(temp_array == 'p'); %%Find the CP position of the image
end

[delta_1_0,rank1] = sort(cong_delta.delta,'descend'); % sort image pairs based on delta tDxC of congruence images
delta_1_x = [delta_1_0 cong_delta.img(rank1,1) cong_delta.num_sub(rank1,1) c(rank1,1)];
delta_2_x = [incong_delta.delta(rank1,1) incong_delta.img(rank1,1) incong_delta.num_sub(rank1,1) c(rank1,1)];

cong_list =  {delta_1_x(c_n_i_n(rank1,1),:); delta_1_x(c_n_i_y(rank1,1),:);delta_1_x(c_y_i_n(rank1,1),:);delta_1_x(c_y_i_y(rank1,1),:)}; % categorize image pairs based on t-test significance
incong_list = {delta_2_x(c_n_i_n(rank1,1),:); delta_2_x(c_n_i_y(rank1,1),:);delta_2_x(c_y_i_n(rank1,1),:);delta_2_x(c_y_i_y(rank1,1),:)};

addpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results\cbrewer');
colours = cbrewer('qual','Pastel2',8);
figure_colours = [colours(8,:); colours(2,:); colours(3,:); 1 1 1];

%congruent
for category = 1:4
    current_cong_list = cong_list{category,1};
    current_incong_list = incong_list{category,1};
    
    if length(current_cong_list(:,1)) <= 5
        big_img_num = 1;
    else
        big_img_num = round(length(current_cong_list(:,1)),-1)/5;       
    end    
    
for i = 1:big_img_num
    f = figure('Color',figure_colours(category,:));
    set(gcf,'InvertHardCopy','off');
for p = 1:5
    current_cong= imresize(imread(fullfile(folder1,theFiles1(current_cong_list((i-1).*5+p,2)).name)),[150 150]); % load and resize image to square shape
    ax1= axes('Position',[0.015+(p-1).*0.19 0.52 0.2 0.2]); % fix axes position in figure window
    image(ax1,current_cong); % present image
    
    % mark out the critical objects
    hold on
    if current_cong_list((i-1).*5+p,4) == 1
        % outline the critical object on the image in red, using the x and y
        % coordinates of the object area
        plot(ax1,[1 50 50 1 1],[50 50 1 1 50],'r-','LineWidth',1.2); 
    elseif current_cong_list((i-1).*5+p,4) == 2
        plot(ax1,[50 100 100 50 50],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 3
        plot(ax1,[100 150 150 100 100],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 4
        plot(ax1,[1 50 50 1 1],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 5
        plot(ax1,[50 100 100 50 50],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 6
        plot(ax1,[100 150 150 100 100],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 7
        plot(ax1,[1 50 50 1 0],[150 150 100 100 150],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 8
        plot(ax1,[50 100 100 50 50],[150 150 100 100 150],'r-','LineWidth',1.2);
    else
        plot(ax1,[100 150 150 100 100],[150 150 100 100 150],'r-','LineWidth',1.2);
    end
    hold off
    axis square
    
    if exp == 1 % show number of subjects allocated to image for Exp 1
        title({['\Delta','cong = ', num2str(round(current_cong_list((i-1).*5+p,1),1)),'[',num2str(current_cong_list((i-1).*5 + p,3)),']'],...
            ['\Delta','incong = ', num2str(round(current_incong_list((i-1).*5+p,1),1)),'[',num2str(current_incong_list((i-1).*5 + p,3)),']']});
    else
        title({['\Delta','cong = ',  num2str(round(current_cong_list((i-1).*5+p,1),1))],...
            ['\Delta','incong = ', num2str(round(current_incong_list((i-1).*5+p,1),1))]});
    end
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    current_incong= imresize(imread(fullfile(folder2,theFiles2(current_cong_list((i-1).*5+p,2)).name)),[150 150]);
    ax2= axes('Position',[0.015+(p-1).*0.19 0.3 0.2 0.2]);
    image(ax2,current_incong);
    
    % mark out the critical objects
    hold on
    if current_cong_list((i-1).*5+p,4) == 1
        plot(ax2,[1 50 50 1 1],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 2
        plot(ax2,[50 100 100 50 50],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 3
        plot(ax2,[100 150 150 100 100],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 4
        plot(ax2,[1 50 50 1 1],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 5
        plot(ax2,[50 100 100 50 50],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 6
        plot(ax2,[100 150 150 100 100],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 7
        plot(ax2,[1 50 50 1 1],[150 150 100 100 150],'r-','LineWidth',1.2);
    elseif current_cong_list((i-1).*5+p,4) == 8
        plot(ax2,[50 100 100 50 50],[150 150 100 100 150],'r-','LineWidth',1.2);
    else
        plot(ax2,[100 150 150 100 100],[150 150 100 100 150],'r-','LineWidth',1.2);
    end
    hold off
    axis square
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    if (i-1).*5+p+1 > length(current_cong_list(:,1))
        break
    end
end
    if category == 1 && i == 1
        id = i;
    else
        id = id+1;
    end
    if id < 10
        filename = ['exp_' num2str(exp) 'cong_0', num2str(id),'.jpg'];
    else
        filename = ['exp_' num2str(exp) 'cong_', num2str(id), '.jpg'];
    end
     saveas(gcf,filename);
% if big_img_num ~= 2
 close(f);
% end
end




end


end
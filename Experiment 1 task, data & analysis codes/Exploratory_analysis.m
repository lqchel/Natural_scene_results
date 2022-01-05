%% get data
% Results = importdata('Pooled Results 2.mat'); 
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results'));
Results = get_data2(2);
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


% load image list
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\congruent cropped'; 
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\incongruent cropped';
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);



%% original vs. modified: tDxC
% s = zeros(length(img_id),1);
% for p = 1:max(img_id)
%     s(p,:) = sum(Results(:,11) == p & Results(:,6) == 1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% congruent tDxC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_1 = zeros(length(img_id),3);
%size = zeros(80,1);
cong_significance = zeros(80,1);
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
delta_2 = zeros(length(img_id),3);
%condition_diff = zeros(80,1);
incong_significance = zeros(80,1);
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
%% calculate weibull statistics

%%%%%%%% weibull parameters for image contrast
b = 1;
c = 1;
for img = 1:max(img_id) 
    for condition = 1:2
        % load and resize image to square shape, convert to grayscale
        if condition == 1
            current_img= rgb2gray(imresize(imread(fullfile(folder1,theFiles1(img).name)),[150 150])); 
        else
            current_img= rgb2gray(imresize(imread(fullfile(folder2,theFiles2(img).name)),[150 150]));
        end
        gradient_mag = imgradient(current_img); 
        gradient_vector = reshape(gradient_mag,[150.*150 1]) + 0.0001; % add 0.0001 to the gradient magnitude matrix so they are all positive values
        pd = fitdist(gradient_vector,'weibull');
        dist_values(c,:) = [pd.A pd.B]; 
        c = c+1;
        scale_values(condition) = pd.A; % record beta parameters for current image pair
          
       clear current_img
       clear pd
       clear gradient_mag
       clear gradient_vector
    end
    
    weibull_stat(b,:) = [img scale_values];
    b = b+1;
    
end

data_w = table([cong_delta.delta; incong_delta.delta],[weibull_stat(:,2);weibull_stat(:,2)],categorical([zeros(length(img_id),1); ones(length(img_id),1)]),'VariableNames',...
    {'delta','weibull','congruence'});
lm1 = fitlme(data_w,'delta~weibull*congruence + (1|weibull) + (1|congruence)')
lm_w = fitlme(data_w,'delta ~ weibull:congruence + congruence + (1|weibull) + (1|congruence)')
lm_c = fitlme(data_w,'delta ~ weibull:congruence + weibull +  (1|weibull) + (1|congruence)')
lm_in = fitlme(data_w,'delta ~ weibull + congruence + (1|weibull) + (1|congruence)')

weibull_effect = compare(lm_w, lm1)
congruence_effect = compare(lm_c,lm1)
interaction_effect = compare(lm_in,lm1)

subplot(1,2,1)
dot_colours = cbrewer('qual','Set2',8);
line_colour = cbrewer('qual','Set1',8);
w_effects = fixedEffects(lm1);
scatter(weibull_stat(:,2),cong_delta.delta,'MarkerEdgeColor',dot_colours(3,:),'MarkerFaceColor',dot_colours(3,:),'MarkerFaceAlpha',0.7);
hold on
scatter(weibull_stat(:,2),incong_delta.delta,'MarkerEdgeColor',dot_colours(2,:),'MarkerFaceColor',dot_colours(2,:),'MarkerFaceAlpha',0.7);
plot([0 180],[w_effects(1,:) 180.*w_effects(2,:) + w_effects(1,:)],'Color',line_colour(2,:),'LineWidth',2.3);
plot([0 180],[w_effects(1,:)+w_effects(3,:) 180.*(w_effects(2,:) + w_effects(4,:)) + w_effects(1,:) + w_effects(3,:)],'Color',line_colour(1,:),'LineWidth',2.3);
hold off
set(gca,'FontName','Arial','FontSize',16);
ylim([-6 6]),xlim([0 180]);
xlabel('Weibull Scale Parameter Value');
ylabel('\DeltatDxC');


%%%% plot beta and gemma values of weibull distribution

% scatter(dist_values(:,1),dist_values(:,2));
% h4 = lsline;
% h4.Color = 'r';
% h4.LineWidth = 1.2;
% [r,p] = corrcoef(dist_values(:,1),dist_values(:,2));
% title(['r = ' num2str(round(r(1,2),2)) ', p < .0001']);
% set(gca,'FontName','Arial','FontSize',12);
% xlabel('Scale (Beta)');
% ylabel('Shape (Gamma)');

%% plot original color image, gradient magnitude image, histogram and fitted weibull distribution

% sel_img = [12 1 48 54 97];
sel_img = 5;
figure;
for i = 1:length(sel_img)
    current_colour_img =  imresize(imread(fullfile(folder1,theFiles1(sel_img(i)).name)),[150 150]); %%% congruent
    gradient_mag = imgradient(rgb2gray(current_colour_img));
    gradient_vector = reshape(gradient_mag,[150.*150 1]) + 0.0001; % add 0.0001 to the gradient magnitude matrix so they are all positive values
    pd = fitdist(gradient_vector,'weibull');
    subplot(1,2,1),imshow(uint8(gradient_mag));
    subplot(1,2,2);
    histogram(gradient_mag,256,'Normalization','pdf')
    hold on
    x = linspace(0,max(gradient_vector));
    plot(x,pdf(pd,x),'LineWidth',3);
    xlim([0 max(gradient_vector)]);
    hold off
    title(['\beta' ' = ' num2str(round(pd.A,2))]);
    ylabel('Probability');
    set(gca,'FontName','Arial','FontSize',14)
    axis square
end

figure;
for i = 1:length(sel_img)
    current_colour_img =  imresize(imread(fullfile(folder2,theFiles2(sel_img(i)).name)),[150 150]); %%% congruent
    gradient_mag = imgradient(rgb2gray(current_colour_img));
    gradient_vector = reshape(gradient_mag,[150.*150 1]) + 0.0001; % add 0.0001 to the gradient magnitude matrix so they are all positive values
    pd = fitdist(gradient_vector,'weibull');
    subplot(1,2,1),imshow(uint8(gradient_mag));
    subplot(1,2,2);
    histogram(gradient_mag,256,'Normalization','pdf')
    hold on
    x = linspace(0,max(gradient_vector));
    plot(x,pdf(pd,x),'LineWidth',3);
    xlim([0 max(gradient_vector)]);
    hold off
    title(['\beta' ' = ' num2str(round(pd.A,2))]);
    xlabel('Edge Strength'), ylabel('Probability');
    set(gca,'FontName','Arial','FontSize',14);
    axis square
end
%% calculate sailiency difference statistics for each image
%%%% require SaliencyToolbox 2.3 by Itti et al
%%%% download link: http://www.saliencytoolbox.net/doc/index.html
addpath('C:\Users\liang\Documents\Experiment Codes\SaliencyToolbox');
b = 1;
final_salmap = [];
saliency_stat = [];
loc_mat = [1 2 3; 4 5 6; 7 8 9];

for img = 1:length(img_id) 
    % get difference matrix between congruent and incongruent images, by
    % subtracting incongruent image from congruent image, and get the
    % absolute value - imshow this matrix will reveal the area of change
    
     for condition = 1:2
        % load and resize image to square shape, convert to grayscale
        if condition == 1
            load_img= imresize(imread(fullfile(folder1,theFiles1(img).name)),[150 150]); 
        else
            load_img= imresize(imread(fullfile(folder2,theFiles2(img).name)),[150 150]);
        end
        current_img = initializeImage(load_img);
        params = defaultSaliencyParams;
        salmap = makeSaliencyMap(current_img,params);
        bigMap = imresize(salmap.data,current_img.size(1:2));
%         bigMap(bigMap<0) = 0;
        final_salmap(:,:,condition) = 255*mat2gray(bigMap); %% here convert the saliency map to grayscale image
    end
%         
% %         clear load_img
% %         clear current_img
% %         clear params
% %         clear salmap
% %         clear bigMap
% %         clear final_salmap
% %         clear crit_object_area
% %         clear crit_object_saliency

    [x,y] = find(loc_mat == cong_delta.location(img,:));
    range_rows = (x-1).*50+1:x.*50;
    range_cols = (y-1).*50+1:y.*50;
    saliency_difference = final_salmap(:,:,1) - final_salmap(:,:,2);
%     mean_saliency = mean(final_salmap, 3);
    accepted_difference = saliency_difference(range_rows,range_cols);
%     accepted_mean = mean_saliency(x,y);
    saliency_stat(b,:) = [img sum(abs(accepted_difference),'all')];
    b = b+1;
    
end
data = table([cong_delta.delta; incong_delta.delta],[saliency_stat(:,2);saliency_stat(:,2)],categorical([ones(length(img_id),1); 1+ones(length(img_id),1)]),'VariableNames',...
    {'delta','saliency','congruence'});
lm2 = fitlme(data,'delta~saliency*congruence + (1|saliency) + (1|congruence)');
lm_s = fitlme(data,'delta~congruence + saliency:congruence + (1|saliency) + (1|congruence)');
lm_c = fitlme(data,'delta~saliency + saliency:congruence + (1|saliency) + (1|congruence)');
lm_int = fitlme(data,'delta~saliency + congruence + (1|saliency) + (1|congruence)')

sal_effect = compare(lm_s,lm2)
cong_effect = compare(lm_c,lm2)
int_effect = compare(lm_int,lm2)


subplot(1,2,2)
dot_colours = cbrewer('qual','Set2',8);
line_colour = cbrewer('qual','Set1',8);
s_effects = fixedEffects(lm2);
s_effects(2,:) = 0.00001;
x_max = max(saliency_stat(:,2))+10000;
scatter(saliency_stat(:,2),cong_delta.delta,'MarkerEdgeColor',dot_colours(3,:),'MarkerFaceColor',dot_colours(3,:),'MarkerFaceAlpha',0.7);
hold on
scatter(saliency_stat(:,2),incong_delta.delta,'MarkerEdgeColor',dot_colours(2,:),'MarkerFaceColor',dot_colours(2,:),'MarkerFaceAlpha',0.7);
plot([0 x_max],[s_effects(1,:) x_max.*s_effects(2,:) + s_effects(1,:)],'Color',line_colour(2,:),'LineWidth',2.3);
plot([0 x_max],[s_effects(1,:)+s_effects(3,:) x_max.*(s_effects(2,:) + s_effects(4,:)) + s_effects(1,:) + s_effects(3,:)],'Color',line_colour(1,:),'LineWidth',2.3);

hold off
set(gca,'FontName','Arial','FontSize',16);
ylim([-6 6]),xlim([0 x_max]);
xlabel('|\DeltaSaliency(Cong - Incong)|');
ylabel('\DeltatDxC');



%% plotting images and saliency
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\SaliencyToolbox'));
final_salmap = [];

for img = 5:5
images{1} = imresize(imread(fullfile(folder1,theFiles1(img).name)),[150 150]);
images{2} = imresize(imread(fullfile(folder2,theFiles2(img).name)),[150 150]);
images{3} = abs(images{1} - images{2});
images_title = {'Congruent','Incongruent','|\Delta(Cong - Incong)|'};
rgb_difference = sum(images{3},3);
rgb_diff_quantile = quantile(rgb_difference, [0.95 1],  'all');
[row1,column1] = find(rgb_difference > rgb_diff_quantile(1,:) & rgb_difference < rgb_diff_quantile(2,:));
object_boundary = boundary(row1,column1);

figure;
for condition = 1:3
% original + difference images
ax1= axes('Position',[0.015 + (condition-1).*0.3 0.52 0.3 0.3]); 
image(ax1,images{condition});
hold on
plot(ax1, column1(object_boundary),row1(object_boundary),'r-','LineWidth',1.2);
hold off
axis square

set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
xticks([]),yticks([]);
title(images_title{condition});
end

% get saliency map
final_salmap = [];
for congruence = 1:2
    current_img = initializeImage(images{congruence});
    params = defaultSaliencyParams;
    salmap = makeSaliencyMap(current_img,params);
    bigMap = imresize(salmap.data,current_img.size(1:2));
    final_salmap{congruence} = 255*mat2gray(bigMap); %% here convert the saliency map to grayscale image
end
sal_difference = abs(final_salmap{1}-final_salmap{2});
% sal_quantile = quantile(sal_difference, [0.95 1],  'all');
final_salmap{3} = sal_difference;

% % plot saliency map and difference image

% [row2,column2] = find(sal_difference > sal_quantile(1,:) & sal_difference < sal_quantile(2,:)); % find the area of difference, by returning the row and column indices of pixels within that area 
% sal_boundary = boundary(row2,column2);
stat = sum(sal_difference(row1,column1),'all');
disp(round(stat,0))
for condition = 1:3
    ax2 = axes('Position',[0.015+(condition-1).*0.3 0.1 0.3 0.3]);
    image(ax2,final_salmap{condition});
    colormap gray
%     if condition == 3
%         hold on
%         plot(ax2,column1(object_boundary),row1(object_boundary),'r-','LineWidth',1.2);
%         hold off
%     end
    axis square
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    title(images_title{condition});
end
% if condition == 3
%     image(final_salmap{3})
%     colormap gray
%     set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
% end
% axis square
end


%% calculate critical object sizes
b = 1;
for img = 1:length(img_id) 
    difference_mat = abs(imresize(imread(fullfile(folder1,theFiles1(img_id(img)).name)),[150 150]) - imresize(imread(fullfile(folder2,theFiles2(img_id(img)).name)),[150 150]));
    difference = sum(difference_mat,3);
    current_object_size = (sum(sum(difference > 100)))/(150^2);   
    object_size(b,:) = [img current_object_size];
    b = b+1;
    
    clear current_object_size
    clear difference_mat
end

data_o = table([cong_delta.delta; incong_delta.delta],[object_size(:,2);object_size(:,2)],categorical([zeros(length(img_id),1); ones(length(img_id),1)]),'VariableNames',...
    {'delta','size','congruence'});
lm3 = fitlme(data_o,'delta~size*congruence + (1|size) + (1|congruence)')
lm_o = fitlme(data_o,'delta ~ size:congruence + congruence + (1|size) + (1|congruence)')
lm_c = fitlme(data_o,'delta ~ size:congruence + size +  (1|size) + (1|congruence)')
lm_in = fitlme(data_o,'delta ~ size + congruence + (1|size) + (1|congruence)')

size_effect = compare(lm_o, lm3)
congruence_effect = compare(lm_c,lm3)
interaction_effect = compare(lm_in,lm3)

subplot(1,2,1)

% legend({'Congruent','Incongruent','Image Stat Effect, Cong','Image Stat Effect, Incong'},'box','off')


%% plotting object size examples
b = 1;
selected_img = [5 42 35];
%selected_img = 5;
image_titles = {'Congruent','Incongruent','RGB difference','Detected object area',['\DeltaRGB histogram']};

 for img = 1:length(selected_img)
    images = {imresize(imread(fullfile(folder1,theFiles1(selected_img(img)).name)),[150 150]) imresize(imread(fullfile(folder2,theFiles2(selected_img(img)).name)),[150 150])};
    images{3} = abs(images{1} - images{2});
    difference = sum(images{3},3);
    current_object_size = (sum(sum(difference > 100)))/(150^2);   
    object_size(b,:) = [img current_object_size];
    b = b+1;
    figure('Color','white');
    
    for i = 1:5
        if i < 4
            ax1= axes('Position',[0.015+(i-1).*0.16 0.49 0.2 0.2]); % fix axes position in figure window
            image(ax1,images{i});
            xticks([]),yticks([]);
            set(gca, 'XColor','white','YColor','white')
            box off;
        elseif i == 4
            ax1= axes('Position',[0.015+(i-1).*0.16 0.49 0.2 0.2]); % fix axes position in figure window
            [row1,column1] = find(difference>100);
            image(ax1,images{3});
            hold on
            scatter(column1,row1,3,'MarkerFaceColor','red','MarkerEdgeColor','red');
            hold off
            xticks([]),yticks([]);
            set(gca, 'XColor','white','YColor','white')
            box off
        else
            ax2 = axes('Position',[0.02+(i-1).*0.172 0.49 0.2 0.2]); % fix axes position in figure window)
            [Y1,edges] = histcounts(difference,150^2);
            Y1_cumulative  = cumsum(Y1)/150^2;
            plot(ax2,[Y1_cumulative 1],edges,'LineWidth',2);
            xlim([0.92 1]), ylim([0 750]);
            yticks([0 100 300 500 700]);
            
            if img == length(selected_img)
            xlabel('Percentage');
            end
            
            ylabel(['\Delta RGB']);
            hold on
            plot([0.92 1],[100 100], 'r--','LineWidth', 1.5);
            hold off
            text(0.94,650,['Size = ' num2str(round(object_size(img,2),2))],'FontSize',14)

        end
        
        if img == 1
            title(image_titles{i});
        end
        axis square
        set(gca,'FontName','Arial','FontSize',12);
    end
    
    clear current_object_size
    clear difference_mat
    clear Y1
    clear Y1_cumulative
    clear edges
    clear images
end



%%% uncomment for plotting mean RGB difference across image pairs
% Y1_cumulative  = cumsum(mean(Y1,1));
% mean_edges = mean(edges,1);
% figure('Color','white');
% plot(mean_edges,[Y1_cumulative 150^2],'LineWidth',1.2);
% xlabel('Mean RGB difference across images'), ylabel('Counts');
% set(gca,'FontName','Arial','FontSize',14);


%%
subplot(1,2,1),scatter(cong_delta.eccentricity, object_size(:,2),'MarkerEdgeColor','blue','MarkerFaceColor','blue','jitter','on','jitterAmount',0.2,'MarkerFaceAlpha',0.3);
xlim([-0.5 2.5]),xticks([0 1 2]);
set(gca,'XTickLabel',{'F','P','P-P'},'FontSize',14,'Box','off');
ylabel('Object size (proportion in image)')

subplot(1,2,2),scatter(cong_delta.eccentricity,cong_delta.delta,'filled','MarkerFaceAlpha',0.5,'jitter','on','jitterAmount',0.2);
hold on
scatter(incong_delta.eccentricity,incong_delta.delta,'filled','MarkerFaceAlpha',0.5,'jitter','on','jitterAmount',0.2);
hold off
xlim([-0.5 2.5]),xticks([0 1 2]);
set(gca,'XTickLabel',{'F','P','P-P'},'FontSize',14,'Box','off');
ylabel(['\Delta' 'tDxC']);
legend({'Congruent','Incongruent'},'Box','off');

[h,p] = corrcoef(cong_delta.eccentricity, object_size(:,2))

%% statistical analysis
big_data_matrix = [repmat([cong_delta.img saliency_stat(:,2) object_size(:,2)],[2,1]) [weibull_stat(:,2); weibull_stat(:,3)]...
    [cong_delta.delta; incong_delta.delta]];
data = array2table(big_data_matrix, 'VariableNames',{'image','saliency','size','weibull','delta'});
data = addvars(data,categorical([zeros(136,1); ones(136,1)]), categorical(repmat(cong_delta.eccentricity,[2,1])),'NewVariableNames',...
    {'congruence','eccentricity'});

% size
lm1 = fitlme(data, 'delta ~ size*congruence + (1|image)+(congruence|image)')
lm2 = fitlme(data, 'delta ~ size + congruence + (1|image) + (congruence|image)')
interaction = compare(lm2,lm1)

lm3 = fitlme(data, 'delta ~ size + size:congruence + (1|image) + (congruence|image)')
congruence_effect = compare(lm3,lm1)

lm4 = fitlme(data, 'delta ~ congruence + size:congruence + (1|image) + (congruence|image)')
size_effect = compare(lm4,lm1)

cong_data = table(cong_delta.delta, object_size(:,2),object_size(:,1),'VariableNames',{'delta','size','image'});
lm1_cong = fitlme(cong_data,'delta ~ size + (1|image)')
lm1_cong_reduced = fitlme(cong_data,'delta ~ 1 + (1|image)')
cong_size_effect = compare(lm1_cong_reduced, lm1_cong)

incong_data = table(incong_delta.delta, object_size(:,2),object_size(:,1),'VariableNames',{'delta','size','image'});
lm1_incong = fitlme(incong_data, 'delta ~ size + (1|image)')
lm1_incong_reduced = fitlme(incong_data, 'delta ~ 1 + (1|image)');
incong_size_effect = compare(lm1_incong_reduced, lm1_incong)


%
dot_colours = cbrewer('qual','Set2',8);
line_colour = cbrewer('qual','Set1',8);
o_effects = fixedEffects(lm1);
x_max = max(object_size(:,2)) + 0.01;
scatter(object_size(:,2),cong_delta.delta,'MarkerEdgeColor',dot_colours(3,:),'MarkerFaceColor',dot_colours(3,:),'MarkerFaceAlpha',0.7);
hold on
scatter(object_size(:,2),incong_delta.delta,'MarkerEdgeColor',dot_colours(2,:),'MarkerFaceColor',dot_colours(2,:),'MarkerFaceAlpha',0.7);
plot([0 x_max],[o_effects(1,:) x_max.*o_effects(2,:) + o_effects(1,:)],'Color',line_colour(2,:),'LineWidth',2.3);
plot([0 x_max],[o_effects(1,:)+o_effects(3,:) x_max.*(o_effects(2,:) + o_effects(4,:)) + o_effects(1,:) + o_effects(3,:)],'Color',line_colour(1,:),'LineWidth',2.3);
hold off
set(gca,'FontName','Arial','FontSize',16);
ylim([-6 6]),xlim([0 x_max]);
xlabel('Object Size (Proportion in Image)');
ylabel('\DeltatDxC');


% weibull
%%%%

lm6 = fitlme(data, 'delta ~ size*weibull + (1|image)')
lm7 = fitlme(data, 'delta ~ size + weibull + (1|image)')
sw_interaction = compare(lm7,lm6)

lm8 = fitlme(data, 'delta ~ size + (1|image)')
lm9 = fitlme(data, 'delta ~ weibull + (1|image)')
weibull_effect = compare(lm8,lm7)
size_effect = compare(lm9,lm7)

lm10 = fitlme(data, 'delta ~ size*saliency + (1|image)')
lm11 = fitlme(data, 'delta ~ size + saliency + (1|image)')
lm12 = fitlme(data, 'delta ~ saliency + (1|image)')
ss_interaction = compare(lm11,lm10)
sal_effect = compare(lm8,lm11)
size_effect = compare(lm12,lm11)

fitlme(data, 'delta ~ size + weibull + saliency')

%% display images based on tDxC
% load image list
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

[delta_1_0,rank1] = sort(delta_1(:,1),'descend');
delta_1_x = [delta_1_0 delta_1(rank1,2) delta_1(rank1,3) c(rank1,1)];
delta_2_x = [delta_2(rank1,1) delta_2(rank1,2) delta_2(rank1,3) c(rank1,1)];

cong_list =  {delta_1_x(c_n_i_n(rank1,1),:); delta_1_x(c_n_i_y(rank1,1),:);delta_1_x(c_y_i_n(rank1,1),:);delta_1_x(c_y_i_y(rank1,1),:)};
incong_list = {delta_2_x(c_n_i_n(rank1,1),:); delta_2_x(c_n_i_y(rank1,1),:);delta_2_x(c_y_i_n(rank1,1),:);delta_2_x(c_y_i_y(rank1,1),:)};

addpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results\cbrewer');
colours = cbrewer('qual','Pastel2',8);
figure_colours = [1 1 1; colours(2,:); colours(3,:); colours(8,:)];

%congruent
for category = 1:4
    current_cong_list = cong_list{category,1};
    current_incong_list = incong_list{category,1};
    if round(length(current_cong_list(:,1)),-1) > length(current_cong_list(:,1))
        big_img_num = round(length(current_cong_list(:,1)),-1)/5;
    else
        big_img_num = 1;
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
    title({['d cong = ', num2str(round(current_cong_list((i-1).*5 + p,1)),1),'[',num2str(current_cong_list((i-1).*5 + p,3)),']'],...
        ['d incong = ', num2str(round(current_incong_list((i-1).*5 + p,1),1)),'[',num2str(current_incong_list((i-1).*5 + p,3)),']']});
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
        filename = ['cong_0', num2str(id),'.jpg'];
    else
        filename = ['cong_', num2str(id), '.jpg'];
    end
    saveas(gcf,filename);
% if big_img_num ~= 2
 close(f);
% end
end




end


%%
b = 1;
for sub = 1:num_sub
current_subject = Results(Results(:,1) == subject_id(sub),:);
locations = unique(current_subject(:,13));
    for loc = 1:length(locations)
        selected_subset = current_subject(current_subject(:,13) == locations(loc),:);
        Find_N = selected_subset(:,5) ==1;
        Find_CAP = selected_subset(:,4) == 0 & selected_subset(:,5) == 2; 
        Find_IAP = selected_subset(:,4) == 1 & selected_subset(:,5) == 3;
        current_CR = sum(Find_N & selected_subset(:,8) == -1)/sum(Find_N);
        current_hit = sum((Find_CAP|Find_IAP) & selected_subset(:,8) == 1)/sum(Find_CAP|Find_IAP);
        mat(b,:) = [subject_id(sub) locations(loc) current_hit current_CR ];
        b = b+1;

    end

end


AUC_lme([repmat(mat(:,1),[2,1]) [zeros(length(mat),1); ones(length(mat),1)] repmat(mat(:,2),[2,1]) [mat(:,3);mat(:,4)]] ,2,2)
AUC_lme([mat(:,1) mat(:,2) mat(:,3)],1,1)
AUC_lme([mat(:,1) mat(:,2) mat(:,4)],1,1)

for i = 1:3
se_CR(i) = std(mat(mat(:,2)==locations(i),4))/sqrt(240);
mean_CR(i) = mean(mat(mat(:,2)==locations(i),4));
se_hit(i) = std(mat(mat(:,2)==locations(i),3))/sqrt(240);
mean_hit(i) = mean(mat(mat(:,2)==locations(i),3));
end

errorbar([0 6.5 9.2],mean_CR,se_CR,'.-','LineWidth',1.2,'Capsize',10);
hold on
errorbar([0 6.5 9.2],mean_hit,se_hit,'.-','LineWidth',1.2,'Capsize',10);
hold off
legend({'CR','Hit'})
ylabel('Proportion of Responses');
xticks([0 6.5 9.2]);
set(gca,'XTickLabel',{'F','P','P-P'},'FontSize',12,'Box','off');
ylim([0.5 1]);
xlabel('Patch location');

[h,p] = ttest(mat(mat(:,2) == locations(1),3), mat(mat(:,2) == locations(1),4))
[h,p] = ttest(mat(mat(:,2) == locations(2),3), mat(mat(:,2) == locations(2),4))
[h,p] = ttest(mat(mat(:,2) == locations(3),3), mat(mat(:,2) == locations(3),4))


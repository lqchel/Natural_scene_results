%% get data
% Results = importdata('Pooled Results 2.mat'); 
Results = importdata('Exp2_data.mat');
Results = Results(Results(:,11)~=0,:);
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

%%%% uncomment for exp 1 results
% Results(:,11) = Results(:,2);
% img_id = unique(Results(:,11));
% ecc_level_1 = Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8;
% ecc_level_2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 2;
% ecc_levels = zeros(length(Results),1)+ ecc_level_1 + ecc_level_2;
% Results = [Results ecc_levels];

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

    clear temp
end

%%%%%%%%%%%%%%%%%% tDxC tables for each condition %%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

cong_delta = table(delta_1(:,1),delta_1(:,2),delta_1(:,3),delta_1(:,4),delta_1(:,5),'VariableNames',{'delta','img','num_sub','significance','eccentricity'});
incong_delta = table(delta_2(:,1),delta_2(:,2),delta_2(:,3),delta_2(:,4),delta_1(:,5),'VariableNames',{'delta','img','num_sub','significance','eccentricity'});


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

% 
%         histogram(gradient_mag,256,'Normalization','pdf')
%         hold on
%         x = linspace(0,max(gradient_vector));
%         plot(x,pdf(pd,x),'LineWidth',3);
%         xlim([0 max(gradient_vector)]);
%         pause(0.2)
%         hold off
          
       clear current_img
       clear pd
       clear gradient_mag
       clear gradient_vector
    end
    
    weibull_stat(b,:) = [img scale_values];
    b = b+1;
    
end

subplot(2,3,1)
[r,p] = corrcoef(weibull_stat(:,2),cong_delta.delta);
scatter(weibull_stat(:,2),cong_delta.delta);
h1 = lsline;
h1.Color = 'r';
h1.LineWidth = 1.2;
title(['r = ' num2str(round(r(1,2),2)) ', p = ' num2str(round(p(1,2),2)) '*']);
set(gca,'FontName','Arial','FontSize',12);
ylim([-4 6]),xlim([0 180]);
xlabel('Weibull Scale Parameter');
ylabel(['\Delta' 'Cong']);

subplot(2,3,2)
[r,p] = corrcoef(weibull_stat(:,3),incong_delta.delta);
scatter(weibull_stat(:,3),incong_delta.delta)
h2 = lsline;
h2.Color = 'r';
h2.LineWidth = 1.2;
title(['r = ' num2str(round(r(1,2),2)) ', p = ' num2str(round(p(1,2),2))]);
set(gca,'FontName','Arial','FontSize',12);
ylim([-4 6]),xlim([0 180]);
xlabel('Weibull Scale Parameter');
ylabel(['\Delta' 'Incong']);

subplot(2,3,3)
scatter(mean(weibull_stat(:,2:3),2), mean([cong_delta.delta incong_delta.delta],2));
h3 = lsline;
h3.Color = 'r';
h3.LineWidth = 1.2;
[r,p] = corrcoef(mean(weibull_stat(:,2:3),2), mean([cong_delta.delta incong_delta.delta],2));
title(['r = ' num2str(round(r(1,2),2)) ', p = ' num2str(round(p(1,2),2))]);
set(gca,'FontName','Arial','FontSize',12);
ylim([-4 6]),xlim([0 180]);
xlabel('Weibull Scale Parameter');
ylabel(['mean' '(' '\Delta' 'Cong' '+' '\Delta' 'Incong' ')']);

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

sel_img = [12 1 48 54 97];
for i = 1:length(sel_img)
    current_colour_img =  imresize(imread(fullfile(folder1,theFiles1(sel_img(i)).name)),[150 150]); 
    gradient_mag = imgradient(rgb2gray(current_colour_img));
    gradient_vector = reshape(gradient_mag,[150.*150 1]) + 0.0001; % add 0.0001 to the gradient magnitude matrix so they are all positive values
    pd = fitdist(gradient_vector,'weibull');
    
    figure;
    subplot(1,3,1),imshow(current_colour_img);
    subplot(1,3,2),imshow(uint8(gradient_mag));
    subplot(1,3,3);
    histogram(gradient_mag,256,'Normalization','pdf')
    hold on
    x = linspace(0,max(gradient_vector));
    plot(x,pdf(pd,x),'LineWidth',3);
    xlim([0 max(gradient_vector)]);
    hold off
    title(['\beta' ' = ' num2str(round(pd.A,2)) ' \gamma' ' = ' num2str(round(pd.B,2))]);
    xlabel('Edge Strength'), ylabel('Probability');
    axis square
end

%% calculate critical object sizes

b = 1;
for img = 1:length(img_id) 

    difference_mat = abs(imresize(imread(fullfile(folder1,theFiles1(img).name)),[150 150]) - imresize(imread(fullfile(folder2,theFiles2(img).name)),[150 150]));
    difference_mat = sum(difference_mat,3);
    current_object_size = sum(sum(difference_mat > 0))/150^2;   
    object_size(b,:) = [img current_object_size];
    clear current_object_size
    clear difference_mat
    b = b+1;
end


% subplot(2,3,1)
% [r,p] = corrcoef(object_size(:,2),cong_delta.delta);
% scatter(object_size(:,2),cong_delta.delta);
% h1 = lsline;
% h1.Color = 'r';
% h1.LineWidth = 1.2;
% title(['Congruent, r = ' num2str(round(r(1,2),2)) ', p = ' num2str(round(p(1,2),2))]);
% set(gca,'FontName','Arial','FontSize',12);
% ylim([-4 6]),xlim([0 0.9]);
% xlabel('Critical object size (% in image)');
% ylabel('\Delta' 'tDxC');
% 
% subplot(2,3,2)
% [r2,p2] = corrcoef(object_size(:,2),incong_delta.delta);
% scatter(object_size(:,2),incong_delta.delta);
% h2 = lsline;
% h2.Color = 'r';
% h2.LineWidth = 1.2;
% title(['Incongruent, r = ' num2str(round(r2(1,2),2)) ', p = ' num2str(round(p2(1,2),2))]);
% set(gca,'FontName','Arial','FontSize',12);
% ylim([-4 6]),xlim([0 0.9]);
% xlabel('Critical object size (% in image)');
% ylabel('\Delta tDxC');
% 
% subplot(2,3,3)
% scatter([object_size(:,2); object_size(:,2)], [cong_delta.delta; incong_delta.delta]);
% h3 = lsline;
% h3.Color = 'r';
% h3.LineWidth = 1.2;
% [r3,p3] = corrcoef([object_size(:,2); object_size(:,2)], [cong_delta.delta; incong_delta.delta]);
% title(['All images, r = ' num2str(round(r3(1,2),2)) ', p = ' num2str(round(p3(1,2),2))]);
% set(gca,'FontName','Arial','FontSize',12);
% ylim([-4 6]),xlim([0 0.9]);
% xlabel('Critical object size (% in image)');
% ylabel('\Delta tDxC');

% 
% 
% 
%% calculate sailiency difference statistics for each image
%%%% require SaliencyToolbox 2.3 by Itti et al
%%%% download link: http://www.saliencytoolbox.net/doc/index.html

b = 1;
for img = 1:length(img_id) 
    % get difference matrix between congruent and incongruent images, by
    % subtracting incongruent image from congruent image, and get the
    % absolute value - imshow this matrix will reveal the area of change
    
     difference_mat = abs(imresize(imread(fullfile(folder1,theFiles1(img).name)),[150 150]) - imresize(imread(fullfile(folder2,theFiles2(img).name)),[150 150]));
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

    difference = abs(final_salmap(:,:,1)-final_salmap(:,:,2));
    range = quantile(difference, [0.95 1],  'all');
    accepted_difference = difference(difference>range(1,:) & difference < range(2,:));
    saliency_stat(b,:) = [img sum(accepted_difference)];
    b = b+1;
    
end
 
subplot(2,3,4)
[r,p] = corrcoef(saliency_stat(:,2),cong_delta.delta);
scatter(saliency_stat(:,2),cong_delta.delta);
h1 = lsline;
h1.Color = 'r';
h1.LineWidth = 1.2;
title(['r = ' num2str(round(r(1,2),2)) ', p < .0001***']);
set(gca,'FontName','Arial','FontSize',12);
ylim([-4 6]),xlim([0 max(saliency_stat(:,2))+ 10000]);
xlabel('|\DeltaSaliency(Cong - Incong)|');
ylabel(['\Delta' 'Cong']);

subplot(2,3,5)
[r2,p2] = corrcoef(saliency_stat(:,2),incong_delta.delta);
scatter(saliency_stat(:,2),incong_delta.delta);
h2 = lsline;
h2.Color = 'r';
h2.LineWidth = 1.2;
title(['r = ' num2str(round(r2(1,2),2)) ', p < .0001***']);
set(gca,'FontName','Arial','FontSize',12);
ylim([-4 6]),xlim([0 max(saliency_stat(:,2))+ 10000]);
xlabel('|\DeltaSaliency(Cong - Incong)|');
ylabel(['\Delta' 'Incong']);

subplot(2,3,6)
[r3,p3] = corrcoef(saliency_stat(:,2),mean([cong_delta.delta incong_delta.delta],2))
scatter(saliency_stat(:,2),mean([cong_delta.delta incong_delta.delta],2));
h3 = lsline;
h3.Color = 'r';
h3.LineWidth = 1.2;
title(['r = ' num2str(round(r3(1,2),2)) ', p < .0001***']);
set(gca,'FontName','Arial','FontSize',12);
ylim([-4 6]),xlim([0 max(saliency_stat(:,2))+ 10000]);
xlabel('|\DeltaSaliency(Cong - Incong)|');
ylabel(['mean' '(' '\Delta' 'Cong' '+' '\Delta' 'Incong' ')']);


%% plotting images and saliency
for img = 1:2
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
% if condition == 3
%     hold on
%     plot(ax1, column1(object_boundary),row1(object_boundary),'r-','LineWidth',1.2);
%     hold off
% end
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
sal_quantile = quantile(sal_difference, [0.95 1],  'all');
final_salmap{3} = sal_difference;

% plot saliency map and difference image
[row2,column2] = find(sal_difference > sal_quantile(1,:) & sal_difference < sal_quantile(2,:)); % find the area of difference, by returning the row and column indices of pixels within that area 
sal_boundary = boundary(row2,column2);
stat = sum(sal_difference(sal_difference > sal_quantile(1,:) & sal_difference < sal_quantile(2,:)));
disp(round(stat,0))
for condition = 1:3
    ax2 = axes('Position',[0.015+(condition-1).*0.3 0.1 0.3 0.3]);
    image(ax2,final_salmap{condition});
    colormap gray
    if condition == 3
        hold on
        plot(ax2,column2(sal_boundary),row2(sal_boundary),'r-','LineWidth',1.2);
        hold off
    end
    axis square
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    title(images_title{condition});
end
end
%% plot results plotting the lines
colours = cbrewer('qual', 'Set1', 8);
[Y1,edges] = histcounts(delta_1(:,1),max(img_id));
Y1_cumulative  = cumsum(Y1);
plot(edges,[Y1_cumulative max(img_id)],'LineWidth',1.2,'Color',colours(2,:));

hold on
[Y2,edges] = histcounts(delta_2(:,1),max(img_id));
Y2_cumulative  = cumsum(Y2);
plot(edges,[Y2_cumulative max(img_id)],'LineWidth',1.2,'Color',colours(1,:));
hold off
box off

xlabel(['\Delta','tD�C (Original - Modified)']), ylabel('Cumulative counts');
set(gca,'FontName','Arial','FontSize',12);
xlim([-7 7]);
legend({'Congruent','Incongruent'},'Box','off','Location','northwest');

%% display images based on delta rcBxC
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
    c(img) = find(temp_array == 'p'); %%Find the CP position of the image
end

[delta_1_x,rank1] = sort(delta_1(:,1),'descend');
% [delta_2_x,rank2] = sort(delta_2(:,1),'descend');


delta_1_x = [delta_1_x delta_1(rank1,2) delta_1(rank1,3)];
delta_2_x = [delta_2(rank1,1) delta_2(rank1,2) delta_2(rank1,3)];

delta_1_x = [delta_1_x(:,1:3) c(rank1,1)];
delta_2_x = [delta_2_x(:,1:3) c(rank1,1)];

% display (axis definition [left bottom width height])

%congruent
for i = 1:28
for p = 1:5
    current_cong= imresize(imread(fullfile(folder1,theFiles1(delta_1_x((i-1).*5+p,2)).name)),[150 150]); % load and resize image to square shape
    ax1= axes('Position',[0.015+(p-1).*0.19 0.52 0.2 0.2]); % fix axes position in figure window
    image(ax1,current_cong); % present image
    
    % mark out the critical objects
    hold on
    if delta_1_x((i-1).*5+p,4) == 1
        % outline the critical object on the image in red, using the x and y
        % coordinates of the object area
        plot(ax1,[1 50 50 1 1],[50 50 1 1 50],'r-','LineWidth',1.2); 
    elseif delta_1_x((i-1).*5+p,4) == 2
        plot(ax1,[50 100 100 50 50],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 3
        plot(ax1,[100 150 150 100 100],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 4
        plot(ax1,[1 50 50 1 1],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 5
        plot(ax1,[50 100 100 50 50],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 6
        plot(ax1,[100 150 150 100 100],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 7
        plot(ax1,[1 50 50 1 0],[150 150 100 100 150],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 8
        plot(ax1,[50 100 100 50 50],[150 150 100 100 150],'r-','LineWidth',1.2);
    else
        plot(ax1,[100 150 150 100 100],[150 150 100 100 150],'r-','LineWidth',1.2);
    end
    hold off
    axis square
    title({['d cong = ', num2str(round(delta_1_x((i-1).*5 + p,1)),1),'[',num2str(delta_1_x((i-1).*5 + p,3)),']'],...
        ['d incong = ', num2str(round(delta_2_x((i-1).*5 + p,1),1)),'[',num2str(delta_1_x((i-1).*5 + p,3)),']']});
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    
    current_incong= imresize(imread(fullfile(folder2,theFiles2(delta_1_x((i-1).*5+p,2)).name)),[150 150]);
    ax2= axes('Position',[0.015+(p-1).*0.19 0.3 0.2 0.2]);
    image(ax2,current_incong);
    
    % mark out the critical objects
    hold on
    if delta_1_x((i-1).*5+p,4) == 1
        plot(ax2,[1 50 50 1 1],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 2
        plot(ax2,[50 100 100 50 50],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 3
        plot(ax2,[100 150 150 100 100],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 4
        plot(ax2,[1 50 50 1 1],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 5
        plot(ax2,[50 100 100 50 50],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 6
        plot(ax2,[100 150 150 100 100],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 7
        plot(ax2,[1 50 50 1 1],[150 150 100 100 150],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 8
        plot(ax2,[50 100 100 50 50],[150 150 100 100 150],'r-','LineWidth',1.2);
    else
        plot(ax2,[100 150 150 100 100],[150 150 100 100 150],'r-','LineWidth',1.2);
    end
    hold off
    axis square
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    if i == 28 && p == 1
        break
    end
end

    if i < 10
        filename = ['cong_0', num2str(i),'.jpg'];
    else
        filename = ['cong_', num2str(i), '.jpg'];
    end
    saveas(gcf,filename);
clf;
end


%%
%incongruent
for i = 1:16

for p = 1:5
    current_cong= imresize(imread(fullfile(folder1,theFiles1(delta_2_x((i-1).*5+p,2)).name)),[150 150]);
    ax1= axes('Position',[0.015+(p-1).*0.19 0.52 0.2 0.2]);
    image(ax1,current_cong);
    axis square
    title({['delta = ', num2str(round(delta_2_x((i-1).*5 + p,1),1)), ', ', num2str(round(delta_2_x((i-1).*5 + p,3),1)), ' dva']});
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    current_incong= imresize(imread(fullfile(folder2,theFiles2(delta_2_x((i-1).*5+p,2)).name)),[150 150]);
    ax2= axes('Position',[0.015+(p-1).*0.19 0.3 0.2 0.2]);
    image(ax2,current_incong);
    axis square
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);

end

    if i < 10
        filename = ['incong_0', num2str(i),'.jpg'];
    else
        filename = ['incong_', num2str(i), '.jpg'];
    end
    saveas(gcf,filename);
    clf;
end

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
plot([0 0],y_lim,'--','Color',line_colour(3,:),'LineWidth',1.5); % add x = 0 line
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
set(gca,'FontName','Arial','FontSize',16,'Box','off');
xlim(x_lim), ylim(y_lim);
[r,p] = corrcoef(cong_delta.delta(~c_n_i_n,1), incong_delta.delta(~c_n_i_n,1))
[r,p] = corrcoef(cong_delta.delta, incong_delta.delta);

%% create scatterplot, by stimulus locations

plot([0 0],y_lim,'--','Color',line_colour(3,:),'LineWidth',1.5);
hold on
plot(x_lim,[0 0],'--','Color',line_colour(3,:),'LineWidth',1.5);
plot(x_lim,x_lim,'--','Color',line_colour(3,:),'LineWidth',1.5);
scatter(cong_delta.delta(cong_delta.eccentricity == 0),incong_delta.delta(incong_delta.eccentricity == 0),[],'red','filled','LineWidth',1.2);
scatter(cong_delta.delta(cong_delta.eccentricity == 1),incong_delta.delta(incong_delta.eccentricity == 1),[],'blue','filled','LineWidth',1.2);
scatter(cong_delta.delta(cong_delta.eccentricity == 2),incong_delta.delta(incong_delta.eccentricity == 2),[],'black','filled','LineWidth',1.2);
scatter(cong_delta.delta(condition_significance < 0.05,1),incong_delta.delta(condition_significance < 0.05,1),95,'ks','LineWidth',0.8);
hold off
xlabel(['\Delta' 'tDxC in congruent images']), ylabel(['\Delta' 'tDxC in incongruent images']);
set(gca,'FontName','Arial','FontSize',16,'Box','off');
xlim(x_lim), ylim(y_lim)

[r,p] = corrcoef([cong_delta.eccentricity; incong_delta.eccentricity],[cong_delta.delta; incong_delta.delta])


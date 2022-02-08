%% get data
% Results = importdata('Pooled Results 2.mat'); 
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results'));
Results = get_data2(2);
Results(:,9) = Results(:,8).*Results(:,9);
Results(:,9) = Results(:,9) - 0.5;
Results(:,9) = Results(:,8).*Results(:,9);
img_id = unique(Results(:,11));
subject_id = unique(Results(:,1));


% % %% uncomment for checking Exp 1 data, first 6 patches of each trial
% % Results = Results(Results(:,3)<= 6,:);

%%%% uncomment for checking Exp 1 data, first 24 trials of each
%%%% participant
% Results = Results(Results(:,2)<= 27,:);

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

%% check data

problem_data = zeros(1,14);
for sub = 1:length(subject_id)
    sub_data = Results(Results(:,1) == subject_id(sub),:);
    current_sub_img = unique(sub_data(:,11));
    for i = 1:length(current_sub_img)
        if sum(sub_data(:,11) == current_sub_img(i) & sub_data(:,6) == 1) == 0
            problem_data = [problem_data; sub_data(sub_data(:,11)== current_sub_img(i),:)];
        end
    end
end
problem_data(1,:) = [];
length(unique(problem_data(2,:)))

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

cong_delta_ = array2table(delta_mt{1},'VariableNames',{'delta','img','num_sub','significance','eccentricity','location'});
incong_delta_ = array2table(delta_mt{2},'VariableNames',{'delta','img','num_sub','significance','eccentricity','location'});

%% determine color scale for scatter plot
c_n_i_n = cong_delta.significance == 0 & incong_delta.significance == 0;
c_n_i_y = cong_delta.significance == 0 & incong_delta.significance == 1;
c_y_i_n = cong_delta.significance == 1 & incong_delta.significance == 0;
c_y_i_y = cong_delta.significance == 1 & incong_delta.significance == 1;

%% compare delta between conditions
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
    condition_significance(img) = anova_test(3,:); % take the interaction term only
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
          
       clear current_img
       clear pd
       clear gradient_mag
       clear gradient_vector
    end
    
    weibull_stat(b,:) = [img scale_values];
    b = b+1;
    
end

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
    figure;
    % plot original color image
    ax1= axes('Position',[0.015 0.49 0.27 0.27]);
    current_colour_img = imresize(imread(fullfile(folder1,theFiles1(sel_img(i)).name)),[150 150]); %%% congruent
    image(ax1,current_colour_img);
    xticks([]),yticks([]);
    set(gca, 'XColor','white','YColor','white')
    box off;
    axis square
    
    % compute and plot gradient magnitude image
    ax2= axes('Position',[0.28 0.49 0.27 0.27]);
    gradient_mag = imgradient(rgb2gray(current_colour_img));
    image(ax2,gradient_mag);
    xticks([]),yticks([]);
    set(gca, 'XColor','white','YColor','white') % get rid of x and y axes
    box off;
    colormap gray % show gradient image in black and white
    axis square


    % compute and plot contrast + weibull distributions
    ax3 = axes('Position',[0.63 0.49 0.27 0.27]);
    gradient_vector = reshape(gradient_mag,[150.*150 1]) + 0.0001; % reshape the gradient mag to a single column to make histogram, and add 0.0001 to all values so they are all positive values
    pd = fitdist(gradient_vector,'weibull');
    histogram(ax3,gradient_mag,256,'Normalization','pdf')
    hold on
    x = linspace(0,max(gradient_vector));
    plot(x,pdf(pd,x),'LineWidth',4);
    xlim([0 max(gradient_vector)]);
    hold off
    xlim([0 800]);
    a = gca;
    text(300,max(a.YLim).*0.8,['\beta' ' = ' num2str(round(pd.A,2))],'FontSize',18);
    ylabel('Probability');
    set(gca,'FontName','Arial','FontSize',14)
    axis square
    if i == length(sel_img)
        xlabel('Edge Strengh')
    end
end

% figure;
% for i = 1:length(sel_img)
%     current_colour_img =  imresize(imread(fullfile(folder2,theFiles2(sel_img(i)).name)),[150 150]); %%% congruent
%     gradient_mag = imgradient(rgb2gray(current_colour_img));
%     gradient_vector = reshape(gradient_mag,[150.*150 1]) + 0.0001; % add 0.0001 to the gradient magnitude matrix so they are all positive values
%     pd = fitdist(gradient_vector,'weibull');
%     subplot(1,2,1),imshow(uint8(gradient_mag));
%     subplot(1,2,2);
%     histogram(gradient_mag,256,'Normalization','pdf')
%     hold on
%     x = linspace(0,max(gradient_vector));
%     plot(x,pdf(pd,x),'LineWidth',3);
%     xlim([0 max(gradient_vector)]);
%     hold off
%     title(['\beta' ' = ' num2str(round(pd.A,2))]);
%     xlabel('Edge Strength'), ylabel('Probability');
%     set(gca,'FontName','Arial','FontSize',14);
%     axis square
% end
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



%% plotting images and saliency
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\SaliencyToolbox'));
final_salmap = [];
selected_img = [5 42 35];

for img = 1:length(selected_img)
images{1} = imresize(imread(fullfile(folder1,theFiles1(selected_img(img)).name)),[150 150]);
images{2} = imresize(imread(fullfile(folder2,theFiles2(selected_img(img)).name)),[150 150]);
object_boundary = boundary(row1,column1);

figure;
for condition = 1:2
    % original images
    ax1= axes('Position',[0.015+(condition-1).*0.16 0.49 0.2 0.2]); 
    image(ax1,images{condition});
    axis square
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
end

% get saliency map
final_salmap = [];
for congruence = 1:2
    current_img = initializeImage(images{congruence});
    params = defaultSaliencyParams;
    salmap = makeSaliencyMap(current_img,params);
    bigMap = imresize(salmap.data,current_img.size(1:2));
    final_salmap{congruence} = 255*mat2gray(bigMap); %% here convert the saliency map to grayscale image
    
    % show saliency maps of original and modified patches
    ax2 = axes('Position',[0.015+(congruence+1).*0.16 0.49 0.2 0.2]); 
    image(ax2,final_salmap{congruence});
    axis square
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    colormap gray
end

% get saliency difference map
sal_difference = abs(final_salmap{1}-final_salmap{2});
ax3 = axes('Position',[0.015+4.*0.16 0.49 0.2 0.2]);
image(ax3, sal_difference);
axis square
set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
xticks([]),yticks([]);
colormap gray

end


%% calculate critical object sizes
b = 1;
criteria = 10; % set RGB level selecion criteria
object_size = [];
current_object_size = [];

for img = 1:length(img_id) 
    difference_mat = abs(imresize(imread(fullfile(folder1,theFiles1(img_id(img)).name)),[150 150]) - imresize(imread(fullfile(folder2,theFiles2(img_id(img)).name)),[150 150]));
    difference = sum(difference_mat,3);
    current_object_size = (sum(sum(difference > criteria)))/(150^2);  
    object_size(b,:) = [img current_object_size log(current_object_size)];
    b = b+1;
    
    clear current_object_size
    clear difference_mat
end

% figure('Color','white');
% [Y1,edges_1] = histcounts(object_size(:,2),length(img_id));
% Y1_cumulative  = cumsum(Y1);
% subplot(1,2,1),plot(edges_1,[Y1_cumulative max(img_id)],'LineWidth',1.2);
% xlabel('Raw Object Size Data'), ylabel('Count');
% 
% [Y2,edges_2] = histcounts(object_size(:,3),length(img_id));
% Y2_cumulative = cumsum(Y2);
% subplot(1,2,2),plot(edges_2,[Y2_cumulative max(img_id)],'LineWidth',1.2);
% xlabel('Log-transformed Object Size Data'), ylabel('Count');

%% plotting object size examples
b = 1;
% selected_img = [5 35 12 45 107]; 
 %selected_img = 45; 
 selected_img = 5;
object_size = [];
image_titles = {'Congruent','Incongruent','RGB difference','Detected object area',['\DeltaRGB histogram']};
criteria = 10;


 for img = 1:length(selected_img)
    images = {imresize(imread(fullfile(folder1,theFiles1(selected_img(img)).name)),[150 150]) imresize(imread(fullfile(folder2,theFiles2(selected_img(img)).name)),[150 150])};
    images{3} = abs(images{1} - images{2});
    difference = sum(images{3},3);
    current_object_size = (sum(sum(difference > criteria)))/(150^2);   
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
            [row1,column1] = find(difference>criteria);
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
%              xlim([0.75 1]);
              xlim([0.97 1])
%             ylim([0 400]);
%             yticks(0:100:400);
            %%%% uncomment for normal criteria
            yticks([0 50 100])
            ylim([0 100])
            if img == length(selected_img)
            xlabel('Percentage');
            end
            
            ylabel(['\Delta RGB']);
            hold on
              plot([0.97 1],[criteria criteria], 'r--','LineWidth', 1.5);
%               plot([0.75 1],[criteria criteria], 'r--','LineWidth', 1.5);
            hold off
            text(0.975,80,['Size = ' num2str(round(object_size(img,2),3))], 'FontSize',14)
%              text(0.77,650,['Size = ' num2str(round(object_size(img,2),3))], 'FontSize',14)
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





%% statistical analysis
big_data_matrix = [repmat([cong_delta.img saliency_stat(:,2) object_size(:,3)],[2,1]) [weibull_stat(:,2); weibull_stat(:,3)]...
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

cong_data = table(cong_delta.delta, object_size(:,3),object_size(:,1),'VariableNames',{'delta','size','image'});
lm1_cong = fitlme(cong_data,'delta ~ size + (1|image)')
lm1_cong_reduced = fitlme(cong_data,'delta ~ 1 + (1|image)')
cong_size_effect = compare(lm1_cong_reduced, lm1_cong)

incong_data = table(incong_delta.delta, object_size(:,3),object_size(:,1),'VariableNames',{'delta','size','image'});
lm1_incong = fitlme(incong_data, 'delta ~ size + (1|image)')
lm1_incong_reduced = fitlme(incong_data, 'delta ~ 1 + (1|image)');
incong_size_effect = compare(lm1_incong_reduced, lm1_incong)


%
figure('Color','white');
dot_colours = cbrewer('qual','Set2',8);
line_colour = cbrewer('qual','Set1',8);
o_effects_cong = fixedEffects(lm1_cong);
o_effects_incong = fixedEffects(lm1_incong);
x_min = min(object_size(:,3));
x_max = max(object_size(:,3))+0.5;
scatter(object_size(:,3),cong_delta.delta,'MarkerEdgeColor',dot_colours(3,:),'MarkerFaceColor',dot_colours(3,:),'MarkerFaceAlpha',0.7);
hold on
scatter(object_size(:,3),incong_delta.delta,'MarkerEdgeColor',dot_colours(2,:),'MarkerFaceColor',dot_colours(2,:),'MarkerFaceAlpha',0.7);
plot([x_min x_max],[x_min.*o_effects_cong(2,:) + o_effects_cong(1,:) x_max.*o_effects_cong(2,:) + o_effects_cong(1,:)],'Color',line_colour(2,:),'LineWidth',2.3);
plot([x_min x_max],[x_min.*o_effects_incong(2,:) + o_effects_incong(1,:) x_max.*o_effects_incong(2,:) + o_effects_incong(1,:)],'Color',line_colour(1,:),'LineWidth',2.3);
hold off
set(gca,'FontName','Arial','FontSize',16);
ylim([-6 6]),xlim([x_min x_max]);
xlabel('Log-scale Object Size');
ylabel('\DeltatDxC');


% weibull
%%%%

lm6 = fitlme(data, 'delta ~ size*weibull+ (1|image)')
lm7 = fitlme(data, 'delta ~ size + weibull+ (1|image)')
sw_interaction = compare(lm7,lm6)

lm8 = fitlme(data, 'delta ~ size + (1|image)')
lm9 = fitlme(data, 'delta ~ weibull + (1|image)')
weibull_effect = compare(lm8,lm7)
size_effect = compare(lm9,lm7)

% saliency
lm10 = fitlme(data, 'delta ~ size*saliency+ (1|image)')
lm11 = fitlme(data, 'delta ~ size + saliency+ (1|image)')
lm12 = fitlme(data, 'delta ~ saliency + (1|image)')
ss_interaction = compare(lm11,lm10)
sal_effect = compare(lm8,lm11)
size_effect = compare(lm12,lm11)

fitlme(data, 'delta ~ size + weibull + saliency')

% eccentricity 
lm13 = fitlme(data, 'delta ~ size*eccentricity + (1|image)')

%% plot tDxC histogram
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

%% delta tDxC against eccentricity
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

[delta_1_0,rank1] = sort(cong_delta.delta(:,1),'descend'); % sort image pairs based on delta tDxC of congruence images
delta_1_x = [delta_1_0 delta_1(rank1,2) delta_1(rank1,3) c(rank1,1)];
delta_2_x = [delta_2(rank1,1) delta_2(rank1,2) delta_2(rank1,3) c(rank1,1)];

cong_list =  {delta_1_x(c_n_i_n(rank1,1),:); delta_1_x(c_n_i_y(rank1,1),:);delta_1_x(c_y_i_n(rank1,1),:);delta_1_x(c_y_i_y(rank1,1),:)}; % categorize image pairs based on t-test significance
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




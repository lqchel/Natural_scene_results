% this code retrieves the image viewing sizes of online participants, based
% on the data recorded from the calibration procedure.

%% obtain image viewing size for each subject
clear all
batch = 1:4;
group = 1:4;
a = 1;

for b = 1:length(batch)
    for g = 1:length(group)
        folder = ['C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Exp2Results\formal data\b'...
            num2str(batch(b)) '_g' num2str(group(g))];
        data_file_path = fullfile(folder,'*.xlsx');
        theFiles = struct2cell(dir(data_file_path));
        data_file_name(:,1:size(theFiles,2)) = string(theFiles(1,:)); % read data files
        within_group_sub_code = [10 11 12 13 14 15 1 2 3 4 5 6 7 8 9]; % retrieve subject number
        subject_code = (batch(b)-1).*60 + (group(g)-1).*15 + within_group_sub_code;
        
            for i = 1:length(data_file_name)
                sprintf('Reading batch = %d, group = %d, subject = %d', b,g, within_group_sub_code(i))
                % get screen height, which retries the image size in pixels
                % : image size = 75% of screen height
                screen_y = cell2mat(table2array(readtable(data_file_name(:,i),'Range','O250',...
                    'ReadVariableNames',false)));
                screen_y = string(screen_y(1:3));
                screen_y = str2num(screen_y);
                % get viewing distance
                view_distance = table2array(readtable(data_file_name(:,i),'Range','P250',... 
                    'ReadVariableNames',false));
                % get image viewing size in degrees of visual angle:
                % viewing size in radian = 2 x atan(image size in mm/2
                % viewing distance in mm); 1 radian = ~57.30 degrees
                img_size = (2.*atan((0.75.*screen_y)/(view_distance.*2))).*57.3;
                all_view_sizes(a,:) = [b g subject_code(i) img_size];
                a = a+1;
                clear screen_y
                clear view_distance
                clear img_size
            end
    end    
end

[code_sorted, I] = sort(all_view_sizes(:,3));
all_view_sizes_sorted = all_view_sizes(I,:);
save('image_viewing_sizes.mat','all_view_sizes_sorted','-mat');


%% plot cumulative histogram
colours = cbrewer('qual', 'Set1', 8);
[Y1,edges] = histcounts(all_view_sizes(:,4),240);
Y1_cumulative  = cumsum(Y1);
plot(edges,[Y1_cumulative 240],'LineWidth',1.2,'Color',colours(2,:));
ylim([0 240]), yticks([0:40:240]);
xlabel('Image viewing size (dva)'), ylabel('Cumulative counts');
set(gca,'FontName','Arial','FontSize',14);
title(['M = ' num2str(round(mean(all_view_sizes_sorted(:,4)),1)) ', SD = ' num2str(round(std(all_view_sizes_sorted(:,4)),1))]);

%% calculate patch eccentricity for each individual participant, record it in existing data matrix
load('image_viewing_sizes.mat');
data = importdata('Exp2_data.mat');

for sub = 1:length(unique(data(:,1)))
    current_view_size = all_view_sizes_sorted(all_view_sizes_sorted(:,3)== sub,4);
    Results = data(data(:,1) == sub,:);
    % identify fovea, periphery, and para-periphery locations as 0, 1 & 2
    ecc_level_1 = Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8;
    ecc_level_2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 2;
    current_ecc_levels = zeros(length(Results),1)+ ecc_level_1 + ecc_level_2;
    
    % calculate actual eccentricity for each patch locations
    ecc_scale = current_view_size/3;
    actual_ecc_1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* ecc_scale;
    actual_ecc_2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* sqrt(2.*ecc_scale^2);
    current_actual_ecc = zeros(length(Results),1) + actual_ecc_1 + actual_ecc_2;
    
    if sub == 1
        ecc_levels = current_ecc_levels;
        actual_ecc = current_actual_ecc;
    else
        ecc_levels = [ecc_levels; current_ecc_levels];
        actual_ecc = [actual_ecc; current_actual_ecc];
    end
    
    clear Results
    clear ecc_level_1
    clear ecc_level_2
    clear current_ecc_levels
    clear actual_ecc_1
    clear actual_ecc_2
    clear current_actual_ecc
end

data = [data ecc_levels actual_ecc];
save('Exp2_data.mat','data','-mat');
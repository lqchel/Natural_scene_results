% this code finds the null patches that match with the modified patch in
% low-level property difference

%% load data

load('All absent decision.mat');
load('all present patch name.mat');
load('all absent patch name.mat');

select_num = [7 12 4 15 20] - 2;


%% for each of the 5 images, select null patches with least difference on R/G/B layer

for img = 1:5
    % patch responses for current image
    modified_data = data_allab(data_allab.Image == select_num(img) & data_allab.PatchType ~=0 & data_allab.Congruency == 0,:);
    null_data = data_allab(data_allab.Image == select_num(img) & data_allab.PatchType ==0 & data_allab.Congruency == 0 & ...
        data_allab.AbsentLoc == unique(modified_data.AbsentLoc),:);
    
    % patch file names
    null_list = absent_list(data_allab.Image == select_num(img) & data_allab.PatchType ==0 & data_allab.Congruency == 0 & ...
        data_allab.AbsentLoc == unique(modified_data.AbsentLoc),:); 
    original_list = present_list(data_allab.Image == select_num(img) & data_allab.PatchType ==0 & data_allab.Congruency == 0 & ...
        data_allab.AbsentLoc == unique(modified_data.AbsentLoc),:);
    
    % find main difference of modified patch from the original --> which
    % layer is the difference located?
    
     modified_all_layers = [modified_data.R(1,:) modified_data.G(1,:) modified_data.B(1,:)];
    [M, main_difference_layer] = max(modified_all_layers);
    null_all_layers = [null_data.R null_data.G null_data.B];
    null_selected_layer = null_all_layers(:,main_difference_layer);
    null_layer_difference = abs(null_selected_layer - modified_all_layers(:,main_difference_layer));
    
    % find the null patch that has the least difference from modified out of all null patches
    null_name(img) = null_list(null_layer_difference == min(null_layer_difference)); 
    decision(img) = null_data.Decision(null_layer_difference == min(null_layer_difference));
    original_name(img) = original_list(1,:);
    
    clear modified_data
    clear null_data
    clear null_list
    clear original_list
    clear null_all_layers
    clear null_selected_layer
    clear main_difference_layer
    clear null_layer_difference
    
end
disp(null_name)
disp(decision)
disp(original_name)

%% read the patches
folderI = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\incongruent patch";
filePatternI = fullfile(folderI, '*.jpg');
theFilesI = struct2cell(dir(filePatternI));
selectnameI(1:1260,1) = string(theFilesI(1,:));

folderC = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\congruent patch";
filePatternC = fullfile(folderC, '*.jpg');
theFilesC = struct2cell(dir(filePatternC));
selectnameC(1:1260,1) = string(theFilesC(1,:));

folderS = 'D:\New photos\nishimoto patch';
filePatternS = fullfile(folderS, '*.jpg');
theFilesS = struct2cell(dir(filePatternS));
selectnameS(1:7044,1) = string(theFilesS(1,:));

% read null patch

null_patch = 
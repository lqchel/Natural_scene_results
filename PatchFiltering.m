img = initializeImage('img/balloons.png');

params = defaultSaliencyParams(img.size,'dyadic');

pyrs = makeFeaturePyramids(img, 'Color', params)
% create the saliency map
[salmap,salData] = makeSaliencyMap(img,params);

% display the conspicuity maps
% figure('Name','STB: conspicuity maps','NumberTitle','off'); 
% displayMaps({salData.CM},2);

%% patch
load('All absent decisions.mat');
load('all present patch name.mat');
load('all absent patch name.mat');
select_num = [7 12 4 15 20] - 2;


%% for each of the 5 images, select null patches with least difference on R/G/B layer

for img = 1:5
    % patch responses for current image
    modified_data = data_allab(data_allab.Image == select_num(img) & data_allab.PatchType ~=0 & data_allab.Congruency == 0,:);
    null_data = data_allab(data_allab.Image == select_num(img) & data_allab.PatchType ==0 & data_allab.Congruency == 0 & ...
        data_allab.AbsentLoc == unique(modified_data.AbsentLoc) & data_allab.Decision == 3.5,:);
    
    % patch file names
    null_list = absent_list(data_allab.Image == select_num(img) & data_allab.PatchType ==0 & data_allab.Congruency == 0 & ...
        data_allab.AbsentLoc == unique(modified_data.AbsentLoc)& data_allab.Decision ==3.5,:); 
    original_list = present_list(data_allab.Image == select_num(img) & data_allab.PatchType ==0 & data_allab.Congruency == 0 & ...
        data_allab.AbsentLoc == unique(modified_data.AbsentLoc)& data_allab.Decision == 3.5,:);

    % find the null patch that has the least difference from modified out of all null patches
    null_name(img) = null_list(:,1); 
    decision(img) = null_data.Decision(1,:);
    original_name(img) = original_list(:,1);
    
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
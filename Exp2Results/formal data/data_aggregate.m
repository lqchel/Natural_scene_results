clear all

batch = 1:4;
group = 1:4;

for b = 1:length(batch)
    for g = 1:length(group)
        folder = ['C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Exp2Results\formal data\b'...
            num2str(batch(b)) '_g' num2str(group(g))];
        data_file_path = fullfile(folder,'*.xlsx');
        theFiles = struct2cell(dir(data_file_path));
        data_file_name(:,1:size(theFiles,2)) = string(theFiles(1,:));
        within_group_sub_code = [10 11 12 13 14 15 1 2 3 4 5 6 7 8 9];
        subject_code = (batch(b)-1).*60 + (group(g)-1).*15 + within_group_sub_code;

            for i = 1:length(data_file_name)
                current_subject_file = importdata(data_file_name(:,i));
                current_subject_data = current_subject_file.data(current_subject_file.data(:,2)>4, 1:12); % remove 4 practice trials
                sub_catch_accuracy = catch_accuracy(current_subject_data);
%                 if sub_catch_accuracy.Accuracy <= 0.38 % if current participant's catch accuracy is lower than 0.38, then excluded from further analysis
%                     sprintf('Batch = %d, Group = %d, Excluded subject = %d', batch(b), group(g), within_group_sub_code(i))
%                     continue
%                 end
%                 current_file_id = current_subject_data(:,1);
                current_subject_data(:,1) = subject_code(i);
                if i == 1 
                    aggregate_data = current_subject_data;
%                     file_id = current_file_id;
%                 elseif i > 1 && ~exist('aggregate_data')
%                     aggregate_data = current_subject_data;
%                     file_id = current_file_id;
                else
                    aggregate_data = [aggregate_data; current_subject_data];
%                     file_id = [file_id; current_file_id];
                end
                clear current_subject_file
                clear current_subject_data
            end

        data_sorted = [aggregate_data(aggregate_data(:,1)>= subject_code(7) & aggregate_data(:,1)<= subject_code(15),:); ...
            aggregate_data(aggregate_data(:,1)>subject_code(15),:)];
%         file_id_sorted = unique([file_id(aggregate_data(:,1) >= subject_code(7) & aggregate_data(:,1) <= subject_code(15));...
%             file_id(aggregate_data(:,1)>9,:)],'stable');
        
        if g == 1
            batch_data = data_sorted;
        else 
            batch_data = [batch_data; data_sorted];
        end
        
     end
    
    if b == 1
        all_data_matrix = batch_data;
    else
        all_data_matrix = [all_data_matrix; batch_data];
    end
end

%% save file for current group
save_path = sprintf('Exp2_data.mat');
save(save_path,'all_data_matrix','-mat');


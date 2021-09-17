%%% this code checks if each participant in Exp 2 gets unique null patches

function [out] = null_unique(data)
   data = data(data(:,11)~= 0,:);
   subject_id = unique(data(:,1));
   for sub = 1:length(subject_id)
       current_subject_null = data(data(:,1)== subject_id(sub) & data(:,12)~=0,12);
       current_sub_img = data(data(:,1)== subject_id(sub) & data(:,12)~=0,11);
       [img_sequence,order] = sort(current_sub_img);
       if sub == 1
           aggregate_null = current_subject_null(order,:);
       else
           aggregate_null = [aggregate_null current_subject_null];
       end
       clear current_subject_null
       clear current_sub_img
   end
   if size(unique(aggregate_null,'rows'))== size(aggregate_null)
       out = 1; % all image-null patch combinations are unique
   else 
       out = 0; % there are overlapping patches for the same images across participants
   end
end
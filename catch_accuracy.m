%%%%% this function calculates the catch trial accuracy of online
%%%%% experiment participants

function [out] = catch_accuracy(data)
subject_id = unique(data(:,1));
% real catch accuracy
catch_accuracy = table(zeros(length(subject_id),1),zeros(length(subject_id),1),'VariableNames',{'Subject','Accuracy'});
for sub = 1:length(subject_id)
    current_subject = data(data(:,1) == subject_id(sub) & data(:,11) == 0 & data(:,2) >4,:);
    catch_accuracy.Accuracy(sub,:) = sum(current_subject(:,5) == current_subject(:,8)& current_subject(:,6) == current_subject(:,9))/...
        size(current_subject,1);
    catch_accuracy.Subject(sub,:) = subject_id(sub);
end
out = catch_accuracy;
end
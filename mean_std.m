

function [m,sd] = mean_std(data,num_sub)
for sub = 1:num_sub
    R_indv = data(data(:,1)==sub,:);
    m(sub) = mean(R_indv(:,9));
    sd(sub) = std(R_indv(:,9));
end

end
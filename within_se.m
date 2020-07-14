% this function calculates within-subject standard errors 
% for within-subject, multiple-condition data
% A for the data matrix, sub for number of rows --> subject numbers,
% and type for number of conditions --> number of columns

function[se] = within_se(A,sub,type)
indv_mean = nanmean(A,2);
grandmean = nanmean(nanmean(A));
    for s = 1:sub
        B(s,:) = A(s,:)- indv_mean(s,:) + grandmean;
    end
    
    for d = 1:type
        se(:,d)= nanstd(B(:,d))/sqrt(sum(~isnan(B(:,d))));
    end

end
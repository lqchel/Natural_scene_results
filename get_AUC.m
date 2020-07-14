% This function calculates the Type 1 accuracy of any two-alternative force
% choice judgment task with 4 confidence levels(1 = not confidence at all, 4 
% = highly confident). First entry = data matrix for present stimuli,
% second entry for data matrix for absent stimuli, and
% final column is the total number of subjects.


function [AUC] = get_AUC(PresentData, AbsentData, total_num_sub)
matrix1 = zeros(total_num_sub,1);
for sub = 1:total_num_sub
  
    Confidence_N = AbsentData(AbsentData(:,1)==sub,9);
    Confidence_AP = PresentData(PresentData(:,1)==sub,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    for i = 1:9
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i) = Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
        Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
        Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix1(sub,:)= AUC;

clear Confidence_N
clear Confidence_AP

end
end

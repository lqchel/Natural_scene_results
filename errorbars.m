clear all;
Results = importdata('Pooled Results.mat');

%% patch index
Results(:,9) = Results(:,9)-0.5;
Results(:,9) = Results(:,8).* Results(:,9);
% signal present and absent trial classification

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2
Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;

%% location type
central = Results(:,7)==5;
periphery_1= 2.*(Results(:,7)==2|Results(:,7)==4|Results(:,7)==6|Results(:,7)==8);
periphery_2 = 3.* (Results(:,7)==1|Results(:,7)==3|Results(:,7)==7|Results(:,7)==9);
location = central + periphery_1 + periphery_2;

Results = [Results location];

%% AP v.s. N
grandmatrix = zeros(4,3,2);
for cong = 1:2
   if cong == 1
       Find_P = Find_CAP;
       Find_CN = Find_N & Results(:,4) ==0;
   else
       Find_P = Find_IAP;
       Find_CN = Find_N & Results(:,4) ==1;
   end
   
    Results_P = Results(Find_P,:);
    matrix1 = zeros(15,3);
    grandmean1 = mean(Results_P(:,9));
    for sub = 1:15
        R_indv = Results_P(Results_P(:,1)==sub,:);
        indv_m = mean(R_indv(:,9));
        for condition = 1:3
            condition_m = mean(R_indv(R_indv(:,end)== condition,9));
            matrix1(sub,condition)= condition_m- indv_m + grandmean1;
        end
    end

    for condition = 1:3
        conditionmean1(condition) = mean(matrix1(:,condition));
        se1(condition) = std(matrix1(:,condition))/sqrt(15);
    end
    
    clear R_indv
    clear condition_m
    clear indv_m

    Results_N = Results(Find_CN, :);
    matrix2 = zeros(15,3);
    grandmean2 = mean(Results_N(:,9));

    for sub = 1:15
        R_indv = Results_N(Results_N(:,1)==sub,:);
        indv_m = mean(R_indv(:,9));
        for condition = 1:3
            condition_m = mean(R_indv(R_indv(:,end)== condition,9));
            matrix2(sub,condition)= condition_m- indv_m + grandmean2;
        end
    end

    for condition = 1:3
        conditionmean2(condition) = mean(matrix2(:,condition));
        se2(condition) = std(matrix2(:,condition))/sqrt(15);
    end
    grandmatrix(1:4,:,cong)= [conditionmean1; se1; conditionmean2; se2];
    
    clear Results_P
    clear Results_N
    clear matrix1
    clear matrix2
end

for cong = 1:2
    if cong == 1
        congruency = 'congruent';
    else 
        congruency = 'incongruent';
    end
subplot(1,2,cong),errorbar([0 7.74 10.95],grandmatrix(1,:,cong),grandmatrix(2,:,cong),'bo-',...
    'MarkerSize',2, 'LineWidth',1), ylim([-4 4]), xticks([0 7.74 10.95]),xlabel('Eccentricity by degree of visual angle'),...
    ylabel('Confidence x response');
hold on
errorbar([0 7.74 10.95],grandmatrix(3,:,cong),grandmatrix(4,:,cong),'ro-',...
    'MarkerSize',2, 'LineWidth',1);
hold off
title(['AP v.s N for ',congruency,' images']);
legend('AP','N');
end

%% congruent v.s. incongruent

grandmatrix = zeros(4,3,2);

for cong = 1:2
   if cong == 1
       Find_P = Find_Congruent_CP;
       Find_CN = Find_Congruent_IP;
   else
       Find_P = Find_Incongruent_IP;
       Find_CN = Find_Incongruent_CP;
   end
   
    Results_P = Results(Find_P,:);
    matrix1 = zeros(15,3);
    grandmean1 = mean(Results_P(:,9));
    for sub = 1:15
        R_indv = Results_P(Results_P(:,1)==sub,:);
        indv_m = mean(R_indv(:,9));
        for condition = 1:3
            condition_m = mean(R_indv(R_indv(:,end)== condition,9));
            matrix1(sub,condition)= condition_m- indv_m + grandmean1;
        end
    end

    for condition = 1:3
        conditionmean1(condition) = mean(matrix1(:,condition));
        se1(condition) = std(matrix1(:,condition))/sqrt(15);
    end
    
    clear R_indv
    clear condition_m
    clear indv_m

    Results_N = Results(Find_CN, :);
    matrix2 = zeros(15,3);
    grandmean2 = mean(Results_N(:,9));

    for sub = 1:15
        R_indv = Results_N(Results_N(:,1)==sub,:);
        indv_m = mean(R_indv(:,9));
        for condition = 1:3
            condition_m = mean(R_indv(R_indv(:,end)== condition,9));
            matrix2(sub,condition)= condition_m- indv_m + grandmean2;
        end
    end

    for condition = 1:3
        conditionmean2(condition) = mean(matrix2(:,condition));
        se2(condition) = std(matrix2(:,condition))/sqrt(15);
    end
    grandmatrix(1:4,:,cong)= [conditionmean1; se1; conditionmean2; se2];
    
    clear Results_P
    clear Results_N
    clear matrix1
    clear matrix2
end

for cong = 1:2
    if cong == 1
        congruency = 'congruent';
    else 
        congruency = 'incongruent';
    end
subplot(1,2,cong),errorbar([0 7.74 10.95],grandmatrix(1,:,cong),grandmatrix(2,:,cong),'bo-',...
    'MarkerSize',2, 'LineWidth',1), ylim([-4 4]),xticks([0 7.74 10.95]),xlabel('Eccentricity by degree of visual angle'),...
    ylabel('Confidence x response');
hold on
errorbar([0 7.74 10.95],grandmatrix(3,:,cong),grandmatrix(4,:,cong),'ro-',...
    'MarkerSize',2, 'LineWidth',1);
hold off
title(['Critical objects for ',congruency,' images']);
legend('present','absent');
end




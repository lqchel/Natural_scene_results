%%% this code calculates the decision biases across eccentricity!


Results = importdata('Pooled Results 2.mat'); %uncomment for lab testing
%results
% Results = data;
num_sub = 15;
Results(:,9) = Results(:,8).*Results(:,9);

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

%Congruent trial with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1;
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2

Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];

%% present vs null patches

colours = cbrewer('qual', 'Set1', 8); 
grandmatrix1 = zeros(2,3);

% hit and FA rate
matrix1 = zeros(num_sub,3);
location = [0 6.5 9.2];

for loc = 1:3
for sub = 1:num_sub
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N)& Results(:,end)==location(loc),:);

for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        hit_1 = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        FA_1 = sum(Results_P(:,8)==1)/size(Results_P,1);
    end
 
    

end
    matrix1(sub,loc) = -0.5.*(norminv(hit_1)+norminv(FA_1));
    clear Results_P
    clear R_indv
end
    
end

    
grandmatrix1(1,:) = mean(matrix1,1);
grandmatrix1(2,:)= within_se(matrix1,num_sub,3);

%% original vs different, congruent condition

matrix2 = zeros(num_sub,3);
grandmatrix2 = zeros(2,3);
    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            hit_2 = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            FA_2 = sum(Results_P(:,8)==1)/size(Results_P,1);
        end


    clear Results_P
    end

       if FA_2 == 1
           FA_2 = 0.99;
       elseif hit_2 ==1
       
           hit_2 = 0.99;
           
       elseif hit_2 == 0
           hit_2 = 0.01;
       elseif FA_2 == 0
           FA_2 = 0.01;
       end
       
       if sub == 5
           disp(FA_2); disp(hit_2);
       end
           matrix2(sub,loc) = -0.5.*(norminv(hit_2)+norminv(FA_2));

    
    clear R_indv
    clear hit_2
    clear FA_2
    end

    end
grandmatrix2(1,:) = nanmean(matrix2,1);
grandmatrix2(2,:)= within_se(matrix2,num_sub,3);

%% original vs. different, incongruent condition

matrix3 = zeros(num_sub,3);
grandmatrix3 = zeros(2,3);
    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP |Find_Incongruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            hit_3 = sum(Results_P(:,8)==1)/size(Results_P,1);
            
            if hit_3 == 1
                
            else
            end
        else
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            FA_3 = sum(Results_P(:,8)==1)/size(Results_P,1);
        end
  
       
    clear Results_P
    

    end

    if sub == 6 && loc == 1
        disp(hit_3); disp(FA_3);
        disp(norminv(hit_3,0,1)); 
    end
    

      matrix3(sub,loc) = -0.5.*(norminv(hit_3,0,1)+norminv(FA_3,0,1));
     
   
    clear R_indv
    clear hit_3
    clear FA_3
    end

    end
grandmatrix3(1,:) = nanmean(matrix3,1);
grandmatrix3(2,:)= within_se(matrix3,num_sub,3);

%% plot
errorbar([0 6.5 9.2],grandmatrix1(1,:),grandmatrix1(2,:),'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',0.8,'Capsize',10);
hold on
errorbar([0 6.5 9.2],grandmatrix2(1,:),grandmatrix2(2,:),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);

errorbar([0 6.5 9.2],grandmatrix3(1,:),grandmatrix3(2,:),'d--','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);

hold off



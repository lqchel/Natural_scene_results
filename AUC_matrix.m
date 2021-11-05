%%% this function produces Type 1 and Type 2 AUC results in a structure for Qianchen's experiment. data =
%%% data matrix, num_sub = number of subjects in the data matrix.
%%% for this function to work, please make sure to install PsychToolbox (http://psychtoolbox.org/)


function out = AUC_matrix(data,exp)
addpath(genpath('C:\Users\liang\Documents\Experiment Codes\Natural_scene_results\cbrewer')); % add colour palatte package path
if exp == 1
    Results = data;
    location1 = Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8;
    location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 2;
    eccentricity_level = zeros(length(Results),1)+ location1 + location2;
    real_location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
    real_location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
    actual_eccentrcity = zeros(length(Results),1) + real_location1 + real_location2;
    Results = [Results eccentricity_level actual_eccentrcity];
    
    Find_N = Results(:,5) ==1;
    % present patch trials -- signal present for hypo 1
    Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
    Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

    %Congruent trial with Congruent object, and incongruent trial with
    %incongruent object -- signal present for hypo 2
    Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3; %% Results(:,6) == 1 for exp 2
    Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;

    %Incongruent trial with congruent object, congruent trial with incongruent
    %object -- signal absent for hypo 2
    Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
    Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
    
else
    Results = data(data(:,11)~=0,:); %% exclude catch trials    

    Find_N = Results(:,5) ==1;
    % present patch trials -- signal present for hypo 1
    Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
    Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 

    %Congruent trial with Congruent object, and incongruent trial with
    %incongruent object -- signal present for hypo 2
    Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 1; %% Results(:,6) == 1 for exp 2
    Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 1;

    %Incongruent trial with congruent object, congruent trial with incongruent
    %object -- signal absent for hypo 2
    Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 1;
    Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 1;

end

location = [0 1 2];

subject_id = unique(Results(:,1));
num_sub = length(subject_id);

Results(:,9) = Results(:,8).*Results(:,9);

Results_NC = Results(Find_N,:); % find absent patches
Results_APC = Results(Find_IAP|Find_CAP,:); % find present patches

%% across eccentricities, hypothesis 1
b = 1;
for sub = 1:num_sub 

    Results_NC = Results(Find_N,:);
    Confidence_N = Results_NC(Results_NC(:,1)==subject_id(sub),9);
    Confidence_AP = Results_APC(Results_APC(:,1)==subject_id(sub),9);
    
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
    matrix1(b,:)= [sub -1 -1 AUC]; % subject number, eccentricity level, actual eccentricity level, AUC 
    b = b+1;

end

%% hypothesis 1, on each eccentricity level
b = 1;
for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_N,:);
    indvP = Results(Results(:,1)==subject_id(sub) & (Find_CAP|Find_IAP),:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,13)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,13) == location(a),:);
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i); % AUC calculation
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
    matrix2(b,:)= [sub location(a) unique(indvN(indvN(:,13)== location(a),14)) AUC];
    b = b+1;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end
matrix2 = matrix2(~isnan(matrix2(:,end)),:);
sub_AUC.hypo1_Type1 = [matrix1; matrix2];

%% hypothesis 2 AUC
b = 1;
for condition = 1:2
    
    for sub = 1:num_sub
    if condition ==1
        Confidence_P = Results(Find_Congruent_CP & Results(:,1)==subject_id(sub),9);
        Confidence_A = Results(Find_Congruent_IP & Results(:,1)==subject_id(sub),9);
    else
        Confidence_P = Results(Find_Incongruent_IP & Results(:,1)==subject_id(sub),9);
        Confidence_A = Results(Find_Incongruent_CP & Results(:,1)==subject_id(sub),9);
    end

for i = -4:4
    Confidence_APCounts(i+5) = nansum(Confidence_P == -i);
    Confidence_NCounts(i+5) = nansum(Confidence_A == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_A);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_P);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_A);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_P);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

matrix3(b,:) = [sub condition -1 -1 AUC];
b = b+1;

clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end

% across eccentricities

%congruent
b = 1;
for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_Congruent_IP,:);
    indvP = Results(Results(:,1)==subject_id(sub) & Find_Congruent_CP,:); % trial classification
for a = 1:3
    
        indvN_loc = indvN(indvN(:,13)== location(a),:); % select trials on that location
        indvP_loc = indvP(indvP(:,13) == location(a),:);
        if size(indvN_loc,1) == 0 || size(indvP_loc,1) == 0
            continue
        end
        Confidence_N = indvN_loc(:,9);
        Confidence_AP = indvP_loc(:,9);

        for i = -4:4
            Confidence_APCounts(i+5) = nansum(Confidence_AP == -i); % AUC calculation
            Confidence_NCounts(i+5) = nansum(Confidence_N == -i);
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
        matrix4(b,:)= [sub 1 location(a) indvN_loc(1,14) AUC];
        b = b+1;

clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC
    
end
end

matrix4 = matrix4(~isnan(matrix4(:,end)),:);

%incongruent condition

clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff

b = 1;
for sub = 1:num_sub
    indvN = Results(Results(:,1)==subject_id(sub) & Find_Incongruent_CP,:);
    indvP = Results(Results(:,1)==subject_id(sub) & Find_Incongruent_IP,:); % trial classification
for a = 1:3
    indvN_loc = indvN(indvN(:,13)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,13) == location(a),:);
    if size(indvN_loc,1) == 0 || size(indvP_loc,1) == 0
        continue
    end
    Confidence_N = indvN_loc(:,9);
    Confidence_AP = indvP_loc(:,9);
    
    for i = -4:4
        Confidence_APCounts(i+5) = nansum(Confidence_AP == -i); % AUC calculation
        Confidence_NCounts(i+5) = nansum(Confidence_N == -i);
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
    matrix5(b,:)= [sub 2 location(a) indvN_loc(1,14) AUC];
    b = b+1;
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end
matrix5 = matrix5(~isnan(matrix5(:,end)),:);
sub_AUC.hypo2_Type1 = [matrix3; matrix4; matrix5];

%% hypothesis 1 type 2 AUC

Results(:,9) = abs(Results(:,9));
Results_NC = Results(Find_N,:); % all trials with absent test probes
Results_APC = Results(Find_IAP|Find_CAP,:); % all trials with present test probes
Results_Correct = [Results_NC(Results_NC(:,8)==-1,:); Results_APC(Results_APC(:,8)==1,:)];
Results_Incorrect = [Results_NC(Results_NC(:,8)==1,:); Results_APC(Results_APC(:,8)==-1,:)];
b = 1;

for sub = 1:num_sub 
    Confidence_Incorrect = Results_Incorrect(Results_Incorrect(:,1)==subject_id(sub),9);
    Confidence_Correct = Results_Correct(Results_Correct(:,1)==subject_id(sub),9);
    if size(Confidence_Incorrect,1) == 0 || size(Confidence_Correct,1) == 0
        continue
    end
%     Confidence_Incorrect = Results_Incorrect(:,9);
%     Confidence_Correct = Results_Correct(:,9);
    
    for i = -4:-1
        Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
        Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
    end
    for i = 1:4
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i) = Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit];
    Cumulative_FA = [0 Cumulative_FA];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

    matrix6(b,:)= [sub -1 -1 AUC]; % subject number, eccentricity level, actual ecc, AUC
    b = b+1;
%     hit_matrix(sub,:,a) = Cumulative_Hit;
%     FA_matrix(sub,:,a) = Cumulative_FA;
    
    clear Cumulative_Hit
    clear Cumulative_FA
    clear AUC
 
end

clear Confidence_Correct
clear Confidence_Incorrect

%% on each eccentricity levels, hypo1, type 2
b = 1;
for a = 1:3
    for sub = 1:num_sub 
    Confidence_Incorrect = Results_Incorrect(Results_Incorrect(:,1)==subject_id(sub) & Results_Incorrect(:,13)==location(a),9);
    Confidence_Correct = Results_Correct(Results_Correct(:,1)==subject_id(sub)& Results_Correct(:,13)==location(a),9);
    if size(Confidence_Incorrect,1) == 0 || size(Confidence_Correct,1) == 0
        continue
    end
%     Confidence_Incorrect = Results_Incorrect(:,9);
%     Confidence_Correct = Results_Correct(:,9);
    for i = -4:-1
        Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
        Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
    end
    for i = 1:4
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i) = Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit];
    Cumulative_FA = [0 Cumulative_FA];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix7(b,:)= [sub location(a) unique(Results_Correct(Results_Correct(:,1)==subject_id(sub)& Results_Correct(:,13)==location(a),14)) AUC]; % subject number, eccentricity level, actual ecc, AUC
    b = b+1;
%     hit_matrix(sub,:,a) = Cumulative_Hit;
%     FA_matrix(sub,:,a) = Cumulative_FA;
    
    clear Cumulative_Hit
    clear Cumulative_FA
    clear AUC
 
    end
end

clear Confidence_Correct
clear Confidence_Incorrect
matrix7 = matrix7(~isnan(matrix7(:,end)),:);
sub_AUC.hypo1_Type2 = [matrix6; matrix7];

%% hypothesis 2 type 2 AUC 
b = 1;
for condition = 1:2
    
    for sub = 1:num_sub
        if condition ==1
            Results_N = Results(Find_Congruent_IP,:);
            Results_A = Results(Find_Congruent_CP,:);
            Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
            Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
        else
            Results_N = Results(Find_Incongruent_CP,:);
            Results_A = Results(Find_Incongruent_IP,:);
            Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
            Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
        end
        if size(Confidence_Incorrect,1) == 0 || size(Confidence_Correct,1) == 0
            continue
        end
       for i = -4:-1
            Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
            Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
        end
        for i = 1:4
            if i == 1
            Cumulative_NCounts(i) = Confidence_NCounts(i);
            Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
            Confidence_YCounts(i) = Confidence_YCounts(i);
            Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
            else
            Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
            Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
            Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
            Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
            end
        end

        Cumulative_Hit = [0 Cumulative_Hit];
        Cumulative_FA = [0 Cumulative_FA];
        AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
        matrix8(b,:)= [sub condition -1 -1 AUC]; % create matrix for individual AUCs
        b = b+1;
        clear Confidence_Correct
        clear Confidence_Incorrect
        clear Cumulative_NCounts;
        clear Cumulative_APCounts;
        clear Cumulative_Hit;
        clear Cumulative_FA;
        clear AUC;
    end

end
matrix8 = matrix8(~isnan(matrix8(:,end)),:);

%% across eccentricities

%congruent

b = 1;
for condition = 1:3
    
    for sub = 1:num_sub
        Results_N = Results(Find_Congruent_IP & Results(:,13)== location(condition),:);
        Results_A = Results(Find_Congruent_CP & Results(:,13) == location(condition),:);
        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
     if size(Confidence_Incorrect,1) == 0 || size(Confidence_Correct,1) == 0
        continue
     end
   for i = -4:-1
        Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
        Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
    end
    for i = 1:4
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i) = Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit];
    Cumulative_FA = [0 Cumulative_FA];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix9(b,:)= [sub 1 location(condition) unique(Results(Results(:,13)== location(condition) & Results(:,1) == subject_id(sub),14)) AUC]; % create matrix for individual AUCs
    b = b+1;
    
clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end
matrix9 = matrix9(~isnan(matrix9(:,end)),:);

%% incongruent condition
clear Results_N
clear Results_A

b = 1;
for condition = 1:3
    
    for sub = 1:num_sub
        Results_N = Results(Find_Incongruent_CP & Results(:,13)== location(condition),:);
        Results_A = Results(Find_Incongruent_IP & Results(:,13) == location(condition),:);

        Confidence_Correct = [Results_N(Results_N(:,1)== subject_id(sub) & Results_N(:,8)==-1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==1,9)];
        Confidence_Incorrect = [Results_N(Results_N(:,1)== subject_id(sub) &Results_N(:,8)==1,9); Results_A(Results_A(:,1)==subject_id(sub) & Results_A(:,8)==-1,9)];
 
        if size(Confidence_Incorrect,1) == 0 || size(Confidence_Correct,1) == 0
            continue
        end
    for i = -4:-1
        Confidence_YCounts(i+5) = sum(Confidence_Correct == -i);
        Confidence_NCounts(i+5) = sum(Confidence_Incorrect == -i);
    end
    for i = 1:4
        if i == 1
        Cumulative_NCounts(i) = Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i) = Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        else
        Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
        Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_Incorrect);
        Confidence_YCounts(i)= Confidence_YCounts(i-1)+ Confidence_YCounts(i);
        Cumulative_Hit(i) = Confidence_YCounts(i)/length(Confidence_Correct);
        end
    end

    Cumulative_Hit = [0 Cumulative_Hit];
    Cumulative_FA = [0 Cumulative_FA];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    matrix0(b,:)= [sub 2 location(condition) unique(Results(Results(:,13)== location(condition)&Results(:,1) == subject_id(sub),14)) AUC]; % create matrix for individual AUCs
    b = b+1;
clear Confidence_Correct
clear Confidence_Incorrect
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end
matrix0 = matrix0(~isnan(matrix0(:,end)),:);
sub_AUC.hypo2_Type2 = [matrix8;matrix9; matrix0];
out = sub_AUC;

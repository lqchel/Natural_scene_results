%% import data and index different types of patches
Results = importdata('Pooled Results.mat');

Results(:,9) = Results(:,8).*Results(:,9);
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

% type 1 ROC curve analysis
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);

%% hypothesis 1 analysis
% hit and CR 
matrix1 = zeros(15,2);

for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
 
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
    end
    matrix1(sub,condition) = accuracy;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

[h,p,ci,stats] = ttest(matrix1(:,1),matrix1(:,2)) % compare accuracy for present and absent patches

% paste output here














% AUC analysis
Results_NC = Results(Find_N,:); % find absent patchese
Results_APC = Results(Find_IAP|Find_CAP,:); % find present patches
matrix1 = zeros(15,1);
for sub = 1:15 
    Results_NC = Results(Find_N,:);
    Confidence_N = Results_NC(Results_NC(:,1)==sub,9);
    Confidence_AP = Results_APC(Results_APC(:,1)==sub,9);
    
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



end
[h,p,ci,stats]= ttest(matrix1(:,1),0.5) % test AUC against chance

% paste output here














%% hypothesis 2

% AUC calculation
grandmatrix = zeros(15,2);
for condition = 1:2
    
    for sub = 1:15
    if condition ==1
        Confidence_P = Results(Find_Congruent_CP & Results(:,1)==sub,9);
        Confidence_A = Results(Find_Congruent_IP & Results(:,1)==sub,9);
    else
        Confidence_P = Results(Find_Incongruent_IP & Results(:,1)==sub,9);
        Confidence_A = Results(Find_Incongruent_CP & Results(:,1)==sub,9);
    end

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_P == -i);
    Confidence_NCounts(i+5) = sum(Confidence_A == -i);
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

grandmatrix(sub,condition)= AUC;

clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end

% t test
[h,p,ci,stats] = ttest(grandmatrix(:,1),0.5,'Alpha',0.025) % congruent AUC against chance, alpha corrected to 0.025
[h,p,ci,stats] = ttest(grandmatrix(:,2),0.5,'Alpha',0.025) % incongruent AUC against chance, alpha corrected to 0.025
[h,p,ci,stats] = ttest(grandmatrix(:,1),matrix1,'Alpha',0.025) % congruent AUC compared with AUC for hypothesis 1, alpha = 0.025
[h,p,ci,stats] = ttest(grandmatrix(:,2),matrix1,'Alpha',0.025)
[h,p,ci,stats] = ttest(grandmatrix(:,1),grandmatrix(:,2)) % congruent v.s. incongruent AUC, alpha = 0.05

% paste output here


















%% hypothesis 3 analysis

% create column showing eccentricity of each patch
location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];
% hit and CR analysis
matrix1 = zeros(90,4);
location = [0 6.48 9.16];

for loc = 1:3
for sub = 1:15
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N)& Results(:,end)==location(loc),:);

for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
    end
    rowindex = size(matrix1,1)-sum(matrix1(:,1)==0)+1;
    
    matrix1(rowindex,1) = sub;
    matrix1(rowindex,2) = condition - 1;
    matrix1(rowindex,3) = location(loc);
    matrix1(rowindex,4) = accuracy;
    
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

end


% fit lme, one for hit, one for CR
data1 = table(matrix1(matrix1(:,2)==0,1),matrix1(matrix1(:,2)==0,3),matrix1(matrix1(:,2)==0,4),...
    'VariableName',{'Subjects','Location','Accuracy'});
lm1 = fitlme(data1,'Accuracy ~ Location  + (1|Subjects) + (Location|Subjects) ');
lm2 = fitlme(data1,'Accuracy ~ 1+ (1|Subjects) + (Location|Subjects)');
compare(lm2,lm1)
% assumption checks
subplot(1,2,1), plotResiduals(lm1, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm1,'fitted'), title('Residual against predicted values'),set(gca,'FontSize',14);

% paste output here








data2 = table(matrix1(matrix1(:,2)==1,1),matrix1(matrix1(:,2)==1,3),matrix1(matrix1(:,2)==1,4),...
    'VariableName',{'Subjects','Location','Accuracy'});
lm3 = fitlme(data2,'Accuracy ~ Location  + (1|Subjects) + (Location|Subjects) ');
lm4 = fitlme(data2,'Accuracy ~ 1+ (1|Subjects) + (Location|Subjects)');
compare(lm4,lm3)

%assumption checks
subplot(1,2,1), plotResiduals(lm3, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm3,'fitted'), title('Residual against predicted values'),set(gca,'FontSize',14);


%paste output here


% AUC calculation, on each eccentricity levels 
matrix3 = zeros(15, 3);
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_N,:);
    indvP = Results(Results(:,1)==sub & (Find_CAP|Find_IAP),:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
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
    matrix3(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end

% sort AUC, location and individual ID into a table, then fit lme model
AUC = reshape(matrix3,[45,1]);
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[3,1]);
location = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];


data = table(AUC, subjects, location,'VariableName',{'AUC','subject','location'});
lm1 = fitlme(data, 'AUC ~ location + (location|subject) + (1|subject)');
lm2 = fitlme(data, 'AUC ~ 1 +(location|subject) + (1|subject) ');

compare(lm2,lm1)

%assumption checks for fitting lm1
subplot(1,2,1), plotResiduals(lm1, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm1,'fitted'), title('Residual against predicted values'),set(gca,'FontSize',14);



% paste output here





%% hypothesis 4

%congruent
clear indvN
clear indvP
clear indvN_loc
clear indvP_loc

matrix4 = zeros(15, 3);
location = [0 6.48 9.16];
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Congruent_IP,:);
    indvP = Results(Results(:,1)==sub & Find_Congruent_CP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
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
    matrix4(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC
    
end
end


%sort AUC, location and individual ID into a table, then fit lme model
AUC = reshape(matrix4,[45,1]);
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[3,1]);
location = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];
data = table(AUC, subjects, location,'VariableName',{'AUC','subject','location'});

lm1 = fitlme(data, 'AUC ~ location + (location|subject) + (1|subject)');
lm2 = fitlme(data, 'AUC ~ 1 +(location|subject) + (1|subject) ');
% assumption check
subplot(1,2,1), plotResiduals(lm1, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm1,'fitted'), title('Residual against predicted values'),set(gca,'FontSize',14);


compare(lm2,lm1)
% output


%incongruent condition

clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff

matrix5 = zeros(15, 3);

for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Incongruent_CP,:);
    indvP = Results(Results(:,1)==sub & Find_Incongruent_IP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
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
    matrix5(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    
end
end


AUC = [reshape(matrix4,[45,1]); reshape(matrix5,[45,1])];
clear subjects
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[6,1]);
eccentricity = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];
eccentricity = [eccentricity; eccentricity];
condition = [zeros(45,1);ones(45,1)];

data = table(AUC, subjects, eccentricity,condition,'VariableName',{'AUC','subject','eccentricity','condition'});
lm4 = fitlme(data, 'AUC ~ eccentricity*condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml1 = fitlme(data, 'AUC ~ eccentricity:condition + condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml2 = fitlme(data, 'AUC ~ eccentricity + condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml3 = fitlme(data, 'AUC ~ eccentricity:condition + eccentricity + (eccentricity|subject) + (condition|subject)+(1|subject)');

c_ecc = compare(lml1,lm4)
c_interaction = compare(lml2,lm4)
c_condition = compare(lml3,lm4)

%post-hoc tests
data2 = table(AUC(1:45,:),subjects(1:45,:),eccentricity(1:45,:), condition(1:45,:),'VariableName',{'AUC','subject','eccentricity','condition'});
lm5 = fitlme(data2, 'AUC ~ eccentricity + (eccentricity|subject)+(1|subject)')
lml4 = fitlme(data2, 'AUC ~ 1 + (eccentricity|subject)+(1|subject)')

data3 = table(AUC(46:90,:),subjects(46:90,:),eccentricity(46:90,:), condition(46:90,:),'VariableName',{'AUC','subject','eccentricity','condition'});
lm6 = fitlme(data3, 'AUC ~ eccentricity + (eccentricity|subject) + (1|subject)')
lml5 = fitlme(data3, 'AUC ~ 1 + (eccentricity|subject) + (1|subject)')

compare(lml4,lm5)
compare(lml5,lm6)

compare(lm2,lm1)


% assumption checks



% paste output here





%incongruent condition
 
clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff
 
matrix5 = zeros(15, 3);
 
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_Incongruent_CP,:);
    indvP = Results(Results(:,1)==sub & Find_Incongruent_IP,:); % trial classification
for a = 1:3
    
    indvN_loc = indvN(indvN(:,end)== location(a),:); % select trials on that location
    indvP_loc = indvP(indvP(:,end) == location(a),:);
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
    matrix5(sub,a)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
 
    
end
end
 
 
AUC = [reshape(matrix4,[45,1]); reshape(matrix5,[45,1])];
clear subjects
subjects(1:15,1)= 1:1:15;
subjects = repmat(subjects,[6,1]);
eccentricity = [zeros(15,1); 6.48 + zeros(15,1); 9.16 + zeros(15,1)];
eccentricity = [eccentricity; eccentricity];
condition = [zeros(45,1);ones(45,1)];
 
data = table(AUC, subjects, eccentricity,condition,'VariableName',{'AUC','subject','eccentricity','condition'});
lm4 = fitlme(data, 'AUC ~ eccentricity*condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml1 = fitlme(data, 'AUC ~ eccentricity:condition + condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml2 = fitlme(data, 'AUC ~ eccentricity + condition + (eccentricity|subject) + (condition|subject)+(1|subject)');
lml3 = fitlme(data, 'AUC ~ eccentricity:condition + eccentricity + (eccentricity|subject) + (condition|subject)+(1|subject)');
 
c_ecc = compare(lml1,lm4)
c_interaction = compare(lml2,lm4)
c_condition = compare(lml3,lm4)
 
% assumption check for full model
subplot(1,2,1), plotResiduals(lm4, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm4,'fitted'), title('Residual against predicted values'),set(gca,'FontSize',14);

% output


%post-hoc tests, estimating lme model on each condition
data2 = table(AUC(1:45,:),subjects(1:45,:),eccentricity(1:45,:), condition(1:45,:),'VariableName',{'AUC','subject','eccentricity','condition'});
lm5 = fitlme(data2, 'AUC ~ eccentricity + (eccentricity|subject)+(1|subject)')
lml4 = fitlme(data2, 'AUC ~ 1 + (eccentricity|subject)+(1|subject)')

compare(lml4,lm5)

% assumption check
subplot(1,2,1), plotResiduals(lm5, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm5,'fitted'), title('Residual against predicted values'),set(gca,'FontSize',14);

% output

 
data3 = table(AUC(46:90,:),subjects(46:90,:),eccentricity(46:90,:), condition(46:90,:),'VariableName',{'AUC','subject','eccentricity','condition'});
lm6 = fitlme(data3, 'AUC ~ eccentricity + (eccentricity|subject) + (1|subject)')
lml5 = fitlme(data3, 'AUC ~ 1 + (eccentricity|subject) + (1|subject)')
 
compare(lml5,lm6)
 
 
% assumption checks
subplot(1,2,1), plotResiduals(lm6, 'histogram'),title('Residual histogram'),set(gca,'FontSize',14);
subplot(1,2,2), plotResiduals(lm6,'fitted'), title(['Residual against' newline 'predicted values']),set(gca,'FontSize',14);

% output















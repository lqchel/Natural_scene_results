 clear all
Results = importdata('Pooled & Cleaned Results.mat');
Results(:,9) = Results(:,8).*Results(:,9);
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
% 
%% hypo 1, not by congruency
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);
% type 1 ROC curve analysis

Confidence_N = Results_N(:,9);
Confidence_AP = Results_AP(:,9);

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

colours = cbrewer('qual','Set3',8);
h = area(Cumulative_FA, Cumulative_Hit);
h.FaceColor = colours(3,:);
h.EdgeColor = colours(3,:);

hold on
plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'); % plot ROC curve

plot([0 1],[0 1], 'r-'); % plot chance line
hold off
set(gca,'FontSize',12);
axis square
% 
% Confidence_N = [];
% Confidence_AP = [];
% Results_N = [];
% Results_AP = [];
% Cumulative_NCounts = [];
% Cumulative_APCounts = [];
% Cumulative_Hit = [];
% Cumulative_FA = [];
% AUC = [];
% 
% %% hypo 2, not by congruency
% Results_N = Results(Find_Congruent_IP | Find_Incongruent_CP,:);
% Results_AP = Results(Find_Congruent_CP | Find_Incongruent_IP,:);
% 
% Confidence_N = Results_N(:,9);
% Confidence_AP = Results_AP(:,9);
% for i = -4:4
%     Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%     length(Confidence_AP(Confidence_AP == 4))
%     Confidence_NCounts(i+5) = sum(Confidence_N == -i);
% end
% for i = 1:9
%     if i == 1
%     Cumulative_NCounts(i) = Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i) = Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     else
%     Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     end
% end
% 
% Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
% Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
% AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
% 
% subplot(3,2,2),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC curve for critical objects discrimination AUC = ', num2str(AUC)]);
% axis square
% 
% 
% Confidence_N = [];
% Confidence_AP = [];
% Results_N = [];
% Results_AP = [];
% Cumulative_NCounts = [];
% Cumulative_APCounts = [];
% Cumulative_Hit = [];
% Cumulative_FA = [];
% AUC = [];
% 
% 
% %% hypo 1, by congruency
% %Congruent trials
% Results_N = Results(Find_N & Results(:,4)== 0,:);
% Results_AP = Results(Find_CAP,:);
% 
% Confidence_N = Results_N(:,9);
% Confidence_AP = Results_AP(:,9);
% for i = -4:4
%     Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%   
%     Confidence_NCounts(i+5) = sum(Confidence_N == -i);
% end
% for i = 1:9
%     if i == 1
%     Cumulative_NCounts(i) = Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i) = Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     else
%     Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     end
% end
% 
% Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
% Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
% AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
% 
% subplot(3,2,3),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['A/P discrimination in congruent trials AUC = ', num2str(AUC)]);
% axis square
% 
% 
% Confidence_N = [];
% Confidence_AP = [];
% Results_N = [];
% Results_AP = [];
% Cumulative_NCounts = [];
% Cumulative_APCounts = [];
% Cumulative_Hit = [];
% Cumulative_FA = [];
% AUC = [];
% 
% %incongruent trials
% Results_N = Results(Find_N & Results(:,4)==1,:);
% Results_AP = Results(Find_IAP,:);
% 
% Confidence_N = Results_N(:,9);
% Confidence_AP = Results_AP(:,9);
% for i = -4:4
%     Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%     length(Confidence_AP(Confidence_AP == 4))
%     Confidence_NCounts(i+5) = sum(Confidence_N == -i);
% end
% for i = 1:9
%     if i == 1
%     Cumulative_NCounts(i) = Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i) = Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     else
%     Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     end
% end
% 
% Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
% Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
% AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
% 
% subplot(1,2,1),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['A/P discrimination in incongruent trials AUC = ', num2str(AUC)]);
% axis square
% 
% 
% Confidence_N = [];
% Confidence_AP = [];
% Results_N = [];
% Results_AP = [];
% Cumulative_NCounts = [];
% Cumulative_APCounts = [];
% Cumulative_Hit = [];
% Cumulative_FA = [];
% AUC = [];
% 
%% hypo 2, by congruency
%congruent trials
Results_AP = Results(Find_Congruent_CP,:);
Results_N =  Results(Find_Congruent_IP,:);

Confidence_N = Results_N(:,9);
Confidence_AP = Results_AP(:,9);
for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);


h = area(Cumulative_FA, Cumulative_Hit);
h.FaceColor = colours(3,:);
h.EdgeColor = colours(3,:);

hold on
plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'); % plot ROC curve

plot([0 1],[0 1], 'r-'); % plot chance line
hold off
set(gca,'FontSize',12);
axis square

subplot(3,2,5),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['Critical objects discrimination in congruent trials AUC = ', num2str(AUC)]);
axis square


Confidence_N = [];
Confidence_AP = [];
Results_N = [];
Results_AP = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];
% 
% %Incongruent trials
% Results_N = Results(Find_Incongruent_CP,:);
% Results_AP = Results(Find_Incongruent_IP,:);
% 
% Confidence_N = Results_N(:,9);
% Confidence_AP = Results_AP(:,9);
% for i = -4:4
%     Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%     length(Confidence_AP(Confidence_AP == 4))
%     Confidence_NCounts(i+5) = sum(Confidence_N == -i);
% end
% for i = 1:9
%     if i == 1
%     Cumulative_NCounts(i) = Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i) = Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     else
%     Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N);
%     Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP);
%     end
% end
% 
% Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
% Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
% AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
% 
% subplot(3,2,6),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['Critical objects discrimination in incongruent trials AUC = ', num2str(AUC)]);
% axis square
% 
% %% find miss & FA patches
% 
%     folderI = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\incongruent patch";
%     filePatternI = fullfile(folderI, '*.jpg');
%     theFilesI = struct2cell(dir(filePatternI));
%     selectnameI(1:1260,1) = string(theFilesI(1,:));
% 
%     folderC = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\congruent patch";
%     filePatternC = fullfile(folderC, '*.jpg');
%     theFilesC = struct2cell(dir(filePatternC));
%     selectnameC(1:1260,1) = string(theFilesC(1,:));
% 
%     folderS = 'D:\New photos\nishimoto patch';
%     filePatternS = fullfile(folderS, '*.jpg');
%     theFilesS = struct2cell(dir(filePatternS));
%     selectnameS(1:7044,1) = string(theFilesS(1,:));
%     
%     folderM = 'C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\squareimage\congruent cropped';
%     filePatternM = fullfile(folderM, '*.jpg');
%     theFilesM = struct2cell(dir(filePatternM));
%     selectnameM(1:140,1) = string(theFilesM(1,:));
%     
%     
%     folderN = 'C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\squareimage\incongruent cropped';
%     filePatternN = fullfile(folderN, '*.jpg');
%     theFilesN = struct2cell(dir(filePatternN));
%     selectnameN(1:140,1) = string(theFilesN(1,:));
%     
% 
% % for congruent images
% CAP_error = Find_CAP & Results(:,9) == -4;
% CN_error = Find_N & Results(:,4) == 0 & Results(:,9) == 4;
% 
% C_error = Results(CAP_error | CN_error,:);
% 
% [C,ia,ic] = unique(C_error(:,[2 5 6 12]),'rows');
% 
% C_error = C_error(ia,:);
% 
% Img_num = unique(C_error(:,2));
% 
% for num = 1: length(Img_num)
%     CurrentImg = C_error(C_error(:,2)== Img_num(num),:);
%     ErrorNum = length(CurrentImg(:,1));
%     
%     if rem(ErrorNum,2) == 0
%         row_num = 2;
%         col_num = ErrorNum/2;
%     else
%         en = ErrorNum +1;
%         row_num = 2;
%         col_num = en/2;
%     end
%     
%     if Img_num(num) < 10
%        Img_path = ['SquareCongruent_00', num2str(Img_num(num)),'.jpg'];
%     else
%        Img_path = ['SquareCongruent_0', num2str(Img_num(num)),'.jpg']; 
%     end
%     
%     Img_path = char(fullfile(folderM, Img_path));
%     Image = imread(Img_path);
%     Image = imresize(Image, [300 300]);
%     panel = uint8(255)+ zeros(300, 325+col_num.*125,3,'uint8');
%     panel(1:300,1:300,:)= Image;
%    
%     for r = 1:ErrorNum
% 
%         if CurrentImg(r,5)~= 1
%            path = ['_',num2str(Img_num(num)),'_',num2str(CurrentImg(r,7)),'_'];
%            patchname = selectnameC(contains(selectnameC(:,1),path),1);
%            fullname = char(fullfile(folderC,patchname));
%         else
%            patchname = selectnameS(CurrentImg(r,12));
%            fullname = char(fullfile(folderS, patchname));
%         end
% 
%        Patch = imread(fullname);
%        Patch = [imresize(Patch, [100 100])] ;
% 
%        if r <= col_num
%            panel(1:100,325+(r-1).*125:325+(r-1).*125+99,:) = Patch;
%        else
%            panel(201:300,325+(r-col_num-1).*125:325+(r-col_num-1).*125+99,:) = Patch;
%        end
%        
%     end
% 
%     if Img_num(num) < 10
%     savepath = ['C_0', num2str(Img_num(num)),'.jpg'];
%     else
%     savepath = ['C_', num2str(Img_num(num)),'.jpg']; 
%     end
% 
% imwrite(panel, savepath); 
% 
% 
%     
% end


%for incongruent images

% IAP_error = Find_IAP & Results(:,9) == -4;
% IN_error = Find_N & Results(:,4) == 1 & Results(:,9) == 4;
% 
% I_error = Results(IAP_error | IN_error,:);
% 
% [I,ia,ic] = unique(I_error(:,[2 5 6 12]),'rows');
% 
% I_error = I_error(ia,:);
% 
% Img_num = unique(I_error(:,2));
% 
% for num = 1: length(Img_num)
%     CurrentImg = I_error(I_error(:,2)== Img_num(num),:);
%     ErrorNum = length(CurrentImg(:,1));
%     
%     if rem(ErrorNum,2) == 0
%         row_num = 2;
%         col_num = ErrorNum/2;
%     else
%         en = ErrorNum +1;
%         row_num = 2;
%         col_num = en/2;
%     end
%     
%     if Img_num(num) < 10
%        Img_path = ['SquareIncongruent_00', num2str(Img_num(num)),'.jpg'];
%     else
%        Img_path = ['SquareIncongruent_0', num2str(Img_num(num)),'.jpg']; 
%     end
%     
%     Img_path = char(fullfile(folderN, Img_path));
%     Image = imread(Img_path);
%     Image = imresize(Image, [300 300]);
%     panel = uint8(255)+ zeros(300, 325+col_num.*125,3,'uint8');
%     panel(1:300,1:300,:)= Image;
%    
%     for r = 1:ErrorNum
% 
%         if CurrentImg(r,5)~= 1
%            path = ['_',num2str(Img_num(num)),'_',num2str(CurrentImg(r,7)),'_'];
%            patchname = selectnameI(contains(selectnameI(:,1),path),1);
%            fullname = char(fullfile(folderI,patchname));
%         else
%            patchname = selectnameS(CurrentImg(r,12));
%            fullname = char(fullfile(folderS, patchname));
%         end
% 
%        Patch = imread(fullname);
%        Patch = [imresize(Patch, [100 100])] ;
% 
%        if r <= col_num
%            panel(1:100,325+(r-1).*125:325+(r-1).*125+99,:) = Patch;
%        else
%            panel(201:300,325+(r-col_num-1).*125:325+(r-col_num-1).*125+99,:) = Patch;
%        end
%        
%     end
% 
%     if Img_num(num) < 10
%     savepath = ['I_0', num2str(Img_num(num)),'.jpg'];
%     else
%     savepath = ['I_', num2str(Img_num(num)),'.jpg']; 
%     end
% 
% imwrite(panel, savepath); 
% 
% 
%     
% end

%% central v.s. periphery analysis
central = Results(:,7)==5;
periphery_1= 2.*(Results(:,7)==2|Results(:,7)==4|Results(:,7)==6|Results(:,7)==8);
periphery_2 = 3.* (Results(:,7)==1|Results(:,7)==3|Results(:,7)==7|Results(:,7)==9);
location = central + periphery_1 + periphery_2;

Results = [Results location];

locType = {'Centre','Periphery 1','Periphery 2'};
Results_N = Results(Find_N,:);
Results_AP = Results(Find_IAP|Find_CAP,:);

% hypo 1, not by congruency
for a = 1:3
    
Results_N_loc = Results_N(Results_N(:,end)==a,:);
Results_AP_loc = Results_AP(Results_AP(:,end)==a,:);

Confidence_N = Results_N_loc(:,9);
Confidence_AP = Results_AP_loc(:,9);

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N_loc);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP_loc);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N_loc);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP_loc);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(1,3,a),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['ROC for AP/N patch in',locType(a),'AUC = ', num2str(AUC)]);
axis square

Confidence_N = [];
Confidence_AP = [];
Results_N_loc = [];
Results_AP_loc = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];




end

%Incongruent trials
Results_N = Results(Find_Incongruent_CP,:);
Results_AP = Results(Find_Incongruent_IP,:);

for a = 1:3
Results_N_loc = Results_N(Results_N(:,end)==a,:);
Results_AP_loc = Results_AP(Results_AP(:,end)==a,:);

Confidence_N = Results_N_loc(:,9);
Confidence_AP = Results_AP_loc(:,9);

for i = -4:4
    Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
    length(Confidence_AP(Confidence_AP == 4))
    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
end
for i = 1:9
    if i == 1
    Cumulative_NCounts(i) = Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N_loc);
    Cumulative_APCounts(i) = Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP_loc);
    else
    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Results_N_loc);
    Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
    Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Results_AP_loc);
    end
end

Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

subplot(1,3,a),plot(Cumulative_FA, Cumulative_Hit, 'b-o','LineWidth',1), xlabel('FA'),ylabel('Hit'),title(['Critical objects discrimination in',...
    locType(a),'AUC = ', num2str(AUC)]);
axis square

Confidence_N = [];
Confidence_AP = [];
Results_N_loc = [];
Results_AP_loc = [];
Cumulative_NCounts = [];
Cumulative_APCounts = [];
Cumulative_Hit = [];
Cumulative_FA = [];
AUC = [];
end


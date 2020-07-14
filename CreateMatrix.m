%% locate image files
% clear all;

folderI = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\incongruent patch";
filePatternI = fullfile(folderI, '*.jpg');
theFilesI = struct2cell(dir(filePatternI));
selectnameI(1:1260,1) = string(theFilesI(1,:));

folderC = "C:\Users\liang\OneDrive\Documents\honours\research project\Experiment\New photos\congruent patch";
filePatternC = fullfile(folderC, '*.jpg');
theFilesC = struct2cell(dir(filePatternC));
selectnameC(1:1260,1) = string(theFilesC(1,:));

folderS = 'D:\New photos\nishimoto patch';
filePatternS = fullfile(folderS, '*.jpg');
theFilesS = struct2cell(dir(filePatternS));
selectnameS(1:7044,1) = string(theFilesS(1,:));
%%
% Exp_Results = importdata('Pooled Results.mat'); 
% Exp_Results(:,9) = Exp_Results(:,8).* Exp_Results(:,9);
% 
% %% calculate the length of the big matrix
% Results = importdata('Pooled Results.mat');
% for i = 1:80
%     CR = Results(Results(:,2) == i+2 & Results(:,6)==3,:);
%     index(i) = CR(1,7)==1|CR(1,7)==3|CR(1,7)==5|CR(1,7)==7|CR(1,7)==9;
% end
% 
% full_length = [sum(index==0).*4 + sum(index==1).*5] .*240;
% 
% %create empty matrix
% bigmatrix = zeros(full_length,7);
% 
% %% get the matrix
% for q = 3:82   
%       disp(q); 
% % select data of current image
% Results = Exp_Results(Exp_Results(:,2)== q,:);
%     
% Find_N = Results(:,5) ==1;
% 
% % present patch trials -- signal present for hypo 1
% Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
% Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 
% 
% %Congruent trial with Congruent object, and incongruent trial with
% %incongruent object -- signal present for hypo 2
% Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
% Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;
% 
% %Incongruent trial with congruent object, congruent trial with incongruent
% %object -- signal absent for hypo 2
% Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
% Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
% 
% %select all Present patches for this image
% Present = Results(Find_CAP|Find_IAP,:);
% location = unique(Present(:,7));
% 
%     for a = 1:length(location)
% 
%         mat = zeros(240,7);
%         % select critical absent patch, compare it with responses to 
%         % present patch of current location
%             Confidence_AP = Results(Find_CAP & Results(:,7)==location(a),9);
%             Confidence_N = Results(Find_Congruent_IP,9);
%             num_cong = size(Confidence_N,1);
%             
%             for i = -4:4
%                 Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%                 Confidence_NCounts(i+5) = sum(Confidence_N == -i);
%             end
%             for i = 1:9
%                 if i == 1
%                 Cumulative_NCounts(i) = Confidence_NCounts(i);
%                 Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
%                 Cumulative_APCounts(i) = Confidence_APCounts(i);
%                 Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
%                 else
%                 Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%                 Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
%                 Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%                 Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
%                 end
%             end
% 
%             Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
%             Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
%             
%             mat(1:num_cong,3) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
%             mat(1:num_cong,1) = q-2;
%             
%             % record present and absent patch location
%             mat(1:num_cong,4) = location(a);
%             mat(1:num_cong,5) = unique(Results(Find_Congruent_IP,7));
%             
%             
%             present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(location(a)),'_']));
%             present_name = char(fullfile(folderC,present_name));
%             present_patch = imread(present_name);
% 
%             absent_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
%             absent_name = char(fullfile(folderI,absent_name));
%             absent_patch = imread(absent_name);
% 
%            present_patch = imresize(present_patch,[147 147]);
%            absent_patch = imresize(absent_patch,[147 147]);
% 
%            PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
%            mat(1:15,2) = abs(PD);
%             % record comparison type, critical present vs critical absent in congruent conidition =
%             % 1, incongruent = 2, others = 0.
%             mat(1:num_cong,6) = 1;
%             % record decision x confidence
%             mat(1:num_cong,7) = (abs(Results(Find_Congruent_IP,9))-0.5).*Results(Find_Congruent_IP,8);
%            
%             
%             
%          clear Confidence_APCounts
%          clear Confidence_NCounts
%          clear Confidence_N
%        
% 
%             Confidence_AP = Results(Find_IAP & Results(:,7)==location(a),9);
%             Confidence_N = Results(Find_Incongruent_CP,9);
%             
%             for i = -4:4
%                 Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%                 Confidence_NCounts(i+5) = sum(Confidence_N == -i);
%             end
%             for i = 1:9
%                 if i == 1
%                 Cumulative_NCounts(i) = Confidence_NCounts(i);
%                 Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
%                 Cumulative_APCounts(i) = Confidence_APCounts(i);
%                 Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
%                 else
%                 Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%                 Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
%                 Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%                 Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
%                 end
%             end
% 
%             Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
%             Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
%             mat(num_cong+1:15,3) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
%             mat(num_cong+1:15,1) = q-2;
% 
%             % record present and absent patch location
%             mat(num_cong+1:15,4) = location(a);
%             mat(num_cong+1:15,5) = unique(Results(Find_Incongruent_CP,7));
%             mat(num_cong+1:15,6) = 2;
%             mat(num_cong+1:15,7) = (abs(Results(Find_Incongruent_CP,9))-0.5).*Results(Find_Incongruent_CP,8);
% 
%             
%          clear Confidence_APCounts
%          clear Confidence_NCounts
%          clear Confidence_N
%        
%          
%        % find responses to all null patches, and make pair-wise comparison
%        
%             Absent = Results(Find_N,:);
% 
%             for p = 1:225
% 
%                 Confidence_N = Absent(p,9);
% 
%                 for i = -4:4
%                     Confidence_APCounts(i+5) = sum(Confidence_AP == -i);
%                     Confidence_NCounts(i+5) = sum(Confidence_N == -i);
%                 end
%                 for i = 1:9
%                     if i == 1
%                     Cumulative_NCounts(i) = Confidence_NCounts(i);
%                     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
%                     Cumulative_APCounts(i) = Confidence_APCounts(i);
%                     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
%                     else
%                     Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
%                     Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
%                     Cumulative_APCounts(i)= Cumulative_APCounts(i-1)+ Confidence_APCounts(i);
%                     Cumulative_Hit(i) = Cumulative_APCounts(i)/length(Confidence_AP);
%                     end
%                 end
% 
%                 Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
%                 Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
%                 mat(p+15,3) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
%                 mat(p+15,1) = q-2;
%                 mat(p+15,4) = location(a);
%                 mat(p+15,5) = Absent(p,7);
%                 mat(p+15,6) = 0;
%                 mat(p+15,7) = (abs(Absent(p,9))-0.5).*Absent(p,8);
%                 
%             % calculate patch difference
%             present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(location(a)),'_']));
%             present_name = char(fullfile(folderC,present_name));
%             present_patch = imread(present_name);
% 
%             absent_name = selectnameS(Absent(p,12),:);
%             absent_name = char(fullfile(folderS,absent_name));
%             absent_patch = imread(absent_name);
% 
%            present_patch = imresize(present_patch,[147 147]);
%            absent_patch = imresize(absent_patch,[147 147]);
% 
%            PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
%            mat(p+15,2) = abs(PD);
%             end
%            
%             % record all AUC data for present patch in current location in the
%             % big matrix
%             
%             combined_mat((a-1).*240 + 1: a.*240,:) = mat;       
%             
%     end
%     
% 
% 
% if q == 3
%     bigmatrix(1:length(combined_mat),:) = combined_mat;
%     bigmatrix_index = length(combined_mat);
% else
%     bigmatrix(bigmatrix_index + 1: bigmatrix_index + length(combined_mat),:) = combined_mat;
%     bigmatrix_index = bigmatrix_index + length(combined_mat);
% end
% 
% clear Results
% clear Present
% clear Absent
% 
% clear Confidence_APCounts
% clear Confidence_NCounts
% clear Confidence_N
% clear combined_mat
% 
% end
% 
% data1 = table(bigmatrix(:,1),bigcolumn(:,1),bigmatrix(:,3),bigmatrix(:,4),bigmatrix(:,5),bigcolumn(:,2),bigmatrix(:,7),'VariableNames',{'Image','RGB','AUC','PresentLoc','AbsentLoc','PatchType','Decision'});
%  save_path = ['Second order data matrix.mat'];
%  save(save_path,'data1','-mat');
% 
% %% correcting columns
% 
% Results = importdata('Pooled Results.mat');
% for i = 1:80
%     CR = Results(Results(:,2) == i+2 & Results(:,6)==3,:);
%     index(i) = CR(1,7)==1|CR(1,7)==3|CR(1,7)==5|CR(1,7)==7|CR(1,7)==9;
% end
% 
% full_length = [sum(index==0).*4 + sum(index==1).*5] .*240;
% 
% %create empty matrix
% bigcolumn = zeros(full_length,2);
% 
% for q = 3:82   
%       disp(q); 
% % select data of current image
% Results = Exp_Results(Exp_Results(:,2)== q,:);
%     
% Find_N = Results(:,5) ==1;
% 
% % present patch trials -- signal present for hypo 1
% Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
% Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 
% 
% %Congruent trial with Congruent object, and incongruent trial with
% %incongruent object -- signal present for hypo 2
% Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
% Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;
% 
% %Incongruent trial with congruent object, congruent trial with incongruent
% %object -- signal absent for hypo 2
% Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
% Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
% 
% %select all Present patches for this image
% Present = Results(Find_CAP|Find_IAP,:);
% location = unique(Present(:,7));
% 
%     for a = 1:length(location)
% 
%         mat = zeros(240,2);
%         % select critical absent patch, compare it with responses to 
%         % present patch of current location
%             
%             
%             present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(location(a)),'_']));
%             present_name = char(fullfile(folderC,present_name));
%             present_patch = imread(present_name);
% 
%             absent_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
%             absent_name = char(fullfile(folderI,absent_name));
%             absent_patch = imread(absent_name);
% 
%            present_patch = imresize(present_patch,[147 147]);
%            absent_patch = imresize(absent_patch,[147 147]);
% 
%            PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
%            mat(1:15,1) = abs(PD);
%            
%            if location(a) == unique(Results(Find_Congruent_IP,7))
%                mat(1:15,2) = 2; % critical object comparison
%            else
%                mat(1:15,2) = 1; % critical absent vs non-critical present
%            end
%             % record comparison type, critical present vs critical absent in congruent conidition =
%             % 1, incongruent = 2, others = 0.           
%             
%             
%          clear Confidence_APCounts
%          clear Confidence_NCounts
%          clear Confidence_N
%        
%          
%        % find responses to all null patches, and make pair-wise comparison
%        
%             Absent = Results(Find_N,:);
% 
%             for p = 1:225
% 
%                 
%             % calculate patch difference
%             present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(location(a)),'_']));
%             present_name = char(fullfile(folderC,present_name));
%             present_patch = imread(present_name);
% 
%             absent_name = selectnameS(Absent(p,12),:);
%             absent_name = char(fullfile(folderS,absent_name));
%             absent_patch = imread(absent_name);
% 
%            present_patch = imresize(present_patch,[147 147]);
%            absent_patch = imresize(absent_patch,[147 147]);
% 
%            PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
%            mat(p+15,1) = abs(PD);
% 
%             end
%            
%             % record all AUC data for present patch in current location in the
%             % big matrix
%             
%             combined_mat((a-1).*240 + 1: a.*240,:) = mat;       
%             
%     end
%     
% 
% 
% if q == 3
%     bigcolumn(1:length(combined_mat),:) = combined_mat;
%     bigcolumn_index = length(combined_mat);
% else
%     bigcolumn(bigcolumn_index + 1: bigcolumn_index + length(combined_mat),:) = combined_mat;
%     bigcolumn_index = bigcolumn_index + length(combined_mat);
% end
% 
% clear Results
% clear Present
% clear Absent
% 
% clear Confidence_APCounts
% clear Confidence_NCounts
% clear Confidence_N
% clear combined_mat
% 
% end

% bigcolumn = zeros(full_length,1);
% Exp_Results = importdata('Pooled Results.mat'); 
% Exp_Results(:,9) = Exp_Results(:,8).* (Exp_Results(:,9)-0.5);
% 
% for q = 3:82  
%   
% % select data of current image
% 
% Results = Exp_Results(Exp_Results(:,2)== q,:);
%     
% Find_N = Results(:,5) ==1;
% 
% % present patch trials -- signal present for hypo 1
% Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; 
% Find_IAP = Results(:,4) == 1 & Results(:,5) == 3; 
% 
% %Congruent trial with Congruent object, and incongruent trial with
% %incongruent object -- signal present for hypo 2
% Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
% Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;
% 
% %Incongruent trial with congruent object, congruent trial with incongruent
% %object -- signal absent for hypo 2
% Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
% Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;
% 
% %select all Present patches for this image
% Present = Results(Find_CAP|Find_IAP,:);
% location = unique(Present(:,7));
% 
%     for a = 1:length(location)
% 
%         mat = zeros(226,1);
% 
%             % record comparison type, critical present vs critical absent =
%             % 2, others = 1.
%         
%        mat(1,1) = Results(:,9);
%         % calcualte and record patch difference
% 
% 
%          clear Confidence_APCounts
%          clear Confidence_NCounts
%          clear Confidence_N
%          clear Absent
%         
%          
%        % find responses to all null patches, and make pair-wise comparison
%        Absent = Results(Find_N,:);
%        for p = 1:225    
%            mat(p+1,1) = Results(p,9);
%        end    
%             % record all AUC data for present patch in current location in the
%             % big matrix
%             
%             combined_mat((a-1).*226 + 1: a.*226,:) = mat;       
%             
%     end
%     
% if q == 3
%     bigcolumn(1:length(combined_mat),:) = combined_mat;
%     bigcolumn_index = length(combined_mat);
% else
%     bigcolumn(bigcolumn_index + 1: bigcolumn_index + length(combined_mat),:) = combined_mat;
%     bigcolumn_index = bigcolumn_index + length(combined_mat);
% end
% 
% clear Results
% clear Present
% clear Absent
% 
% clear Confidence_APCounts
% clear Confidence_NCounts
% clear Confidence_N
% clear combined_mat
% 
% end
% 
% bigmatrix = [bigmatrix bigcolumn];
% % 
% % save_path = ['AUC and RGB matrix.mat'];
% % save(save_path,'data','-mat');
% 
% 
%% take absent patch decision only

Exp_Results = importdata('Pooled Results.mat');
Exp_Results(:,9) = (Exp_Results(:,9)-0.5).*Exp_Results(:,8);

full_length = (1+15).*80.*15;

%create empty matrix
matrix = zeros(full_length,1);


for q = 3:82   
      disp(q); 
    % select data of current image
Results = Exp_Results(Exp_Results(:,2)== q,:);

%patch index
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

Find_N = Results(:,5) ==1;
    
% create empty matrix
mat = zeros(240,1);
     
%  select critical absent patch, compare it with responses to 
        % present patch of current location
        
%             Confidence_N = Results(Find_Congruent_IP,9);
            sub_cong = Results(Find_Congruent_IP,1);
            num_cong = size(sub_cong,1);
            
            mat(1:num_cong,1) = sub_cong;
            
            % record present and absent patch location
%             mat(1:15,3) = unique(Results(Find_Congruent_IP,7));
%             mat(1:num_cong,4) = Confidence_N;
            
            
%             present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
%             present_name = char(fullfile(folderC,present_name));
%             present_patch = imread(present_name);
% 
%             absent_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
%             absent_name = char(fullfile(folderI,absent_name));
%             absent_patch = imread(absent_name);
% 
%            present_patch = imresize(present_patch,[147 147]);
%            absent_patch = imresize(absent_patch,[147 147]);
% 
%            PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
%            mat(1:15,2) = abs(PD);
           
%             % record comparison type, critical present vs critical absent in congruent conidition =
%             % 1, incongruent = 2, others = 0.
%             mat(1:num_cong,5) = 1;
%      
%             
%            % record eccentricity 
%            if unique(Results(Find_Congruent_IP,7))== 1||unique(Results(Find_Congruent_IP,7))==3||unique(Results(Find_Congruent_IP,7))==7||unique(Results(Find_Congruent_IP,7))==9
%                mat(1:15,6) = 9.20;
%            elseif unique(Results(Find_Congruent_IP,7)) == 2||unique(Results(Find_Congruent_IP,7))==4||unique(Results(Find_Congruent_IP,7))==6||unique(Results(Find_Congruent_IP,7))==8
%                mat(1:15,6) = 6.50;
%            else 
%                mat(1:15,6)= 0;
%            end
            
            
         clear Confidence_APCounts
         clear Confidence_NCounts
         clear Confidence_N
         clear present_name
         clear absent_name
         clear present_patch
         clear absent_patch       

%             Confidence_N = Results(Find_Incongruent_CP,9);
            sub_incong = Results(Find_Incongruent_CP,1);
            mat(num_cong+1:15,1) = sub_incong;

            
         clear Confidence_APCounts
         clear Confidence_NCounts
         clear Confidence_N
   

% find responses to all null patches


        Absent_current = Results(Find_N,:);
        
        %select absent patch of current location
       
       for p = 1:size(Absent_current,1)

                mat(p+15,1) = Absent_current(p,1);
%                 mat(p+15,3) = Absent_current(p,7);
%                 mat(p+15,4) = Absent_current(p,9);
%                 
%             % calculate patch difference
%             
%             if Absent_current(p,4) == 0
%             present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(Absent_current(p,7)),'_']));
%             present_name = char(fullfile(folderC,present_name));
%             present_patch = imread(present_name);
%             else
%             present_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(Absent_current(p,7)),'_']));
%             present_name = char(fullfile(folderI,present_name));
%             present_patch = imread(present_name);
%             end
%             
%             absent_name = selectnameS(Absent_current(p,12),:);
%             absent_name = char(fullfile(folderS,absent_name));
%             absent_patch = imread(absent_name);
% 
%            present_patch = imresize(present_patch,[147 147]);
%            absent_patch = imresize(absent_patch,[147 147]);
% 
%            PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
%            mat(p+15,2) = abs(PD);
%            
%            %record patch type
%            mat(p+15,5) = 0;
%            
%            % record eccentricity
%            
%            
%            if Absent_current(p,7) == 1||Absent_current(p,7)==3||Absent_current(p,7)==7||Absent_current(p,7)==9
%                mat(p+15,6) = 9.20;
%            elseif Absent_current(p,7) == 2||Absent_current(p,7)==4||Absent_current(p,7)==6||Absent_current(p,7)==8
%                mat(p+15,6) = 6.50;
%            else 
%                mat(p+15,6)= 0;
%            end
%            
       end
%       

      
            clear Absent_current

    


if q == 3
    matrix(1:size(mat,1),:) = mat;
    matrix_index = size(mat,1);
else
    matrix(matrix_index + 1: matrix_index + size(mat,1),:) = mat;
    matrix_index = matrix_index + size(mat,1);
end


clear Results
clear Absent


clear Confidence_APCounts
clear Confidence_NCounts
clear Confidence_N
clear present_name
clear absent_name
clear present_patch
clear absent_patch

end

% data_allab = table(matrix(:,1),matrix(:,2),matrix(:,3),matrix(:,4),matrix(:,5),matrix(:,6),'VariableNames',{'Image','RGB','AbsentLoc','Decision','PatchType','Eccentricity'});
data_allab = addvars(data_allab,matrix, 'NewVariableNames',{'Subject'});
save_path = 'All absent decisions.mat';
save(save_path,'data_allab','-mat');

%% this calculates the luminance and contrast

Exp_Results = importdata('Pooled Results.mat');
Exp_Results(:,9) = (Exp_Results(:,9)-0.5).*Exp_Results(:,8);

full_length = (1+15).*80.*15;

%create empty matrix
matrix = zeros(full_length,2);


for q = 3:82   
      disp(q); 
    % select data of current image
Results = Exp_Results(Exp_Results(:,2)== q,:);

%patch index
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

Find_N = Results(:,5) ==1;
    
% create empty matrix
mat = zeros(240,2);


            % get modified patches, calculate luminance and contrast
            % difference

            present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
            present_name = char(fullfile(folderC,present_name));
            present_patch = imread(present_name);

            absent_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
            absent_name = char(fullfile(folderI,absent_name));
            absent_patch = imread(absent_name);

           present_patch = imresize(present_patch,[147 147]);
           absent_patch = imresize(absent_patch,[147 147]);

           % luminance
           PD = mean2(rgb2gray(present_patch))-mean2(rgb2gray(absent_patch));
           mat(1:15,1) = abs(PD);
           
           % contrast
          
           contrast_absent = (rgb2gray(absent_patch)-min(rgb2gray(absent_patch)))./...
               (max(rgb2gray(absent_patch))-min(rgb2gray(absent_patch)));
           contrast_present = (rgb2gray(present_patch)-min(rgb2gray(present_patch)))./...
               (max(rgb2gray(present_patch))-min(rgb2gray(present_patch)));
           contrast_difference = sum(sum(contrast_absent))-sum(sum(contrast_present));
           mat(1:15,2) =abs(contrast_difference);

    
         clear present_name
         clear absent_name
         clear present_patch
         clear absent_patch       
         clear PD
         clear contrast_difference
         clear contrast_absent
         clear contrast_present
 
% find responses to all null patches


        Absent_current = Results(Find_N,:);
        
        %select absent patch of current location
       
       for p = 1:size(Absent_current,1)

            % calculate patch difference
            
            if Absent_current(p,4) == 0
            present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(Absent_current(p,7)),'_']));
            present_name = char(fullfile(folderC,present_name));
            present_patch = imread(present_name);
            else
            present_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(Absent_current(p,7)),'_']));
            present_name = char(fullfile(folderI,present_name));
            present_patch = imread(present_name);
            end
            
            absent_name = selectnameS(Absent_current(p,12),:);
            absent_name = char(fullfile(folderS,absent_name));
            absent_patch = imread(absent_name);

           present_patch = imresize(present_patch,[147 147]);
           absent_patch = imresize(absent_patch,[147 147]);
           
           % luminance
           PD = mean2(rgb2gray(present_patch))-mean2(rgb2gray(absent_patch));
           mat(p+15,1) = abs(PD);
           
           % contrast
          
           contrast_absent = (rgb2gray(absent_patch)-min(rgb2gray(absent_patch)))./...
               (max(rgb2gray(absent_patch))-min(rgb2gray(absent_patch)));
           contrast_present = (rgb2gray(present_patch)-min(rgb2gray(present_patch)))./...
               (max(rgb2gray(present_patch))-min(rgb2gray(present_patch)));
           contrast_difference = sum(sum(contrast_absent))-sum(sum(contrast_present));
           mat(p+15,2) =abs(contrast_difference);
           
         clear present_name
         clear absent_name
         clear present_patch
         clear absent_patch       
         clear PD
         clear contrast_difference
         clear contrast_absent
         clear contrast_present
           
       end
      

      
            clear Absent_current

    


if q == 3
    matrix(1:size(mat,1),:) = mat;
    matrix_index = size(mat,1);
else
    matrix(matrix_index + 1: matrix_index + size(mat,1),:) = mat;
    matrix_index = matrix_index + size(mat,1);
end


clear Results
clear Absent


clear Confidence_APCounts
clear Confidence_NCounts
clear Confidence_N
clear present_name
clear absent_name
clear present_patch
clear absent_patch

end

data_allab = addvars(data_allab,matrix(:,1),matrix(:,2),'NewVariableNames',{'Luminance','Contrast'});
save('All absent decisions.mat','data_allab','-mat');

%% Find patch filename
Exp_Results = importdata('Pooled Results.mat');
Exp_Results(:,9) = (Exp_Results(:,9)-0.5).*Exp_Results(:,8);

full_length = (1+15).*80.*15;

for q = 3:82   
      disp(q); 
    % select data of current image
Results = Exp_Results(Exp_Results(:,2)== q,:);

%patch index
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

Find_N = Results(:,5) ==1;
    
% create empty matrix
patch_type = data_allab.PatchType(data_allab.Image == q-2);
patch_type = patch_type(1:15,:);

%  select critical absent patch, compare it with responses to 
            

            current_present_list(1:sum(patch_type == 1),1) = selectnameC(contains(selectnameC,['_', num2str(q), '_',...
                num2str(unique(Results(Find_Congruent_IP,7))),'_']));
            current_present_list(sum(patch_type==1)+1:15,1) = selectnameI(contains(selectnameI,['_', num2str(q), '_',...
                num2str(unique(Results(Find_Congruent_IP,7))),'_']));

            current_absent_list(1:sum(patch_type == 1),1) = selectnameI(contains(selectnameI,['_', num2str(q), '_',...
                num2str(unique(Results(Find_Congruent_IP,7))),'_']));
            current_absent_list(sum(patch_type==1)+1:15,1) = selectnameC(contains(selectnameC,['_', num2str(q), '_',...
                num2str(unique(Results(Find_Congruent_IP,7))),'_'])); 


% find all null patches


        Absent_current = Results(Find_N,:);
        
        %select absent patch of current location
       
       for p = 1:size(Absent_current,1)

               
%            % record present patch corresponding to the null patch
               if Absent_current(p,4) == 0
                   current_present_list(p+15,1) = selectnameC(contains(selectnameC,['_', num2str(q), '_',...
                num2str(unique(Results(Find_Congruent_IP,7))),'_']));
               else
                   current_present_list(p+15,1,:) = selectnameI(contains(selectnameI,['_', num2str(q), '_',...
                num2str(unique(Results(Find_Congruent_IP,7))),'_']));
               end
%             
              % record null patch name  
%             
                current_absent_list(p+15,1)= selectnameS(Absent_current(p,12),:);
%            
       end
%       

      
            clear Absent_current

    


if q == 3
    present_list = current_present_list;
    absent_list = current_absent_list;
else
    present_list = [present_list;current_present_list];
    absent_list = [absent_list;current_absent_list];

end

clear current_absent_list
clear current_present_list
clear patch_type

end
save('all absent patch name.mat','absent_list','-mat');
save('all present patch name.mat','present_list','-mat');

%% calculate saliency

Saliency = mean([zscore(data_allab.RGB) zscore(data_allab.Contrast) zscore(data_allab.Luminance)],2);

data_allab = addvars(data_allab,Saliency,'After','Contrast','NewVariableNames',{'Saliency'});

% patchtype = 0, min(null patch Saliency - present patch saliency),
% Eccentricity == 9.2; 
save('All absent decisions.mat','data_allab','-mat');

%% record congruency

Exp_Results = importdata('Pooled Results.mat');
Exp_Results(:,9) = (Exp_Results(:,9)-0.5).*Exp_Results(:,8);

full_length = (1+15).*80.*15;

for q = 3:82   
      disp(q); 
    % select data of current image
Results = Exp_Results(Exp_Results(:,2)== q,:);

%patch index
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

Find_N = Results(:,5) ==1;
    
% create empty matrix
patch_type = data_allab.PatchType(data_allab.Image == q-2);
patch_type = patch_type(1:15,:);

%  select critical absent patch, compare it with responses to 
            
modified_congruency = patch_type == 2;


% find all null patches


null_congruency = Results(Find_N,4);

if q == 3
    Congruency = [modified_congruency; null_congruency];
else
    Congruency = [Congruency; modified_congruency; null_congruency];
end

clear null_congruency
clear modified_congruency

end

data_allab = addvars(data_allab,Congruency,'After','Image','NewVariableNames','Congruency');

%% this calculates differences in R, G, B

Exp_Results = importdata('Pooled Results.mat');
Exp_Results(:,9) = (Exp_Results(:,9)-0.5).*Exp_Results(:,8);

full_length = (1+15).*80.*15;

%create empty matrix
matrix = zeros(full_length,3);


for q = 3:82   
      disp(q); 
    % select data of current image
Results = Exp_Results(Exp_Results(:,2)== q,:);

%patch index
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

Find_N = Results(:,5) ==1;
    
% create empty matrix
mat = zeros(240,3);

    

            % get modified patches, calculate luminance and contrast
            % difference

            present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
            present_name = char(fullfile(folderC,present_name));
            present_patch = imread(present_name);

            absent_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(unique(Results(Find_Congruent_IP,7))),'_']));
            absent_name = char(fullfile(folderI,absent_name));
            absent_patch = imread(absent_name);

           present_patch = imresize(present_patch,[147 147]);
           absent_patch = imresize(absent_patch,[147 147]);
           
     for layer = 1:3
         
           % differences on each layer
           PD = sum(sum(present_patch(:,:,layer)-absent_patch(:,:,layer)));
           mat(1:15,layer) = abs(PD);

     end
         clear present_name
         clear absent_name
         clear present_patch
         clear absent_patch       
         clear PD
         clear contrast_difference
         clear contrast_absent
         clear contrast_present
 
% find responses to all null patches


        Absent_current = Results(Find_N,:);
        
        %select absent patch of current location
       
       for p = 1:size(Absent_current,1)

            % load patches
            
            if Absent_current(p,4) == 0
            present_name = selectnameC(contains(selectnameC,['_', num2str(q), '_',num2str(Absent_current(p,7)),'_']));
            present_name = char(fullfile(folderC,present_name));
            present_patch = imread(present_name);
            else
            present_name = selectnameI(contains(selectnameI,['_', num2str(q), '_',num2str(Absent_current(p,7)),'_']));
            present_name = char(fullfile(folderI,present_name));
            present_patch = imread(present_name);
            end
            
            absent_name = selectnameS(Absent_current(p,12),:);
            absent_name = char(fullfile(folderS,absent_name));
            absent_patch = imread(absent_name);

           present_patch = imresize(present_patch,[147 147]);
           absent_patch = imresize(absent_patch,[147 147]);
           
           % calculate patch difference on each of R, G and B
           for layer = 1:3
               PD = sum(sum(present_patch(:,:,layer)-absent_patch(:,:,layer)));
               mat(p+15,layer) = abs(PD);
           end
           
         clear present_name
         clear absent_name
         clear present_patch
         clear absent_patch       
         clear PD
         clear contrast_difference
         clear contrast_absent
         clear contrast_present
           
       end
      
    
      
            clear Absent_current

    


if q == 3
    matrix(1:size(mat,1),:) = mat;
    matrix_index = size(mat,1);
else
    matrix(matrix_index + 1: matrix_index + size(mat,1),:) = mat;
    matrix_index = matrix_index + size(mat,1);
end


clear Results
clear Absent


clear Confidence_APCounts
clear Confidence_NCounts
clear Confidence_N
clear present_name
clear absent_name
clear present_patch
clear absent_patch

end

 data_allab = addvars(data_allab, matrix(:,1),matrix(:,2),matrix(:,3),'After','Saliency','NewVariableNames',{'R','G','B'});
% 
 save('All absent decisions.mat','data_allab','-mat');

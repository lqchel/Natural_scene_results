clear all;

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



%% calculate the mean AUC for each image, find out the min and max

for q = 3:82
    
Results = importdata('Pooled Results.mat');  
Results(:,9) = Results(:,8).* Results(:,9);
Results = Results(Results(:,2)== q,:);
    
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

Present = Results(Find_CAP|Find_IAP,:);

location = unique(Present(:,7));
% RGB = unique(Absent(:,11));
% RGB = sort(RGB,'ascend');

    for a = 1:length(location)

        mat = zeros(226,1);

            Absent = Results(Find_Congruent_IP|Find_Incongruent_CP,:);
            Confidence_AP = Present(Present(:,7)==location(a),9);
            Confidence_N = Absent(:,9);
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
            mat(1,1) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);


         clear Confidence_APCounts
         clear Confidence_NCounts
         clear Confidence_N
         clear Absent

            Absent = Results(Find_N,:);

            for p = 1:225

                Confidence_N = Absent(p,9);

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
                mat(p+1,1) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

            end
            combined_mat((a-1).*226 + 1: a.*226,:) = mat;       

    end
    
mAUC(q-2,1) = mean(combined_mat);

clear Results
clear Present
clear Absent

clear Confidence_APCounts
clear Confidence_NCounts
clear Confidence_N

end

[M,I] = min(mAUC,[],1,'linear');
imgN(1) = I+2;
[M,I] = sort(mAUC,'ascend');
imgN(2) = find(I == 40,1) + 2;
[M,I] = max(mAUC,[],1,'linear');
imgN(3) = I+2;

% 
% histogram(mAUC);
% m = round(mean(mAUC),2);
% title(['Distribution of Image AUC, mean = ', num2str(m)]);


clear combined_mat
clear Results

colours = cbrewer('qual', 'Set1', 8); 
for n = 1:3
    
imgNum = imgN(n);
Results = importdata('Pooled Results.mat');
Results(:,9) = Results(:,8).* Results(:,9);

Results = Results(Results(:,2)== imgNum,:);

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


Present = Results(Find_CAP|Find_IAP,:);

location = unique(Present(:,7));
% RGB = unique(Absent(:,11));
% RGB = sort(RGB,'ascend');


for a = 1:length(location)
    
    mat = zeros(226,3);
    Absent = Results(Find_Congruent_IP|Find_Incongruent_CP,:);
        Confidence_AP = Present(Present(:,7)==location(a),9);
        Confidence_N = Absent(:,9);
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
        mat(1,1) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
        
        if Absent(1,7) == location(a)
            loc_sq = a;
            crit_loc = (a-1).*226+1;
        end
        
        % calcualte patch difference
            present_name = selectnameC(contains(selectnameC,['_', num2str(imgNum), '_',num2str(location(a)),'_']));
            present_name = char(fullfile(folderC,present_name));
            present_patch = imread(present_name);

            absent_name = selectnameI(contains(selectnameI,['_', num2str(imgNum), '_',num2str(Absent(1,7)),'_']));
            absent_name = char(fullfile(folderI,absent_name));
            absent_patch = imread(absent_name);

           present_patch = imresize(present_patch,[147 147]);
           absent_patch = imresize(absent_patch,[147 147]);

           PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
           mat(1,2) = abs(PD);
           
           if length(location) == 4
               mat(1,3) = 6.18;
           elseif length(location) == 5 && location(a) == 5
               mat(1,3) = 0;
           else
               mat(1,3) = 9.16;
           end
           
     clear Confidence_APCounts
     clear Confidence_NCounts
     clear Confidence_N
     clear PD
     clear Absent
        
        Absent = Results(Find_N,:);

        for p = 1:225

            Confidence_N = Absent(p,9);

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
            mat(p+1,1) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

            %calculate patch difference
                present_name = selectnameC(contains(selectnameC,['_', num2str(imgNum), '_',num2str(location(a)),'_']));
                present_name = char(fullfile(folderC,present_name));
                present_patch = imread(present_name);
          
                absent_name = selectnameS(Absent(p,12),:);
                absent_name = char(fullfile(folderS,absent_name));
                absent_patch = imread(absent_name);

               present_patch = imresize(present_patch,[147 147]);
               absent_patch = imresize(absent_patch,[147 147]);

               PD = sum(sum(sum(double(present_patch) - double(absent_patch))));
               mat(p+1,2) = abs(PD);
               mat(p+1,3) = location(a);

        end

       
        [B,I] = sort(mat(2:end, 2),'ascend');
        d = mat(2:end,1);
        mat(2:end,1) = d(I);
        mat(2:end,2) = B;
        combined_mat((a-1).*226 + 1: a.*226,:) = [mat reshape(0:1:225,[226,1])];       
        
end
AUC = mean(combined_mat(:,1));

%% get a heatmap
% data = table(combined_mat(:,1),combined_mat(:,3),combined_mat(:,4),'VariableNames',{'AUC','location','RGB'});
% subplot(2,1,2), h = heatmap(data, 'RGB', 'location','ColorVariable','AUC');
% h.XLabel = 'RGB Difference';
% h.YLabel = 'Present Patch Location';
% % ax = gca;
% % ax.Color = 'white';
% title(['individual present-absent patch comparison AUC for image ', imgNum]);

%% get a scatterplot

if a == 4 
    
% subplot(1,3,n), h = gscatter(combined_mat(:,2),combined_mat(:,1),combined_mat(:,3),[colours(1:2,:);colours(5,:);colours(8,:)],'.',[12 12]);
subplot(1,3,n), h = plot(combined_mat(:,2),combined_mat(:,1),'.','MarkerEdgeColor',colours(5,:),'MarkerFaceColor',colours(5,:),'MarkerSize',10);
xlabel('RGB Difference');
ylabel('AUC');
title(['Patch AUC for image ', num2str(imgNum),' mean = ', num2str(round(AUC,2))]);
xlim([0 15.*10^6]);
ylim([0 1]);
hold on
subplot(1,3,n),plot(combined_mat(crit_loc,2),combined_mat(crit_loc,1),'s','MarkerEdgeColor',colours(6,:),'LineWidth',2,'MarkerFaceColor',colours(5,:),...
    'MarkerSize',12);
hold off
legend({'Eccentricity = 6.18','Critical objects'});
   
else
    
subplot(1,3,n), h = gscatter(combined_mat(:,2),combined_mat(:,1),combined_mat(:,3),colours(1:2,:),'.',[12 12]);
xlabel('RGB Difference');
ylabel('AUC');
title(['Patch AUC for image ', num2str(imgNum),' mean = ', num2str(round(AUC,2))]);
xlim([0 15.*10^6]);
ylim([0 1]);
legend({'Eccentricity = 0','Eccentricity = 9.16'});
hold on

if combined_mat(crit_loc,3) == 0
    loc_sq = 1;
else
    loc_sq = 2;
end

subplot(1,3,n),plot(combined_mat(crit_loc,2),combined_mat(crit_loc,1),'s','MarkerEdgeColor',colours(6,:),'LineWidth',2,'MarkerFaceColor',colours(loc_sq,:),...
    'MarkerSize',12,'DisplayName','Critical objects');

hold off

end

clear crit_loc
clear loc_sq
clear Absent
clear Results
clear combined_mat

end

subplot(2,1,2), plot(1:80,mAUC,'bo-');
ylim([0.5 1]), xlabel('Image number'), ylabel('AUC');
title('AUC of image 1-80, new calculation method');




%% compute bits per second

Idv_bits = zeros(80,1);
AUC = importdata('AUC and RGB matrix.mat');
index = importdata('FourOrFive.mat');

for d = 1:80
    accuracy = mean(AUC(AUC(:,1)==d,3));
    
    if index(d)== 0
    Idv_bits(d,1) = [(226 + 4) .* log2(2.* accuracy)]/0.133;
    else
    Idv_bits(d,1) = [(226 + 5) .* log2(2.* accuracy)]/0.133;
    end
end

subplot(1,2,1), histogram(Idv_bits), title(['Information transfer rate, mean = ', num2str(round(mean(Idv_bits),2))]);
xlabel('Bits per second'), ylabel('Frequency');
set(gca,'FontSize',12);
subplot(1,2,2), plot(1:80,Idv_bits,'bo-');
xlabel('Image number'), ylabel('Bits/sec'),title('Information transfer rate of image 1-80');
set(gca,'FontSize',12);

%% new analysis, on null patches of subject 1, for image 3
Results = importdata('Pooled Results.mat');  
Results(:,9) = Results(:,8).* Results(:,9);
Results = Results(Results(:,2)== 3,:);
    
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

%% Congruent condition
 

Present = Results(Find_CAP,:);

modified = Results(Find_Congruent_IP,:);
location = unique(Present(:,7));

mat = zeros(length(location)+1,16);



    for a = 1:length(location)

%first calculate original vs. modified
            Absent = Results(Find_Congruent_IP,:);
            Confidence_AP = Present(Present(:,7)==location(a),9);
            Confidence_N = Absent(:,9);
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
            mat(a,16) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);


         clear Confidence_APCounts
         clear Confidence_NCounts
         clear Confidence_N
         clear Absent
 % then, original vs. null
 
            Absent = Results(Find_N,:);

            for p = 1:15

                Confidence_N = Absent(p,9);

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
                mat(a,p) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

            end  
            
        clear p
         clear Confidence_APCounts
         clear Confidence_NCounts
         clear Confidence_N
         clear Cumulative_NCounts
         clear Cumulative_FA
         clear Cumulative_APCounts
         clear Cumulative_Hit
         
            
   % Finally, modified vs. null
            Confidence_modified = modified(:,9);
            for p = 1:15

                Confidence_N = Absent(p,9);
                
                for i = -4:4
                    Confidence_modCounts(i+5) = sum(Confidence_modified == -i);
                    Confidence_NCounts(i+5) = sum(Confidence_N == -i);
                end
                for i = 1:9
                    if i == 1
                    Cumulative_NCounts(i) = Confidence_NCounts(i);
                    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
                    Cumulative_modCounts(i) = Confidence_modCounts(i);
                    Cumulative_Hit(i) = Cumulative_modCounts(i)/length(Confidence_modified);
                    else
                    Cumulative_NCounts(i) = Cumulative_NCounts(i-1) + Confidence_NCounts(i);
                    Cumulative_FA(i) = Cumulative_NCounts(i)/length(Confidence_N);
                    Cumulative_modCounts(i)= Cumulative_modCounts(i-1)+ Confidence_modCounts(i);
                    Cumulative_Hit(i) = Cumulative_modCounts(i)/length(Confidence_modified);
                    end
                end

                Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
                Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
                mat(a+1,p) = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);

            end  
            
    end

    %% 
% heat_axes = axes('Position',[0.2 0.2 0.7 0.7])
AUC_map = heatmap(mat);
axes = gca;

    
    for  a = 1:length(location)
    present_list_sub1(:,a) = selectnameC(contains(selectnameC,['cong_', num2str(3), '_',num2str(location(a)),'_']));
    
    end
   

null_list_sub1 = [absent_list(data_allab.Image == 1 & data_allab.PatchType == 0 & data_allab.Subject == 1); absent_list(1,:)];
modified_list = selectnameS(contains(selectnameS,['incong_',num2str(3),'_',num2str(unique(modified(:,7)))]));

heat_axes_position = axes.Position;

absent_axes = axes('Position',[0.2 0.1 0.8 0.1],'XLim',[0 8060], 'YLim', [0 450]);
present_axes = axes('Position', [0.1 0.1 0.1 0.3], 'XLim',[0 440], 'YLim',[0 8060]);









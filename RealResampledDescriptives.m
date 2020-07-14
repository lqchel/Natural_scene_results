
 clear all
% load data
load('All absent decisions.mat');
yes_no = data_allab.Decision>0;
yes_no = yes_no - (yes_no==0);
data_allab.Decision = abs(data_allab.Decision)+0.5;
data_allab.Decision = data_allab.Decision .* yes_no;

location = [0 6.50 9.20];
decision_option = [-4 -3 -2 -1 1 2 3 4];


%% luminance
for e = 1:3 % uncomment for raw distribution data

% create histogram of modified patch Luminance
h = histogram(data_allab.Luminance(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);

% bin up the null patch responses, based on the modified patch histogram
Null_bin = cell(1,h.NumBins);

for i = 1:h.NumBins
Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.Luminance >= h.BinEdges(i)& data_allab.Luminance <= h.BinEdges(i+1)...
    & data_allab.Eccentricity == location(e),:);
end

% bootstrapping

for n = 1:200
    sampled_trials = zeros(sum(h.NumBins .* h.Values),width(data_allab));
    for c = 1:h.NumBins
        if h.Values(c) ~= 0
            rand_index = randperm(height(Null_bin{:,c})); 
            old_sq = Null_bin{:,c};
            new_sq = old_sq(rand_index,:);
            selected = new_sq(1:h.Values(:,c),:);
            
            if c == 1
                sampled_trials = selected;
            else
                sampled_trials = [sampled_trials; selected];
            end
        end
    end
    
    for sub = 1:15
        subject_trials = sampled_trials(sampled_trials.Subject == sub,:);
     
        for i = 1:8
            sampled_matrix(sub,i) = size(subject_trials(subject_trials.Decision == decision_option(i),:),1)/...
            length(subject_trials.Decision);
        end
        clear subject_trials;
   
    end
    all_matrix(:,:,n) = sampled_matrix;
    clear sampled_matrix;
end
lum_locmeans(:,:,e) = mean(all_matrix,3);

clear h;
clear all_matrix;
end
clear e;

% separating eccentricity
% 
% lum_locmedians = zeros(200,3);
% 
% for e = 1:3
% 
% % create histogram of modified patch Luminance
% h = histogram(data_allab.Luminance(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);
% 
% % bin up the null patch responses, based on the modified patch histogram
% Null_bin = cell(1,h.NumBins);
% 
% for i = 1:h.NumBins
% Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.Luminance >= h.BinEdges(i)& data_allab.Luminance <= h.BinEdges(i+1)...
%     & data_allab.Eccentricity == location(e),:);
% end
% 
% % bootstrapping
% 
% btp_trials = cell(200,20);
% 
% 
% for n = 1:200
%     for c = 1:h.NumBins
%         if h.Values(c) ~= 0
%             rand_index = randperm(height(Null_bin{:,c})); 
%             old_sq = Null_bin{:,c};
%             new_sq = old_sq(rand_index,:);
%             selected = new_sq(1:h.Values(:,c),:);
%             if c == 1
%                 btp_data = selected.Decision;
%             else
%                 btp_data = [btp_data;selected.Decision];
%             end
%         end
%         clear selected
%         clear new_sq
%         clear old_sq
%     end
%     lum_locmedians(n,e) = median(btp_data);
% end    
% clear h;
% end
% clear e;

% convert data histogram into line-plottable
% 
% lum_cross_ecc = mean(lum_locmedians,2);
% d = histogram(lum_cross_ecc,20);
% for i = 1:20
%     for s = 1:d.NumBins
%         
%     end
% end


%% Contrast
for e = 1:3 %% uncomment for raw distribution data

% create histogram of modified patch Luminance
h = histogram(data_allab.Contrast(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);

% bin up the null patch responses, based on the modified patch histogram
Null_bin = cell(1,h.NumBins);

for i = 1:h.NumBins
Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.Contrast >= h.BinEdges(i)& data_allab.Contrast <= h.BinEdges(i+1)...
    & data_allab.Eccentricity == location(e),:);
end

% bootstrapping

for n = 1:200
    sampled_trials = zeros(sum(h.NumBins .* h.Values),width(data_allab));
    for c = 1:h.NumBins
        if h.Values(c) ~= 0
            rand_index = randperm(height(Null_bin{:,c})); 
            old_sq = Null_bin{:,c};
            new_sq = old_sq(rand_index,:);
            selected = new_sq(1:h.Values(:,c),:);
            
            if c == 1
                sampled_trials = selected;
            else
                sampled_trials = [sampled_trials; selected];
            end
        end
    end
    
    for sub = 1:15
        subject_trials = sampled_trials(sampled_trials.Subject == sub,:);
        for i = 1:8
            sampled_matrix(sub,i) = size(subject_trials(subject_trials.Decision == decision_option(i),:),1)/...
            length(subject_trials.Decision);
        end
        clear subject_trials
    end
    all_matrix(:,:,n) = sampled_matrix;
   clear sampled_matrix 
end

crt_locmeans(:,:,e) = mean(all_matrix,3);

clear h;
clear all_matrix
end
clear e;
% 
% crt_locmedian = zeros(200,3);
% 
% for e = 1:3
% 
% % create histogram of modified patch Luminance
% h = histogram(data_allab.Contrast(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);
% 
% % bin up the null patch responses, based on the modified patch histogram
% Null_bin = cell(1,h.NumBins);
% 
% for i = 1:h.NumBins
% Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.Contrast >= h.BinEdges(i)& data_allab.Contrast <= h.BinEdges(i+1)...
%     & data_allab.Eccentricity == location(e),:);
% end
% 
% % bootstrapping
% 
% 
% for n = 1:200
%     for c = 1:h.NumBins
%         if h.Values(c) ~= 0
%             
%             rand_index = randperm(height(Null_bin{:,c})); 
%             old_sq = Null_bin{:,c};
%             new_sq = old_sq(rand_index,:);
%             selected = new_sq(1:h.Values(:,c),:);
%             
%             if c == 1
%                 btp_data = selected.Decision;
%             else
%                 btp_data = [btp_data; selected.Decision];
%             end
%             
%         end
%         clear selected
%         clear old_sq
%         clear new_sq
%     end
%     crt_locmedian(n,e) = median(btp_data);
%     clear btp_data
%     
% end    
% clear h;
% end
% clear e;


%% RGB
for e = 1:3

% create histogram of modified patch Luminance
h = histogram(data_allab.RGB(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);

% bin up the null patch responses, based on the modified patch histogram
Null_bin = cell(1,h.NumBins);

for i = 1:h.NumBins
Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.RGB >= h.BinEdges(i)& data_allab.RGB <= h.BinEdges(i+1)...
    & data_allab.Eccentricity == location(e),:);
end

% bootstrapping

for n = 1:200
    sampled_trials = zeros(sum(h.NumBins .* h.Values),width(data_allab));
    for c = 1:h.NumBins
        if h.Values(c) ~= 0
            rand_index = randperm(height(Null_bin{:,c})); 
            old_sq = Null_bin{:,c};
            new_sq = old_sq(rand_index,:);
            selected = new_sq(1:h.Values(:,c),:);
            
            if c == 1
                sampled_trials = selected;
            else
                sampled_trials = [sampled_trials; selected];
            end
        end
    end
    
    for sub = 1:15
        subject_trials = sampled_trials(sampled_trials.Subject == sub,:);
        for i = 1:8
            sampled_matrix(sub,i) = size(subject_trials(subject_trials.Decision == decision_option(i),:),1)/...
            length(subject_trials.Decision);
        end
        clear subject_trials
    end
    all_matrix(:,:,n) = sampled_matrix;
   clear sampled_matrix 
end

RGB_locmeans(:,:,e) = mean(all_matrix,3);

clear h;
clear all_matrix
end
clear e;

% RGB_locmedians= zeros(200,3);
% 
% for e = 1:3
% 
% % create histogram of modified patch Luminance
% h = histogram(data_allab.RGB(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);
% 
% % bin up the null patch responses, based on the modified patch histogram
% Null_bin = cell(1,h.NumBins);
% 
% for i = 1:h.NumBins
% Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.RGB >= h.BinEdges(i)& data_allab.RGB <= h.BinEdges(i+1)...
%     & data_allab.Eccentricity == location(e),:);
% end
% 
% % bootstrapping
% 
% btp_trials = cell(200,20);
% 
% 
% for n = 1:200
%     for c = 1:h.NumBins
%         if h.Values(c) ~= 0
%             rand_index = randperm(height(Null_bin{:,c})); 
%             old_sq = Null_bin{:,c};
%             new_sq = old_sq(rand_index,:);
%             selected = new_sq(1:h.Values(:,c),:);
%             if c == 1
%                 btp_data = selected.Decision;
%             else
%                 btp_data = [btp_data;selected.Decision];
%             end
%         end
%         clear selected
%         clear old_sq
%         clear new_sq
%     end
%     RGB_locmedians(n,e) = median(btp_data);
%     clear btp_data
% end    
% clear h;
% end
% clear e;


%% check whether responses with three different matched perceptual properties are the same

RGB_submeans = median(mean(RGB_locmeans,3));
RGB_se = within_se(median(RGB_locmeans,3),15,8);

crt_submeans = mean(mean(crt_locmeans,3));
crt_se = within_se(mean(crt_locmeans,3),15,8);

lum_submeans = mean(mean(lum_locmedians,3));
lum_se = within_se(mean(lum_locmedians,3),15,8);

subplot(1,4,1), errorbar(decision_option, RGB_submeans,RGB_se,'.-','MarkerSize',2,'LineWidth',1);
hold on
errorbar(decision_option,crt_submeans,crt_se,'.-','MarkerSize',2,'LineWidth',1);
errorbar(decision_option,lum_submeans,lum_se,'.-','MarkerSize',2,'LineWidth',1);
hold off
legend({'RGB','Contrast','Luminance'});
xticks(decision_option);
title('Across eccentricities');
set(gca,'FontName','Arial','FontSize',12,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'});
ylim([0 0.7]);
xlim([-4 4]);
ylabel('Percentage of responses');

location = [0 6.5 9.2];

for loc = 1:3
RGB_submeans_loc = mean(RGB_locmeans(:,:,loc));
RGB_se_loc = within_se(RGB_locmeans(:,:,loc),15,8);

crt_submeans_loc = mean(crt_locmeans(:,:,loc));
crt_se_loc = within_se(crt_locmeans(:,:,loc),15,8);

lum_submeans_loc = mean(lum_locmedians(:,:,loc));
lum_se_loc = within_se(lum_locmedians(:,:,loc),15,8);

subplot(1,4,loc+1), errorbar(decision_option, RGB_submeans,RGB_se_loc,'.-','MarkerSize',2,'LineWidth',1);
hold on
errorbar(decision_option,crt_submeans,crt_se,'.-','MarkerSize',2,'LineWidth',1);
errorbar(decision_option,lum_submeans,lum_se,'.-','MarkerSize',2,'LineWidth',1);
hold off
xticks(decision_option);
set(gca,'FontName','Arial','FontSize',12,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'});
ylim([0 0.7]);
xlim([-4 4]);

if loc ~= 1
    yticks([]);
end

title(['Eccentricity = ', num2str(location(loc))]);
end

%% collapse across three properties

all_property_means = zeros(15,8,3);

for e = 1:3
current_matrix = zeros(15,8,3);
current_matrix(:,:,1) = RGB_locmeans(:,:,e);
current_matrix(:,:,2) =  crt_locmeans(:,:,e);
current_matrix(:,:,3) = lum_locmeans(:,:,e);
all_property_means(:,:,e) = mean(current_matrix,3);
end


across_ecc = mean(all_property_means,3);
across_se = std(across_ecc)/sqrt(15);
across_mean = mean(across_ecc);

loc1_se = std(all_property_means(:,:,1))/sqrt(15);
loc1_mean = mean(all_property_means(:,:,1));

loc2_se = std(all_property_means(:,:,2))/sqrt(15);
loc2_mean = mean(all_property_means(:,:,2));

loc3_se = std(all_property_means(:,:,3))/(15);
loc3_mean = mean(all_property_means(:,:,3));

clf

%% descriptives of real data, hypo 1
Results = importdata('Pooled Results 2.mat');
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
location = [0 6.5 9.2];

%% hypothesis 1
%proportion of responses for absent patches, on each Decision x Confidence
%category for each of 15 subjects

pcgN_matrix = zeros(15,9);
for sub = 1:15
for i = 1:9
    pcgN_matrix(sub,i) = size(Results(Find_N & Results(:,9)== i-5 & Results(:,1)==sub,:),1)/length(Results(Find_N & Results(:,1)==sub,:));
end
end

pcgN_matrix = [pcgN_matrix(:,1:4) pcgN_matrix(:,6:9)];

% mean proportion of responses across 15 participants on each response category
m_N = mean(pcgN_matrix);

% standard erors
se_N = std(pcgN_matrix)/sqrt(15);


%do the same for present patches
pcgAP_matrix = zeros(15,9);

for sub = 1:15
for i = 1:9
    pcgAP_matrix(sub,i) = size(Results((Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
        length(Results((Find_CAP|Find_IAP) & Results(:,1)== sub,:));
end
end
pcgAP_matrix = [pcgAP_matrix(:,1:4) pcgAP_matrix(:,6:9)];
m_AP = mean(pcgAP_matrix);
se_AP = std(pcgAP_matrix)/sqrt(15);


%% addtional LME analysis
sub_num = reshape([1:15],15,1);

for response_cat = 1:8
    if response_cat == 1
       data_m = [pcgAP_matrix(:,response_cat) sub_num zeros(15,1)+decision_option(response_cat)]
    else
       current_mat = [pcgAP_matrix(:,response_cat) sub_num zeros(15,1)+decision_option(response_cat)];
       data_m = [data_m; current_mat];
    end
    clear current_mat
end
data_q1 = table(data_m(:,1),data_m(:,2),data_m(:,3),'VariableNames',{'Percentage','Subject','Category'});

lme1 = fitlme(data_q1,'Percentage ~ Category^2 + (1|Subject)+ (Category|Subject)');
lme2 = fitlme(data_q1,'Percentage ~ Category + (1|Subject) + (Category|Subject)');
q1_effect = compare(lme2,lme1);


%% plot line graph
colours = cbrewer('qual', 'Set1', 8); 
out = figure;
subplot(3,4,1),errorbar([1:8],m_AP,se_AP,'.-','MarkerSize',12,'Color',colours(5,:),'MarkerEdgeColor',colours(5,:),'MarkerFaceColor',colours(5,:),'LineWidth',1);

hold on
errorbar(1:8,m_N,se_N,'.-','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta','LineWidth',1);
hold off

set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 0.7]);
xlim([0.5 8.5]);
title('Across eccentricities','FontName','Arial');

%% hypothesis 1 on eccentricity 


for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N(sub,i) = size(Results(Results(:,13)==location(loc)& Find_N & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            length(Results(Results(:,end)==location(loc)&Find_N & Results(:,1)==sub,:));
        pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)&(Find_IAP|Find_CAP) & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            length(Results(Results(:,13)==location(loc)&(Find_CAP|Find_IAP) & Results(:,1)== sub,:));
        end
    end

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];

if loc == 1
    pcg_N1 = pcg_N;
elseif loc == 2
    pcg_N2 = pcg_N;
else
    pcg_N3 = pcg_N;
end

pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];
 
percentage_N = mean(pcg_N,1);
percentage_AP = mean(pcg_AP,1);

se_AP = std(pcg_AP)/sqrt(15);
se_N =  std(pcg_N)/sqrt(15); 


m_All = [reshape(percentage_AP,[8,1]) reshape(percentage_N,[8,1])];

% addtional LME analysis
sub_num = reshape([1:15],15,1);

for response_cat = 1:8
    if response_cat == 1
       data_m = [pcg_AP(:,response_cat) sub_num zeros(15,1)+decision_option(response_cat)];
    else
       current_mat = [pcg_AP(:,response_cat) sub_num zeros(15,1)+decision_option(response_cat)];
       data_m = [data_m; current_mat];
    end
    clear current_mat
end
data_q1 = table(data_m(:,1),data_m(:,2),data_m(:,3),'VariableNames',{'Percentage','Subject','Category'});

lme1 = fitlme(data_q1,'Percentage ~ Category^2 + (1|Subject)+ (Category|Subject)');
lme2 = fitlme(data_q1,'Percentage ~ Category + (1|Subject) + (Category|Subject)');
disp(location(loc))
q1_effect = compare(lme2,lme1)

subplot(3,4,loc+1), errorbar([1:8],percentage_N,se_N,'.-', 'MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','Color','magenta',...
    'LineWidth',1);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
ylim([0 0.7]);
xlim([0.5 8.5]);
% xticks([]);
hold on
errorbar(1:8,percentage_AP,se_AP, '.-','MarkerSize',12,'MarkerEdgeColor', colours(5,:),'MarkerFaceColor',colours(5,:),'Color',colours(5,:),'LineWidth',1);
hold off

% if loc ~= 1
% yticks([]);
% end


    if loc ==1
    title('Eccentricity = 0 dva','FontWeight','normal');
    ylabel(' ');
    elseif loc == 2
    title('Eccentricity = 6.5 dva','FontWeight','normal');
    else
    title('Eccentricity = 9.2 dva','FontWeight','normal');
    end
    legend('off');
 
clear pcg_AP
clear pcg_N
clear percentage_N
clear percentage_AP
clear std_hypo3
clear m_hypo3
clear percentage_All
clear lme1
clear lme2 
clear q1_effect
clear data_q1

end

%% additional analysis on hypo1, absent patches
[h1,p1] = kstest2(sum(pcg_N1(:,1:4),2),sum(pcg_N2(:,1:4),2))
[h2,p2] = kstest2(sum(pcg_N2(:,1:4),2),sum(pcg_N3(:,1:4),2))
[h3,p3] = kstest2(sum(pcg_N1(:,1:4),2),sum(pcg_N3(:,1:4),2))

%% hypo 2

for i = 1:9
   
    for sub = 1:15
        percentage_CP(sub,i) = size(Results(Results(:,1)== sub & Find_Congruent_CP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Congruent_CP,:));
        percentage_CI(sub,i) = size(Results(Results(:,1)== sub &Find_Congruent_IP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Congruent_IP,:));
        percentage_IP(sub,i) =  size(Results(Results(:,1)== sub &Find_Incongruent_IP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Incongruent_IP,:));
        percentage_II(sub,i) =  size(Results(Results(:,1)== sub & Find_Incongruent_CP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== sub & Find_Incongruent_CP,:));
    end
end


percentage_CP = [percentage_CP(:,1:4) percentage_CP(:,6:9)];
percentage_CI = [percentage_CI(:,1:4) percentage_CI(:,6:9)];
percentage_IP = [percentage_IP(:,1:4) percentage_IP(:,6:9)];
percentage_II = [percentage_II(:,1:4) percentage_II(:,6:9)];

se_CP = std(percentage_CP)/sqrt(15);
se_CI = std(percentage_CI)/sqrt(15);
se_IP = std(percentage_IP)/sqrt(15);
se_II = std(percentage_II)/sqrt(15);


percentage_CP = mean(percentage_CP,1);
percentage_CI = mean(percentage_CI,1);
percentage_IP = mean(percentage_IP,1);
percentage_II = mean(percentage_II,1);

subplot(3,4,5),errorbar(1:8,percentage_CP,se_CP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
hold on
errorbar(1:8,percentage_CI,se_CI,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
    'LineWidth',1);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
errorbar(1:8, across_mean,across_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
ylim([0 0.7]),xlim([0.5 8.5]);
hold off

subplot(3,4,9),errorbar(1:8,percentage_IP,se_IP,'.--','MarkerSize',12,'Color',colours(2,:),'MarkerEdgeColor',colours(2,:),'LineWidth',1);
hold on
errorbar(1:8,percentage_II,se_II,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerEdgeColor',colours(1,:),'LineWidth',1.2);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');
errorbar(1:8, across_mean,across_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1.2);
xlim([0.5 8.5]),ylim([0 0.7]);
hold off

%% some 2-way anova
for cat = 1:8
cong_column = [percentage_CP(:,cat); percentage_CI(:,cat)];
incong_column = [percentage_IP(:,cat);percentage_II(:,cat)];
CI_mat = [cong_column incong_column];

p(cat,:) = anova2(CI_mat,15);

clear CI_mat
clear cong_column
clear incong_column
end

test_results = table(p(:,1),p(:,2),p(:,3),'VariableNames',{'modification','congruency','interaction'})

%% hypo2, across eccentricity
% across_ecc = mean(all_property_means,3);
% across_se = within_se(across_ecc,15,8);
% across_mean = mean(across_ecc);
% 
% loc1_se = within_se(all_property_means(:,:,1),15,8);
% loc1_mean = mean(all_property_means(:,:,1));
% 
% loc2_se = within_se(all_property_means(:,:,2),15,8);
% loc2_mean = mean(all_property_means(:,:,2));
% 
% loc3_se = within_se(all_property_means(:,:,3),15,8);
% loc3_mean = mean(all_property_means(:,:,3));


for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N(sub,i) = size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,1)==sub,:),1);
        pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,1)== sub,:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = nanmean(pcg_N,1);
percentage_AP = nanmean(pcg_AP,1);
se_AP = std(pcg_AP)/sqrt(15);
se_N = std(pcg_N)/sqrt(15);
    
subplot(3,4,loc+5),errorbar(1:8,percentage_AP,se_AP,'.-','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
hold on
errorbar(1:8,percentage_N,se_N,'.-','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
    'LineWidth',1);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');

if loc == 1
    errorbar(1:8, loc1_mean,loc1_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
elseif loc == 2
    errorbar(1:8,loc2_mean,loc2_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
elseif loc == 3
    errorbar(1:8,loc3_mean,loc3_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
end

ylim([0 0.7]),xlim([0.5 8.5]);
hold off
% 
% if loc ~= 1
%     yticks([]);
% end



    
clear pcg_AP
clear pcg_N
clear percentage_N
clear percentage_AP
clear std_hypo3
clear m_hypo3
clear percentage_All

end



for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N(sub,i) = size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,1)==sub,:),1);
        pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,1)== sub,:),1);
        end
    end
    

pcg_N = [pcg_N(:,1:4) pcg_N(:,6:9)];
pcg_AP = [pcg_AP(:,1:4) pcg_AP(:,6:9)];

percentage_N = nanmean(pcg_N,1);
percentage_AP = nanmean(pcg_AP,1);
se_AP = std(pcg_AP)/sqrt(15);
se_N = std(pcg_N)/sqrt(15);

subplot(3,4,loc+9),errorbar(1:8,percentage_AP,se_AP,'.--','MarkerSize',12,'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',...
    colours(2,:),'LineWidth',1);
hold on
errorbar(1:8,percentage_N,se_N,'.--','MarkerSize',12,'Color',colours(1,:),'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),...
    'LineWidth',1);
set(gca,'XTick',[1:1:8],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'FontSize', 12,'FontName','Arial','Box','off');

if loc == 1
    errorbar(1:8, loc1_mean,loc1_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
elseif loc == 2
    errorbar(1:8,loc2_mean,loc2_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
elseif loc == 3
    errorbar(1:8,loc3_mean,loc3_se,'.--','MarkerSize',12,'Color','magenta','MarkerEdgeColor','magenta','MarkerFaceColor','magenta',...
    'LineWidth',1);
end

ylim([0 0.7]),xlim([0.5 8.5]);
hold off
    
% if loc ~= 1
%     yticks([]);
% end

 

% clear pcg_AP
% clear pcg_N
% clear percentage_N
% clear percentage_AP
% clear std_hypo3
% clear m_hypo3
% clear percentage_All

end


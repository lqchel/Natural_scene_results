% additional anova and t-tests

%Results = [pilot.b2_g1;pilot.b2_g2;pilot.b2_g3]; %% uncomment for experiment 2
%Results = importdata('Exp2_data.mat');
Results = importdata('Pooled Results 2.mat');
%Results = Results(Results(:,11)~=0,:);
Results(:,9) = Results(:,8).*Results(:,9);

%Find trials that presented N patches -- signal absent for hypo 1
%Find trials that presented N patches -- signal absent for hypo 1
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

subject_id = unique(Results(:,1));
num_sub = length(subject_id);
location = [0 1 2];

 
%% LME model for whether percentage of no responses is associated with eccentricity
a = 1;

for loc = 1:3
        for sub = 1:num_sub
            img_num = unique(Results(Results(:,1)==subject_id(sub) & Results(:,13) == location(loc),11));
            for img = 1:length(img_num)
                if sum(Results(:,13)==location(loc)& Find_N & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub)) == 0
                    continue
                end
                loc_matrix(a,1) = subject_id(sub);
                loc_matrix(a,2) = img_num(img);
                loc_matrix(a,3) = unique(Results(Results(:,13)==location(loc)& Find_N & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub),14));
                loc_matrix(a,4) = unique(Results(Results(:,13)==location(loc)& Find_N & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub),13));
                loc_matrix(a,5) = size(Results(Results(:,13)==location(loc)& Find_N & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub) & Results(:,8) == -1,:),1)/...
                    size(Results(Results(:,13)==location(loc)&Find_N & Results(:,1)==subject_id(sub) & Results(:,11)== img_num(img),:),1);
                a = a+1;
            end
            clear img_num
        end
end

lme1_data = table(loc_matrix(:,1),loc_matrix(:,2),loc_matrix(:,3),loc_matrix(:,5),'VariableNames',{'subject','image','ecc','percentage'});
lm1 = fitlme(lme1_data,'percentage ~ ecc + (ecc | subject)+ (ecc|image)+(1|subject) + (1|image)');
lm2 = fitlme(lme1_data,'percentage~ 1+ (ecc | subject)+ (ecc|image)+(1|subject) + (1|image)');
compare(lm2,lm1)


%% LME model for whether percentage of yes responses is associated with eccentricity

a = 1;
Find_AP = Find_IAP | Find_CAP;

for loc = 1:3
        for sub = 1:num_sub
            img_num = unique(Results(Results(:,1)==subject_id(sub) & Results(:,13) == location(loc),11));
            for img = 1:length(img_num)
                if sum(Results(:,13)==location(loc)& Find_AP & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub)) == 0
                    continue
                end
                loc_matrix(a,1) = subject_id(sub);
                loc_matrix(a,2) = img_num(img);
                loc_matrix(a,3) = unique(Results(Results(:,13)==location(loc)& Find_AP & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub),14));
                loc_matrix(a,4) = unique(Results(Results(:,13)==location(loc)& Find_AP & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub),13));
                loc_matrix(a,5) = size(Results(Results(:,13)==location(loc)& Find_AP & Results(:,11)== img_num(img) & Results(:,1)== subject_id(sub) & Results(:,8) == 1,:),1)/...
                    size(Results(Results(:,13)==location(loc)&Find_AP & Results(:,1)==subject_id(sub) & Results(:,11)== img_num(img),:),1);
                a = a+1;
            end
            clear img_num
        end
end

lme2_data = table(loc_matrix(:,1),loc_matrix(:,2),loc_matrix(:,3),loc_matrix(:,5),'VariableNames',{'subject','image','ecc','percentage'});
lm3 = fitlme(lme2_data,'percentage ~ ecc +(ecc | subject)+ (ecc|image)+ (1|subject) + (1|image)');
lm4 = fitlme(lme2_data,'percentage~ 1+ (ecc | subject)+ (ecc|image)+(1|subject) + (1|image)');
compare(lm4,lm3)


% for loc = 1:3
%     for i = 1:9
%         for sub = 1:num_sub
%         pcg_AP(sub,i) = size(Results(Results(:,13)==location(loc)& (Find_IAP | Find_CAP)& Results(:,9)== i-5 & Results(:,1)== subject_id(sub),:),1)/...
%             size(Results(Results(:,13)==location(loc)&(Find_IAP | Find_CAP) & Results(:,1)==subject_id(sub),:),1);
%         end
%     end
% temp1 = reshape([pcg_AP(:,6:9)],[4.*num_sub,1]);
% temp2 = [subject_id; subject_id; subject_id; subject_id];
% current_loc_matrix = [temp2 location(loc)+zeros(4.*num_sub,1) temp1];
% if loc == 1
%     loc_matrix = current_loc_matrix;
% else
%     loc_matrix = [loc_matrix; current_loc_matrix];
% end
% 
% end
% 
% % no responses across eccentricities
% data = table(loc_matrix(:,1),loc_matrix(:,2),loc_matrix(:,3),'VariableNames',{'subject','ecc','percentage'});
% lm1 = fitlme(data,'percentage ~ ecc + (1|subject)');
% lm2 = fitlme(data,'percentage~ 1+ (1|subject)');
% compare(lm2,lm1)

%% comparing between original and modified, for each response category

for i = 1:9
    for sub = 1:num_sub
        percentage_CP(sub,i) = size(Results(Results(:,1)== subject_id(sub) & Find_Congruent_CP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== subject_id(sub) & Find_Congruent_CP,:));
        percentage_CI(sub,i) = size(Results(Results(:,1)== subject_id(sub) &Find_Congruent_IP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== subject_id(sub) & Find_Congruent_IP,:));
        percentage_IP(sub,i) =  size(Results(Results(:,1)== subject_id(sub) &Find_Incongruent_IP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== subject_id(sub) & Find_Incongruent_IP,:));
        percentage_II(sub,i) =  size(Results(Results(:,1)== subject_id(sub) & Find_Incongruent_CP & Results(:,9)== i-5,:),1)/...
            length(Results(Results(:,1)== subject_id(sub) & Find_Incongruent_CP,:));
    end
end


percentage_CP = [percentage_CP(:,1:4) percentage_CP(:,6:9)];
percentage_CI = [percentage_CI(:,1:4) percentage_CI(:,6:9)];
percentage_IP = [percentage_IP(:,1:4) percentage_IP(:,6:9)];
percentage_II = [percentage_II(:,1:4) percentage_II(:,6:9)];

% t tests on each response category

for i = 1:8
    [h1,p1] = ttest(percentage_CP(:,i),percentage_CI(:,i),'alpha',0.003);
    cong_hypo2_t(i) = p1;
    cong_hypo2_s(i) = h1;
    
    [h2,p2] = ttest(percentage_IP(:,i),percentage_II(:,i),'alpha',0.003);
    incong_hypo2_t(i) = p2;
    incong_hypo2_s(i) = h2;
    
    clear p1
    clear p2
    clear h1
    clear h2
end
disp(cong_hypo2_s)
disp(cong_hypo2_t)
disp(incong_hypo2_s)
disp(incong_hypo2_t)

% on each eccentricity

location = [0 1 2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
        pcg_N_cong(sub,i) = size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,1)==subject_id(sub),:),1);
        pcg_AP_cong(sub,i) = size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,1)== subject_id(sub),:),1);
        pcg_N_incong(sub,i) = size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,1)==subject_id(sub),:),1);
        pcg_AP_incong(sub,i) = size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,1)== subject_id(sub),:),1);
        end
    end
    

    pcg_N_cong = [pcg_N_cong(:,1:4) pcg_N_cong(:,6:9)];
    pcg_AP_cong = [pcg_AP_cong(:,1:4) pcg_AP_cong(:,6:9)];
    pcg_N_incong = [pcg_N_incong(:,1:4) pcg_N_incong(:,6:9)];
    pcg_AP_incong = [pcg_AP_incong(:,1:4) pcg_AP_incong(:,6:9)];

    for i = 1:8
        [h1,p1] = ttest(pcg_AP_cong(:,i),pcg_N_cong(:,i),'alpha',0.001);
        cong_ecc_t(loc,i) = p1;
        cong_ecc_s(loc,i) = h1;

        [h2,p2] = ttest(pcg_AP_incong(:,i),pcg_N_incong(:,i),'alpha',0.001);
        incong_ecc_t(loc,i) = p2;
        incong_ecc_s(loc,i) = h2;

        clear p1
        clear p2
        clear h1
        clear h2
    end
    
    clear pcg_N_cong
    clear pcg_AP_cong
    clear pcg_N_incong
    clear pcg_AP_incong

end
disp(cong_ecc_s)
disp(cong_ecc_t)
disp(incong_ecc_s)
disp(incong_ecc_t)

%% comparing specific DxC categories, among different eccentricity levels

location = [0 1 2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:num_sub
        pcg_N_cong(sub,i,loc) = size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Congruent_IP & Results(:,1)==subject_id(sub),:),1);
        pcg_AP_cong(sub,i,loc) = size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Congruent_CP & Results(:,1)== subject_id(sub),:),1);
        pcg_N_incong(sub,i,loc) = size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)& Find_Incongruent_CP & Results(:,1)==subject_id(sub),:),1);
        pcg_AP_incong(sub,i,loc) = size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== subject_id(sub),:),1)/...
            size(Results(Results(:,13)==location(loc)&Find_Incongruent_IP & Results(:,1)== subject_id(sub),:),1);
        end
    end
    
end

pcg_N_cong = pcg_N_cong(:,1,:); % take DxC = -4 only for congruent images
pcg_AP_cong = pcg_AP_cong(:,1,:);
pcg_N_incong = pcg_N_incong(:,9,:); % take DxC = 4 only for incongruent images
pcg_AP_incong =  pcg_AP_incong(:,9,:);

% compare DxC = -4 for eccentricity 0 vs. 1, 0 vs. 2, 1 vs. 2, for
% congruent images
[h1,p1] = kstest2(pcg_N_cong(~isnan(pcg_N_cong(:,:,1)),:,1),pcg_N_cong(~isnan(pcg_N_cong(:,:,2)),:,2),'alpha',0.05/3) % 0 vs. 1
[h2,p2] = kstest2(pcg_N_cong(~isnan(pcg_N_cong(:,:,1)),:,1),pcg_N_cong(~isnan(pcg_N_cong(:,:,3)),:,3),'alpha',0.05/3) % 0 vs. 2
[h3,p3] = kstest2(pcg_N_cong(~isnan(pcg_N_cong(:,:,2)),:,2),pcg_N_cong(~isnan(pcg_N_cong(:,:,3)),:,3),'alpha',0.05/3) % 1 vs. 2

% compare DxC = 4 for eccentricity 0 vs. 1, 0 vs. 2, 1 vs. 2, for
% incongruent images


%%

for sub = 1:num_sub
exp_hit(sub) = size(Results((Find_CAP|Find_IAP) & Results(:,8) == 1 & Results(:,1) == sub,:),1)/ size(Results((Find_CAP | Find_IAP) & Results(:,1) == sub,:),1);
exp_miss(sub) = size(Results(Find_N & Results(:,8) == 1 & Results(:,1) == sub,:),1) / size(Results(Find_N & Results(:,1) == sub,:),1);
end
x = categorical({'Present','Absent'});
x = reordercats(x,{'Present','Absent'});
colours = cbrewer('qual','Pastel2',4);
bar(x,[mean(exp_hit) mean(exp_miss)],0.5,'FaceColor',colours(3,:));
hold on
errorbar([mean(exp_hit) mean(exp_miss)],[std(exp_hit)/sqrt(num_sub) std(exp_miss)/sqrt(num_sub)],'k.', 'LineWidth', 1);
hold off

ylim([0 1]),ylabel('Percentage of present judgments')
set(gca,'FontName','Arial','FontSize',12)


%% each temporal position
exp_patch_number = [21 6];
b = 1;
for exp = 1:2
    if exp == 1
        Results = importdata('Pooled Results 2.mat');
    else
        Results = importdata('Exp2_data.mat');
        Results = Results(Results(:,11)~=0,:);
    end
    Results(:,9) = Results(:,8).*Results(:,9);
    %Find trials that presented N patches -- signal absent for hypo 1
    %Find trials that presented N patches -- signal absent for hypo 1
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

    location = [0 1 2];
    
    for position = 1:exp_patch_number(exp)
        current_hit = size(Results((Find_CAP|Find_IAP) & Results(:,8) == 1 & Results(:,3) == position,:),1) / size(Results((Find_CAP | Find_IAP) & Results(:,3) == position,:),1);
        all_hits(b,:) = [exp position current_hit];
        current_miss = size(Results(Find_N & Results(:,8) == 1 & Results(:,3) == position,:),1) / size(Results(Find_N & Results(:,3) == position,:),1);
        all_misses(b,:) = [exp position current_miss];
        current_mod_FA = size(Results((Find_Congruent_IP|Find_Incongruent_CP) & Results(:,8) == 1 & Results(:,3) == position,:),1) /...
            size(Results((Find_Congruent_IP|Find_Incongruent_CP) & Results(:,3) == position,:),1);
        all_mod_FA(b,:) = [exp position current_mod_FA];
        b = b+1;
    end
clear Results
end

colours_1 = cbrewer('qual','Set2',8);
colours_2 = cbrewer('qual','Set1',4);
%subplot(1,2,1);
plot(all_hits(all_hits(:,1) == 1,2),all_hits(all_hits(:,1) == 1,3),'.-','Color',colours_1(3,:),'LineWidth',1.5);
hold on
plot(all_misses(all_misses(:,1) == 1,2),all_misses(all_misses(:,1) == 1,3),'.-','Color',colours_1(2,:),'LineWidth',1.5);
plot(all_mod_FA(all_mod_FA(:,1) == 1,2),all_mod_FA(all_mod_FA(:,1) == 1,3),'.-','Color',colours_1(4,:),'LineWidth',1.5);
plot(all_hits(all_hits(:,1) == 2,2),all_hits(all_hits(:,1) == 2,3),'d-','Color',colours_2(2,:),'LineWidth',1.5);
plot(all_misses(all_misses(:,1) == 2,2),all_misses(all_misses(:,1) == 2,3),'d-','Color',colours_2(1,:),'LineWidth',1.5);
plot(all_mod_FA(all_mod_FA(:,1) == 2,2),all_mod_FA(all_mod_FA(:,1) == 2,3),'d-','Color','magenta','LineWidth',1.5);
hold off
ylabel('Percentage of present judgements');
set(gca,'FontSize',12);
xlabel('Patch order within trial');
xlim([1 21]);

clear all_hits
clear all_misses











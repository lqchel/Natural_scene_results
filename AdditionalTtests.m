% additional anova and t-tests

%Results = [pilot.b2_g1;pilot.b2_g2;pilot.b2_g3]; %% uncomment for experiment 2
%Results = importdata('Exp2_data.mat');
Results = get_data2(1);


%%%% uncomment for checking Exp 1 data, first 6 patches of each trial
Results = Results(Results(:,3)<= 6,:);

%%%% uncomment for checking Exp 1 data, first 24 trials of each
%%%% participant
% Results = Results(Results(:,2)<= 27,:);


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

 
%% comparing between present + original and null patches, for each response category
Find_A = Find_CAP | Find_IAP;
for loc = 1:4
    for i = 1:9
        for sub = 1:num_sub
            if loc == 1 % caculate percentage of responses to present and null, across all eccentricity levels
            percentage_A(sub,i) = size(Results(Results(:,1)== subject_id(sub) & Find_A & Results(:,9)== i-5,:),1)/...
                length(Results(Results(:,1)== subject_id(sub) & Find_A,:));
            percentage_N(sub,i) = size(Results(Results(:,1)== subject_id(sub) &Find_N & Results(:,9)== i-5,:),1)/...
                length(Results(Results(:,1)== subject_id(sub) & Find_N,:));
            else % caculate percentage of responses to present and null, on each eccentricity level
            percentage_A(sub,i) = size(Results(Results(:,1)== subject_id(sub) & Find_A & Results(:,9)== i-5 & Results(:,13) == location(loc-1),:) ,1)/...
                length(Results(Results(:,1)== subject_id(sub) & Find_A & Results(:,13) == location(loc -1),:));
            percentage_N(sub,i) = size(Results(Results(:,1)== subject_id(sub) &Find_N & Results(:,9)== i-5 & Results(:,13) == location(loc-1),:),1)/...
                length(Results(Results(:,1)== subject_id(sub) & Find_N & Results(:,13) == location(loc -1),:));
            end
        end
    end


    percentage_A = [percentage_A(:,1:4) percentage_A(:,6:9)];
    percentage_N = [percentage_N(:,1:4) percentage_N(:,6:9)];


    % t tests on each response category

    for i = 1:8
        [h1,p1] = ttest(percentage_A(:,i),percentage_N(:,i),'alpha',0.003);
        hypo1_p(loc,i) = p1;
        hypo1_s(loc,i) = h1;
        clear p1
        clear p2
        clear h1
        clear h2
    end
    clear percentage_A
    clear percentage_N
    

end

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











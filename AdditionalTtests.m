% additional anova and t-tests

% load data and indexes
clear all
Results = importdata('Pooled Results.mat');
Results(:,9) = Results(:,8).*Results(:,9);

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


%% hypohtesis 2, congruent

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

% t tests on each response category

for i = 1:8
    [h1,p1] = ttest(percentage_CP(:,i),percentage_CI(:,i),'alpha',0.003);
    cong_hypo2_t(i) = p1;
    
    [h2,p2] = ttest(percentage_IP(:,i),percentage_II(:,i),'alpha',0.003);
    incong_hypo2_t(i) = p2;
    
    clear p1
    clear p1
end

%% on each eccentricity

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];
location = [0 6.5 9.2];

for loc = 1:3
    
    for i = 1:9
        for sub = 1:15
        pcg_N_cong(sub,i) = size(Results(Results(:,end)==location(loc)& Find_Congruent_IP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)& Find_Congruent_IP & Results(:,1)==sub,:),1);
        pcg_AP_cong(sub,i) = size(Results(Results(:,end)==location(loc)&Find_Congruent_CP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)&Find_Congruent_CP & Results(:,1)== sub,:),1);
        pcg_N_incong(sub,i) = size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,9)== i-5 & Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)& Find_Incongruent_CP & Results(:,1)==sub,:),1);
        pcg_AP_incong(sub,i) = size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,9)== i-5& Results(:,1)== sub,:),1)/...
            size(Results(Results(:,end)==location(loc)&Find_Incongruent_IP & Results(:,1)== sub,:),1);
        end
    end
    

    pcg_N_cong = [pcg_N_cong(:,1:4) pcg_N_cong(:,6:9)];
    pcg_AP_cong = [pcg_AP_cong(:,1:4) pcg_AP_cong(:,6:9)];
    pcg_N_incong = [pcg_N_incong(:,1:4) pcg_N_incong(:,6:9)];
    pcg_AP_incong = [pcg_AP_incong(:,1:4) pcg_AP_incong(:,6:9)];

    for i = 1:8
        [h1,p1] = ttest(pcg_AP_cong(:,i),pcg_N_cong(:,i),'alpha',0.001);
        cong_ecc_t(loc,i) = p1;

        [h2,p2] = ttest(pcg_AP_incong(:,i),pcg_N_incong(:,i),'alpha',0.001);
        incong_ecc_t(loc,i) = p2;

        clear p1
        clear p2
    end
    
    clear pcg_N_cong
    clear pcg_AP_cong
    clear pcg_N_incong
    clear pcg_AP_incong

end






Results = importdata('Pooled Results 2.mat'); %uncomment for lab testing
%results

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

%% hypothesis 1

colours = cbrewer('qual', 'Set1', 8); 
grandmatrix = zeros(2,2);
matrix1 = zeros(num_sub,2);

for sub = 1:num_sub
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
% [h,p,ci]= ttest(matrix1(:,1),matrix1(:,2));

% calculate standard errors
grandmatrix(1,:) = mean(matrix1);

for condition = 1:2
grandmatrix(2,condition) = std(matrix1(:,condition))/sqrt(num_sub);
end
patch = categorical({'present','absent'});
patch = reordercats(patch,{'present','absent'});
out = figure;
subplot(1,2,1), bar(patch,grandmatrix(1,:),'BarWidth',0.7), ylabel('%Accurate judgments');
hold on
errorbar(patch,grandmatrix(1,:),grandmatrix(2,:),'k.','LineWidth',1);
hold off

%% hypothesis 1 with eccentricity


colours = cbrewer('qual', 'Set1', 8); 
grandmatrix1 = zeros(2,3,2);

% hit rate
matrix1 = zeros(num_sub,3,2);
location = [0 6.48 9.16];

for loc = 1:3
for sub = 1:num_sub
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
 
    matrix1(sub,loc,condition) = accuracy;
    
    clear Results_P
    clear condition_mean 
    clear accuracy
end
end

end

for condition = 1:2 
    
    grandmatrix1(1,:,condition) = mean(matrix1(:,:,condition));
    grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
    
end

figure(out);
subplot(2,2,1),errorbar(-5,grandmatrix(1,1),grandmatrix(2,1),'d','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8);
hold on
errorbar(-5,grandmatrix(1,2),grandmatrix(2,2),'o','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
    'k','Color','k','LineWidth',0.8);
hold off
ylabel('%Accurate judgments');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.50','9.20'},'FontSize',12);
ylim([0.4 1]);
legend('off');

hold on
errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
    'k','Color','k','LineWidth',0.8,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
legend({'Hit','CR'});

clear grandmatrix1
clear matrix1
clear R_indv
clear grandmatrix
clear Find_patch

%% hypothesis 1 confidence

% grandmatrix = zeros(2,2);
% matrix1 = zeros(num_sub,2);
% 
% for sub = 1:num_sub
%     R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
%     R_indv(:,9) = R_indv(:,9).* R_indv(:,8);
% for condition = 1:2
%     if condition == 1
%         Find_patch = R_indv(:,5)~=1 & R_indv(:,8) ==1;
%         Results_P = R_indv(Find_patch,:);
%         confidence = sum(Results_P(:,9))/size(Results_P,1);
%     else
%         Find_patch = R_indv(:,5)==1 & R_indv(:,8) == -1;
%         Results_P = R_indv(Find_patch,:);
%         confidence = sum(Results_P(:,9))/size(Results_P,1);
%     end
%     
%     matrix1(sub,condition) = confidence;
%     
%     clear Results_P
%     clear condition_mean 
%     clear confidence
%     clear Find_patch
% end
% end
% 
% % calculate standard errors
% grandmatrix(1,:) = mean(matrix1);
% grandmatrix(2,:) = within_se(matrix1,num_sub,2);
% 
% 
% 
% % across eccetricities
% grandmatrix1 = zeros(2,3,2);
% matrix1 = zeros(num_sub,3,2);
% location = [0 6.48 9.16];
% 
% for loc = 1:3
% for sub = 1:num_sub
%     R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N)& Results(:,end)==location(loc),:);
%     R_indv(:,9) = R_indv(:,9).* R_indv(:,8);
% for condition = 1:2
%    if condition == 1
%         Find_patch = R_indv(:,5)~=1 & R_indv(:,8) ==1;
%         Results_P = R_indv(Find_patch,:);
%         confidence = sum(Results_P(:,9))/size(Results_P,1);
%     else
%         Find_patch = R_indv(:,5)==1 & R_indv(:,8) == -1;
%         Results_P = R_indv(Find_patch,:);
%         confidence = sum(Results_P(:,9))/size(Results_P,1);
%     end
%      
%     matrix1(sub,loc,condition) = confidence;
%     
%     clear Results_P
%     clear condition_mean 
%     clear accuracy
% end
% end
% 
% end
% 
% for condition = 1:2 
%     
%     grandmatrix1(1,:,condition) = mean(matrix1(:,:,condition));
%     grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
%     
% end
% 
% figure(out);
% subplot(2,3,2),errorbar(-5,grandmatrix(1,1),grandmatrix(2,1),'d','MarkerSize',8,...
%     'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8);
% hold on
% errorbar(-5,grandmatrix(1,2),grandmatrix(2,2),'o','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
%     'k','Color','k','LineWidth',0.8);
% hold off
% ylabel('Confidence level');
% xlim([-7 11]),xticks([-5 0 6.48 9.16]);
% set(gca,'XTickLabel',{'All','0','6.50','9.20'},'FontSize',12);
% ylim([1 4]),yticks([1:1:4]);
% legend('off');
% 
% hold on
% errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',8,...
%     'MarkerFaceColor',colours(6,:),'MarkerEdgeColor','k','Color','k','LineWidth',0.8,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',8,'MarkerFaceColor',colours(5,:),'MarkerEdgeColor',...
%     'k','Color','k','LineWidth',0.8,'Capsize',10);
% hold off
% legend({'Hit','CR'});
% 
% 
% %confidence plot for regan
% 
% grandmatrix = zeros(2,2);
% matrix1 = zeros(num_sub,2);
% 
% for sub = 1:num_sub
%     R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
%     R_indv(:,9) = abs(R_indv(:,9));
% for condition = 1:2
%     if condition == 1
%         Find_patch = R_indv(:,5)~=1;
%         Results_P = R_indv(Find_patch,:);
%         confidence = sum(Results_P(:,9))/size(Results_P,1);
%     else
%         Find_patch = R_indv(:,5)==1;
%         Results_P = R_indv(Find_patch,:);
%         confidence = sum(Results_P(:,9))/size(Results_P,1);
%     end
%     
%     matrix1(sub,condition) = confidence;
%     
%     clear Results_P
%     clear condition_mean 
%     clear confidence
%     clear Find_patch
% end
% end
% 
% % calculate standard errors
% grandmatrix(1,:) = mean(matrix1);
% 
% 
% for condition = 1:2
% grandmatrix(2,condition) = std(matrix1(:,condition))/sqrt(num_sub);
% end
% 
% patch = categorical({'present','absent'});
% patch = reordercats(patch,{'present','absent'});
% subplot(1,2,2), bar(patch,abs(grandmatrix(1,:)),'BarWidth',0.7), ylabel('%Accurate judgments');
% hold on
% errorbar(patch,abs(grandmatrix(1,:)),grandmatrix(2,:),'k.','LineWidth',1);
% hold off


%% hypothesis 1 AUC

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];

Results_NC = Results(Find_N,:); % find all absent patches
Results_APC = Results(Find_IAP|Find_CAP,:); % find all present patches
matrix1 = zeros(num_sub,1); % create empty matrix for recording subject AUC

for sub = 1:num_sub 
    Results_NC = Results(Find_N,:); % find absent patches for the current subject
    Confidence_N = Results_NC(Results_NC(:,1)==sub,9); 
    Confidence_AP = Results_APC(Results_APC(:,1)==sub,9); % find present patch responses for current subject
    
    % calculate number of response on each response category (-4 to 4), for present
    % (AP) and absent (N) patches
    
    for i = -4:4
        Confidence_APCounts(i+5) = sum(Confidence_AP == -i); 
        Confidence_NCounts(i+5) = sum(Confidence_N == -i);
    end
    
    % calculate cumulative hit and FA on each response category
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

    % calculate AUC of the subject
    Cumulative_Hit = [0 Cumulative_Hit(1:4) Cumulative_Hit(6:9)];
    Cumulative_FA = [0 Cumulative_FA(1:4) Cumulative_FA(6:9)];
    AUC = round(AreaUnderROC([Cumulative_Hit; Cumulative_FA]'),2);
    
    % record AUC of subject
    matrix1(sub,:)= AUC;



end

se1 = std(AUC)/sqrt(num_sub);

% AUC across eccentricities

matrix3 = zeros(num_sub, 3);
for sub = 1:num_sub
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

se2 = within_se(matrix3,num_sub,3);

% plot graph
figure(out)
subplot(2,2,2),errorbar(-5,mean(matrix1),se1,'d','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
ylabel('AUC');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.50','9.20'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar([0 6.48 9.16],mean(matrix3),se2,'d-','MarkerSize',8,...
    'MarkerFaceColor',colours(6,:),'MarkerEdgeColor',colours(5,:),'Color',colours(5,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off




%% hypothesis 2
clear grandmatrix
clear matrix1

% hypothesis 2 eccentricity

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
eccentricity = zeros(length(Results),1)+ location1 + location2;
Results = [Results eccentricity];
    % congruent ---------------------------------------------------------------------------------------------------
    matrix1 = zeros(num_sub,3,2);
    location = [0 6.48 9.16];

    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
        end

        matrix1(sub,loc,condition) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end

    end
    disp(matrix1); pause(1);
    for condition = 1:2 

        grandmatrix1(1,:,condition) = mean(matrix1(:,:,condition));
        grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
        
    end
    


    %incongruent --------------------------------------------------------------------------------
     matrix2 = zeros(num_sub,3,2);

    for loc = 1:3
    for sub = 1:num_sub
        R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP|Find_Incongruent_IP)& Results(:,end)==location(loc),:);

    for condition = 1:2
        if condition == 1
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        else
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==-1)/size(Results_P,1);
        end

        matrix2(sub,loc,condition) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end

    end

    disp(matrix2);pause(1);
    for condition = 1:2 

        grandmatrix2(1,:,condition) = mean(matrix2(:,:,condition));
        grandmatrix2(2,:,condition)= within_se(matrix2(:,:,condition),num_sub,3);
        
    end

    % calculate hit and CR collapsed across eccentricity levels
    all_cong(1,:) = [mean(grandmatrix1(1,:,1)) mean(grandmatrix1(1,:,2))];
    all_cong(2,:) = within_se([mean(matrix1(:,:,1),2) mean(matrix1(:,:,2),2)],num_sub,2);
    all_incong(1,:) = [mean(grandmatrix2(1,:,1)) mean(grandmatrix2(1,:,2))];
    all_incong(2,:) = within_se([mean(matrix2(:,:,1),2) mean(matrix2(:,:,2),2)],num_sub,2);
    
    % plot figures
    % congruent
    
figure(out);

subplot(2,2,3),errorbar(-5,all_cong(1,1),all_cong(2,1),'d','MarkerSize',8,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
hold on
errorbar(-5, all_cong(1,2), all_cong(2,2),'o','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(1,:),'LineWidth',1);
errorbar(-5,all_incong(1,1),all_incong(2,1),'d','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
    colours(2,:),'Color',colours(2,:),'LineWidth',1);
errorbar(-5,all_incong(1,2),all_incong(2,2),'o','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
    colours(1,:),'Color',colours(1,:),'LineWidth',1);
hold off
ylabel('%Accurate judgments');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.50','9.20'},'FontSize',12);
ylim([0 1]);
legend('off');


hold on
errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',6,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',...
    colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix2(1,:,1),grandmatrix2(2,:,1),'d--','MarkerSize',7,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
errorbar([0 6.48 9.16],grandmatrix2(1,:,2),grandmatrix2(2,:,2),'o--','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
    colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
plot([-7 11],[0.5 0.5],'k--');
hold off
legend({'Congruent Hit','Congruent CR','Incongruent Hit','Incongruent CR'});


% % hypothesis 2 confidence
% 
% location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.48;
% location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.16;
% eccentricity = zeros(length(Results),1)+ location1 + location2;
% Results = [Results eccentricity];
% Results(:,9) = abs(Results(:,9));
% 
%     % congruent ---------------------------------------------------------------------------------------------------
%     matrix1 = zeros(num_sub,3,2);
%     location = [0 6.48 9.16];
% 
%     for loc = 1:3
%     for sub = 1:num_sub
%         R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP)& Results(:,end)==location(loc),:);
% 
%     for condition = 1:2
%         if condition == 1
%             Find_patch = R_indv(:,5)==2 & R_indv(:,8)==1;
%             Results_P = R_indv(Find_patch,:);
%             confidence = sum(Results_P(:,9))/size(Results_P,1);
%         else
%             Find_patch = R_indv(:,5)==3 & R_indv(:,8)==-1;
%             Results_P = R_indv(Find_patch,:);
%             confidence = sum(Results_P(:,9))/size(Results_P,1);
%         end
% 
%         matrix1(sub,loc,condition) = confidence;
% 
%         clear Results_P
%         clear condition_mean 
%         clear confidence
%     end
%     end
% 
%     end
% 
%     for condition = 1:2 
% 
%         grandmatrix1(1,:,condition) = nanmean(matrix1(:,:,condition));
%         grandmatrix1(2,:,condition)= within_se(matrix1(:,:,condition),num_sub,3);
%         
%     end
%     
% 
% 
%     %incongruent --------------------------------------------------------------------------------
%      matrix2 = zeros(num_sub,3,2);
% 
%     for loc = 1:3
%     for sub = 1:num_sub
%         R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP|Find_Incongruent_IP)& Results(:,end)==location(loc),:);
% 
%     for condition = 1:2
%         if condition == 1
%             Find_patch = R_indv(:,5)==3 & R_indv(:,8) == 1;
%             Results_P = R_indv(Find_patch,:);
%             confidence = sum(Results_P(:,9))/size(Results_P,1);
%         else
%             Find_patch = R_indv(:,5)==2 & R_indv(:,8)==-1;
%             Results_P = R_indv(Find_patch,:);
%             confidence = sum(Results_P(:,9))/size(Results_P,1);
%         end
% 
%         matrix2(sub,loc,condition) = confidence;
% 
%         clear Results_P
%         clear condition_mean 
%         clear confidence
%     end
%     end
% 
%     end
% 
%     for condition = 1:2 
% 
%         grandmatrix2(1,:,condition) = nanmean(matrix2(:,:,condition));
%         grandmatrix2(2,:,condition)= within_se(matrix2(:,:,condition),num_sub,3);
% 
%     end
% 
%     % calculate hit and CR collapsed across eccentricity levels
%     all_cong(1,:) = [mean(grandmatrix1(1,:,1)) mean(grandmatrix1(1,:,2))];
%     all_cong(2,:) = within_se([mean(matrix1(:,:,1),2) mean(matrix1(:,:,2),2)],num_sub,2);
%     all_incong(1,:) = [mean(grandmatrix2(1,:,1)) mean(grandmatrix2(1,:,2))];
%     all_incong(2,:) = within_se([mean(matrix2(:,:,1),2) mean(matrix2(:,:,2),2)],num_sub,2);
%     
%     % plot figures
%     % congruent
%     
% figure(out);
% subplot(2,3,5),errorbar(-5,all_cong(1,1),all_cong(2,1),'d','MarkerSize',8,...
%     'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1);
% hold on
% errorbar(-5, all_cong(1,2), all_cong(2,2),'o','MarkerSize',6,...
%     'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(1,:),'LineWidth',1);
% errorbar(-5,all_incong(1,1),all_incong(2,1),'d','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
%     colours(2,:),'Color',colours(2,:),'LineWidth',1);
% errorbar(-5,all_incong(1,2),all_incong(2,2),'o','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
%     colours(1,:),'Color',colours(1,:),'LineWidth',1);
% hold off
% ylabel('Confidence level');
% xlim([-7 11]),xticks([-5 0 6.48 9.16]);
% set(gca,'XTickLabel',{'All','0','6.50','9.20'},'FontSize',12);
% ylim([1 4]);
% legend('off');
% 
% 
% hold on
% errorbar([0 6.48 9.16],grandmatrix1(1,:,1),grandmatrix1(2,:,1),'d-','MarkerSize',6,...
%     'MarkerFaceColor',colours(2,:),'MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix1(1,:,2),grandmatrix1(2,:,2),'o-','MarkerSize',6,'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',...
%     colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix2(1,:,1),grandmatrix2(2,:,1),'d--','MarkerSize',7,...
%     'MarkerFaceColor','white','MarkerEdgeColor',colours(2,:),'Color',colours(2,:),'LineWidth',1,'Capsize',10);
% errorbar([0 6.48 9.16],grandmatrix2(1,:,2),grandmatrix2(2,:,2),'o--','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor',...
%     colours(1,:),'Color',colours(1,:),'LineWidth',1,'Capsize',10);
% plot([-7 11],[0.5 0.5],'k--');
% hold off
% legend({'Congruent Hit','Congruent CR','Incongruent Hit','Incongruent CR'});

% hypothesis 2 AUC
grandmatrix = zeros(num_sub,2);
for condition = 1:2
    
    for sub = 1:num_sub
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

% across eccentricities

%congruent

matrix4 = zeros(num_sub, 3);
location = [0 6.48 9.16];
for sub = 1:num_sub
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



%incongruent condition

clear indvN
clear indvP
clear indvN_loc
clear indvP_loc
clear population_m
clear indv_mean
clear indv_dff

matrix5 = zeros(num_sub, 3);

for sub = 1:num_sub
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

figure(out);
subplot(2,2,4),d = errorbar(-5,mean(grandmatrix(:,1)),std(grandmatrix(:,1))/sqrt(num_sub),'d','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
ylabel('AUC');
xlim([-7 11]),xticks([-5 0 6.48 9.16]);
set(gca,'XTickLabel',{'All','0','6.50','9.20'},'FontSize',12);
ylim([0.4 1]);
legend('off');
hold on
errorbar(-5,mean(grandmatrix(:,2)),std(grandmatrix(:,2))/sqrt(num_sub),'d','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.48 9.16],mean(matrix4),within_se(matrix4,num_sub,3),'d-','MarkerSize',6,...
    'MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
errorbar([0 6.48 9.16],mean(matrix5),within_se(matrix5,num_sub,3),'d--','MarkerSize',6,...
    'MarkerFaceColor','white','MarkerEdgeColor',colours(1,:),'Color',colours(2,:),'LineWidth',1);
plot([-7 11],[0.5 0.5],'k--');
legend({'Congruent','Incongruent'});
title('c)');
hold off

%% bar graph of hit, CR and AUC for hypo 1

colours = cbrewer('qual', 'Set1', 8); 
grandmatrix = zeros(2,2);
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
% [h,p,ci]= ttest(matrix1(:,1),matrix1(:,2));

% calculate standard errors
grandmatrix(1,:) = mean(matrix1);

for condition = 1:2
grandmatrix(2,condition) = std(matrix1(:,condition))/sqrt(15);
end

% AUC
matrix3 = zeros(15, 1);
for sub = 1:15
    indvN = Results(Results(:,1)==sub & Find_N,:);
    indvP = Results(Results(:,1)==sub & (Find_CAP|Find_IAP),:); % trial classification
    Confidence_N = indvN(:,9);
    Confidence_AP = indvP(:,9);
    
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
    matrix3(sub,1)= AUC;
    
clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;

    

end
grandmatrix(1,3) = mean(matrix3);
grandmatrix(2,3)= std(matrix3)/sqrt(15);

patch = categorical({'Hit','CR','AUC'});
patch = reordercats(patch,{'Hit','CR','AUC'});

subplot(1,2,1), b = bar(grandmatrix(1,:));
b.FaceColor = 'flat';
b.CData(1:3,:) = [colours(6,:);colours(5,:);colours(6,:)];
hold on
bar(3,grandmatrix(1,3),'EdgeColor',colours(5,:),'FaceColor',colours(6,:),'LineWidth',1.5);
errorbar([1 2 3],grandmatrix(1,:),grandmatrix(2,:),'k.','LineWidth',1.5);
set(gca,'XTickLabel',{'Hit','CR','AUC'},'YLim',[0.5 0.95],'FontName','Arial','FontSize',12);
sigstar({{'Hit','CR'}},0.001);
hold off
title('a)');

[h,p] = ttest(matrix1(:,1),matrix1(:,2));

%% hit, FA and AuC of hypo 2
% hit and FA
    matrix1 = zeros(15,2,2);
    location = [0 6.48 9.16];

  for congruency = 1:2
    for sub = 1:15
        if congruency ==1
            R_indv = Results(Results(:,1)==sub &(Find_Congruent_CP |Find_Congruent_IP),:);
        else
            R_indv = Results(Results(:,1)==sub &(Find_Incongruent_CP |Find_Incongruent_IP),:);
        end
    for condition = 1:2 
        if condition == 1 && congruency == 1
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        elseif condition == 2 && congruency == 1
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        elseif condition == 1 && congruency == 2
            Find_patch = R_indv(:,5)==3;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        elseif condition == 2 && congruency == 2
            Find_patch = R_indv(:,5)==2;
            Results_P = R_indv(Find_patch,:);
            accuracy = sum(Results_P(:,8)==1)/size(Results_P,1);
        end

        matrix1(sub,condition,congruency) = accuracy;

        clear Results_P
        clear condition_mean 
        clear accuracy
    end
    end
    
  end
  
grandmatrix1(1,1:2)= mean(matrix1(:,:,1),1);
grandmatrix1(1,3) = mean(matrix1(:,2,2),1);
grandmatrix1(1,4) = mean(matrix1(:,1,2));
grandmatrix1(2,1:2) = std(matrix1(:,:,1))/sqrt(15);
grandmatrix1(2,3) = std(matrix1(:,2,2))/sqrt(15);
grandmatrix1(2,4) = std(matrix1(:,1,2))/sqrt(15);
  
  
% AUC
matrix2 = zeros(15,2);
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

matrix2(sub,condition)= AUC;

clear Confidence_P
clear Confidence_A
clear Cumulative_NCounts;
clear Cumulative_APCounts;
clear Cumulative_Hit;
clear Cumulative_FA;
clear AUC;
    end

end

for condition = 1:2
grandmatrix1(1,condition + 4) = mean(matrix2(:,condition));
grandmatrix1(2,condition + 4) = std(matrix2(:,condition))/sqrt(15);
end


% plotting the values

% calculate axes
for i = 1:3
    if i == 1
        ax(2.*i-1) = i - 0.3;
        ax(2.*i) = i  + 0.3;
    elseif i == 2
        ax(2.*i-1) = i + 0.5 - 0.3;
        ax(2.*i) = i  + 0.5 + 0.3;
    elseif i == 3
        ax(2.*i-1) = i + 1 - 0.3;
        ax(2.*i) = i  + 1 + 0.3;
    end
end


subplot(1,2,2), d = bar(ax,grandmatrix1(1,:),'EdgeColor','none','BarWidth',0.9);
d.FaceColor = 'flat';
d.CData(1:2,:) = [colours(2,:);colours(1,:)];
d.CData(5,:) = colours(4,:);
hold on
bar(ax(3),grandmatrix1(1,3),'FaceColor','white','EdgeColor',colours(1,:),'LineWidth',1.5,'BarWidth',0.52);
bar(ax(4),grandmatrix1(1,4),'FaceColor','white','EdgeColor',colours(2,:),'LineWidth',1.5,'BarWidth',0.52);
bar(ax(6),grandmatrix1(1,6),'FaceColor','white','EdgeColor', colours(4,:),'LineWidth',1.5,'BarWidth',0.52);
errorbar(ax,grandmatrix1(1,:),grandmatrix1(2,:),'k.','LineWidth',1.2);
title('b)');
set(gca,'XTickLabel',{' ',' ',' ',' ',' ',' '},'YLim',[0.5 0.95],'XLim',[0.2 4.8],'FontName','Arial','FontSize',12);
hold off






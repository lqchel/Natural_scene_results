%%% this function plots the proportion of responses to each DxC category
%%% on each patch number within trial. 

% function [out] = temporal_DxC(Results_P, Results_A) % exp = 1/2, hypo = 1 for present, cong = 0 for hypo 1, 1 for congruent, 2 for incongruent
% % Results_P = Results(Find_Incongruent_IP,:);
% % Results_A = Results(Find_Incongruent_CP,:);
% location = [0 1 2];
% num_patch = 6;
% % 
% for a = 1:3
%     for t = 1:num_patch
%          for i = 1:9
%             pcgP_matrix(i,t,a) = sum(Results_P(:,9)== i-5 & Results_P(:,13) == location(a) & Results_P(:,3) == t)/...
%                 sum(Results_P(:,13) == location(a) & Results_P(:,3) == t);
%             pcgA_matrix(i,t,a) = sum(Results_A(:,9)== i-5 & Results_A(:,13) == location(a) & Results_A(:,3) == t)/...
%                 sum(Results_A(:,13) == location(a) & Results_A(:,3) == t);
% 
%          end
%     end
% end
% pcgP_matrix(5,:,:) = [];
% pcgA_matrix(5,:,:) = [];
% 
% 
% out = figure;
% for a = 1:3
%     figure(out);
%     subplot(2,3,a);
%     imagesc(flip(pcgP_matrix(:,:,a),1),[0 1]);
%     colorbar
%     set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
%     titles = {'Fovea','Periphery','Para-periphery'};
%     title(titles(a));
%     xticks(1:1:6);
%     
%     subplot(2,3,a+3);
%     imagesc(flip(pcgA_matrix(:,:,a),1),[0 1]);
%     colorbar
%     set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
%     xticks(1:1:6);
%     
% end
% 
% end

%function [out1,out2] = temporal_DxC(Results_P, Results_A) % exp = 1/2, hypo = 1 for present, cong = 0 for hypo 1, 1 for congruent, 2 for incongruent
 Results_P = Results(Find_Incongruent_IP,:);
 Results_A = Results(Find_Incongruent_CP,:);
location = [0 1 2];
num_patch = 6;
% 
for a = 1:3
    for t = 1:num_patch
         for i = 1:9
            pcgP_matrix(i,t,a) = sum(Results_P(:,9)== i-5 & Results_P(:,13) == location(a) & Results_P(:,3) == t)/...
                sum(Results_P(:,13) == location(a) & Results_P(:,3) == t);
            pcgA_matrix(i,t,a) = sum(Results_A(:,9)== i-5 & Results_A(:,13) == location(a) & Results_A(:,3) == t)/...
                sum(Results_A(:,13) == location(a) & Results_A(:,3) == t);

         end
    end
end
pcgP_matrix(5,:,:) = [];
pcgA_matrix(5,:,:) = [];

for patch_num = 1:6
    for loc = 1:3
        P_change_from_1(loc,patch_num) = pcgP_matrix(end,patch_num,loc) - pcgP_matrix(end,1,loc);
        A_change_from_1(loc,patch_num) = pcgA_matrix(end,patch_num,loc) - pcgA_matrix(end,1,loc);
    end
end

out1 = figure; % plot present
for a = 1:3
    figure(out1);
    subplot(2,3,a);
    imagesc(flip(pcgP_matrix(:,:,a),1),[0 1]);
    colorbar
    set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
    titles = {'Fovea','Periphery','Para-periphery'};
    title(titles(a));
    xticks(1:1:6);
    
    subplot(2,3,a+3)
    bar([1:6],P_change_from_1(a,:));
    xlim([0.5 7]);
    ylim([round(min(min(P_change_from_1)),1) 0]);
    set(gca,'FontSize',12,'FontName','Arial');
    


end

out2 = figure;
for a = 1:3
    figure(out2)
    subplot(2,3,a);
    imagesc(flip(pcgA_matrix(:,:,a),1),[0 1]);
    colorbar
    set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
    xticks(1:1:6);
    
    subplot(2,3,a+3)
    bar([1:6],P_change_from_1(:,:,a));
end

%end
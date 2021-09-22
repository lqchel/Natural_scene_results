%%% this function plots the proportion of responses to each DxC category
%%% on each patch number within trial. 

function [out] = temporal_DxC(Results_P, Results_A,change) %% change: 0 = plotting DxC proportion on each temporal position only, 1 = plotting DxC proportion change from patch position 1
% Results_P = Results(Find_Incongruent_IP,:);
% Results_A = Results(Find_Incongruent_CP,:);
% change = 1;
location = [0 1 2];
num_patch = 6;

if change == 0
    
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

%%% plot figure
out = figure;
for a = 1:3
    figure(out);
    subplot(2,3,a);
    imagesc(flip(pcgP_matrix(:,:,a),1),[0 1]);
    colorbar
    set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
    titles = {'Fovea','Periphery','Para-periphery'};
    title(titles(a));
    xticks(1:1:6);
    
    subplot(2,3,a+3);
    imagesc(flip(pcgA_matrix(:,:,a),1),[0 1]);
    colorbar
    set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
    xticks(1:1:6);
    
end

elseif change == 1
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

    for a = 1:3
        P_1_mat = pcgP_matrix(:,1,a) + zeros(8,num_patch); %% get proportions in each DxC category for patch number 1
        A_1_mat = pcgA_matrix(:,1,a) + zeros(8,num_patch);
        P_change_from_1(:,:,a) = pcgP_matrix(:,:,a)- P_1_mat; %% get proportion change from 1 for all DxC categories in all temporal positions ,and in all eccentricities
        A_change_from_1(:,:,a) = pcgA_matrix(:,:,a) - A_1_mat;
    end


    %%% plot figure
    out = figure;
    for a = 1:3
        figure(out);
        subplot(2,3,a);
        imagesc(flip(P_change_from_1(:,:,a),1),[-0.3 0]);
        colormap(flip(parula))
        colorbar
        set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
        titles = {'Fovea','Periphery','Para-periphery'};
        title(titles(a));
        xticks(1:1:6);

        subplot(2,3,a+3);
        imagesc(flip(A_change_from_1(:,:,a),1),[-0.3 0]);
        colormap(flip(parula))
        colorbar
        set(gca,'YTick',[1:1:8],'YTickLabel',{'4','3','2','1','-1','-2','-3','-4'},'FontSize', 12,'FontName','Arial');
        xticks(1:1:6);

    end

end

end

 
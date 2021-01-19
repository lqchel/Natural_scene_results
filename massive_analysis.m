%% get data
Results = importdata('Pooled Results.mat');
Results(:,9) = Results(:,9) - 0.5;
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

location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
eccentricity = zeros(length(Results),1)+ location1 + location2;

% load image list
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\congruent cropped'; 
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\incongruent cropped';
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);



%% original vs. modified

%congruent
delta_1 = zeros(80,3);
%size = zeros(80,1);
cong_significance = zeros(80,1);
for img = 3:82
%     current_cong= imread(fullfile(folder1,theFiles1(img-2).name));
%     current_incong= imread(fullfile(folder2,theFiles2(img-2).name));
    temp = Results(Results(:,2)==img & Find_Congruent_CP,9)-Results(Results(:,2)==img & Find_Congruent_IP,9);
    delta_1(img-2,1) = mean(temp);
    delta_1(img-2,2) = img;
    %delta_1(img-2,3) = unique(eccentricity(Results(:,2)==img & Find_Congruent_CP,1));
    delta_1(img-2,3) = sum(Results(:,2)==img & Find_Congruent_CP);
    [h,delta_1(img-2,4)] = ttest(temp,0); 
%     if img == 66
%     disp(Results(Results(:,2)==img & Find_Congruent_CP,9)-Results(Results(:,2)==img & Find_Congruent_IP,9));
%     break
%     end
    %size(img-2,1) = sum(sum(sum(current_cong-current_incong)))/(440^2);
     c(img-2,1) = unique(Results(Results(:,2)==img & Find_Congruent_CP,7));
    clear temp
end

%incongruent
delta_2 = zeros(80,3);
%condition_diff = zeros(80,1);
incong_significance = zeros(80,1);
for img = 3:82
    temp = Results(Results(:,2)==img & Find_Incongruent_IP,9)-Results(Results(:,2)==img & Find_Incongruent_CP,9);
    delta_2(img-2,1) = mean(temp);
    delta_2(img-2,2) = img;
    delta_2(img-2,3) = sum(Results(:,2)==img & Find_Incongruent_IP);
    %condition_diff(img-2,1) = delta_1(img-2,1)-delta_2(img-2,1);
    [h,delta_2(img-2,4)] = ttest(temp,0);
    in(img-2,1) = unique(Results(Results(:,2)==img & Find_Incongruent_IP,7));
    clear temp
end

cong_delta = table(delta_1(:,1),delta_1(:,2),delta_1(:,3),delta_1(:,4),'VariableNames',{'delta','img','participant','significance'});
incong_delta = table(delta_2(:,1),delta_2(:,2),delta_2(:,3),delta_2(:,4),'VariableNames',{'delta','img','participant','significance'});

%% compare delta between conditions
condition_significance = zeros(80,1);
for img = 3:82
    temp1 = Results(Results(:,2)==img & Find_Congruent_CP,9)-Results(Results(:,2)==img & Find_Congruent_IP,9);
    temp2 = Results(Results(:,2)==img & Find_Incongruent_IP,9)-Results(Results(:,2)==img & Find_Incongruent_CP,9);
%     cond_data = table([temp1(:,1);temp2(:,1)],[temp1(:,2);temp2(:,2)],'VariableNames',{'delta','congruence'});
    [h,condition_significance(img-2,1)]= ttest2(temp1,temp2);
    clear temp1
    clear temp2
end

%% plotting the lines
colours = cbrewer('qual', 'Set1', 8);
[Y1,edges] = histcounts(delta_1(:,1),80);
Y1_cumulative  = cumsum(Y1);
plot(edges,[Y1_cumulative 80],'LineWidth',1.2,'Color',colours(2,:));

hold on
[Y2,edges] = histcounts(delta_2(:,1),80);
Y2_cumulative  = cumsum(Y2);
plot(edges,[Y2_cumulative 80],'LineWidth',1.2,'Color',colours(1,:));
hold off
box off

xlabel(['\Delta','rcD×C (Original - Modified)']), ylabel('Cumulative counts');
set(gca,'FontName','Arial','FontSize',12);
xlim([-7 7]);
legend({'Congruent','Incongruent'},'Box','off','Location','northwest');

%% display images based on delta rcBxC
% load image list
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\congruent cropped'; 
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\MassiveReport_Exp2_QL\squareimage\incongruent cropped';
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);

[delta_1_x,rank1] = sort(delta_1(:,1),'descend');
% [delta_2_x,rank2] = sort(delta_2(:,1),'descend');


delta_1_x = [delta_1_x delta_1(rank1,2) delta_1(rank1,3)];
delta_2_x = [delta_2(rank1,1) delta_2(rank1,2) delta_2(rank1,3)];

delta_1_x = [delta_1_x(:,1:3) c(rank1,1)];
delta_2_x = [delta_2_x(:,1:3) in(rank1,1)];

% display (axis definition [left bottom width height])

%congruent
for i = 1:16
for p = 1:5
    current_cong= imresize(imread(fullfile(folder1,theFiles1(delta_1_x((i-1).*5+p,2)).name)),[150 150]); % load and resize image to square shape
    ax1= axes('Position',[0.015+(p-1).*0.19 0.52 0.2 0.2]); % fix axes position in figure window
    image(ax1,current_cong); % present image
    
    % mark out the critical objects
    hold on
    if delta_1_x((i-1).*5+p,4) == 1
        % outline the critical object on the image in red, using the x and y
        % coordinates of the object area
        plot(ax1,[1 50 50 1 1],[50 50 1 1 50],'r-','LineWidth',1.2); 
    elseif delta_1_x((i-1).*5+p,4) == 2
        plot(ax1,[50 100 100 50 50],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 3
        plot(ax1,[100 150 150 100 100],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 4
        plot(ax1,[1 50 50 1 1],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 5
        plot(ax1,[50 100 100 50 50],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 6
        plot(ax1,[100 150 150 100 100],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 7
        plot(ax1,[1 50 50 1 0],[150 150 100 100 150],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 8
        plot(ax1,[50 100 100 50 50],[150 150 100 100 150],'r-','LineWidth',1.2);
    else
        plot(ax1,[100 150 150 100 100],[150 150 100 100 150],'r-','LineWidth',1.2);
    end
    hold off
    axis square
    title({['d cong = ', num2str(round(delta_1_x((i-1).*5 + p,1)),1),'[',num2str(delta_1_x((i-1).*5 + p,3)),']'],...
        ['d incong = ', num2str(round(delta_2_x((i-1).*5 + p,1),1)),'[',num2str(15-delta_1_x((i-1).*5 + p,3)),']']});
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    
    current_incong= imresize(imread(fullfile(folder2,theFiles2(delta_1_x((i-1).*5+p,2)).name)),[150 150]);
    ax2= axes('Position',[0.015+(p-1).*0.19 0.3 0.2 0.2]);
    image(ax2,current_incong);
    
    % mark out the critical objects
    hold on
    if delta_1_x((i-1).*5+p,4) == 1
        plot(ax2,[1 50 50 1 1],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 2
        plot(ax2,[50 100 100 50 50],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 3
        plot(ax2,[100 150 150 100 100],[50 50 1 1 50],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 4
        plot(ax2,[1 50 50 1 1],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 5
        plot(ax2,[50 100 100 50 50],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 6
        plot(ax2,[100 150 150 100 100],[100 100 50 50 100],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 7
        plot(ax2,[1 50 50 1 1],[150 150 100 100 150],'r-','LineWidth',1.2);
    elseif delta_1_x((i-1).*5+p,4) == 8
        plot(ax2,[50 100 100 50 50],[150 150 100 100 150],'r-','LineWidth',1.2);
    else
        plot(ax2,[100 150 150 100 100],[150 150 100 100 150],'r-','LineWidth',1.2);
    end
    hold off
    axis square
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);

end

    if i < 10
        filename = ['cong_0', num2str(i),'.jpg'];
    else
        filename = ['cong_', num2str(i), '.jpg'];
    end
    saveas(gcf,filename);
clf;
end


%%
%incongruent
for i = 1:16

for p = 1:5
    current_cong= imresize(imread(fullfile(folder1,theFiles1(delta_2_x((i-1).*5+p,2)).name)),[150 150]);
    ax1= axes('Position',[0.015+(p-1).*0.19 0.52 0.2 0.2]);
    image(ax1,current_cong);
    axis square
    title({['delta = ', num2str(round(delta_2_x((i-1).*5 + p,1),1)), ', ', num2str(round(delta_2_x((i-1).*5 + p,3),1)), ' dva']});
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    current_incong= imresize(imread(fullfile(folder2,theFiles2(delta_2_x((i-1).*5+p,2)).name)),[150 150]);
    ax2= axes('Position',[0.015+(p-1).*0.19 0.3 0.2 0.2]);
    image(ax2,current_incong);
    axis square
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);

end

    if i < 10
        filename = ['incong_0', num2str(i),'.jpg'];
    else
        filename = ['incong_', num2str(i), '.jpg'];
    end
    saveas(gcf,filename);
    clf;
end

%% create scatterplot
c_n_i_n = cong_delta.significance > 0.05 & incong_delta.significance > 0.05;
c_n_i_y = cong_delta.significance > 0.05 & incong_delta.significance < 0.05;
c_y_i_n = cong_delta.significance < 0.05 & incong_delta.significance > 0.05;
c_y_i_y = cong_delta.significance < 0.05 & incong_delta.significance < 0.05;

temp_colour = cbrewer('qual', 'Pastel2', 8);
sz = 50;

scatter(cong_delta.delta(c_n_i_n,1),incong_delta.delta(c_n_i_n,1),[],temp_colour(8,:),'filled','LineWidth',1.2)
hold on
scatter(cong_delta.delta(c_n_i_y,1),incong_delta.delta(c_n_i_y,1),[],'red','filled','LineWidth',1.2);
scatter(cong_delta.delta(c_y_i_n,1),incong_delta.delta(c_y_i_n,1),[],'blue','filled','LineWidth',1.2);
scatter(cong_delta.delta(c_y_i_y,1),incong_delta.delta(c_y_i_y,1),[],'black','filled','LineWidth',1.2);
scatter(cong_delta.delta(condition_significance < 0.05,1),incong_delta.delta(condition_significance < 0.05,1),85,'ks','LineWidth',0.8);
hold off

xlabel(['\Delta' 'rcDxC in congruent images']), ylabel(['\Delta' 'rcDxC in incongruent images']);
set(gca,'FontName','Arial','FontSize',12);




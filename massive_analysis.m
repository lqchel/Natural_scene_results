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


%% original vs. modified

%congruent
delta_1 = zeros(80,3);
for img = 3:82
    delta_1(img-2,1) = median(Results(Results(:,2)==img & Find_Congruent_CP,9)-Results(Results(:,2)==img & Find_Congruent_IP,9));
    delta_1(img-2,2) = img;
    delta_1(img-2,3) = unique(eccentricity(Results(:,2)==img & Find_Congruent_CP,1));
    c(img-2,1) = unique(Results(Results(:,2)==img & Find_Congruent_CP,7));
end

%incongruent
delta_2 = zeros(80,3);
for img = 3:82
    delta_2(img-2,1) = median(Results(Results(:,2)==img & Find_Incongruent_IP,9)-Results(Results(:,2)==img & Find_Incongruent_CP,9));
    delta_2(img-2,2) = img;
    delta_2(img-2,3) = unique(eccentricity(Results(:,2)==img & Find_Incongruent_IP,1));
    in(img-2,1) = unique(Results(Results(:,2)==img & Find_Incongruent_IP,7));
end

[delta_1_x,rank1] = sort(delta_1(:,1),'descend');
[delta_2_x,rank2] = sort(delta_2(:,1),'descend');

delta_1_x = [delta_1_x delta_1(rank1,2) delta_1(rank1,3)];
delta_2_x = [delta_2_x delta_2(rank2,2) delta_2(rank2,3)];

delta_1_x = [delta_1_x(:,1:3) c(rank1,1)];
delta_2_x = [delta_2_x(:,1:3) in(rank2,1)];
%% plotting the lines
colours = cbrewer('qual', 'Set1', 8);
[Y1,edges] = histcounts(delta_1(:,1),80);
Y1_cumulative  = cumsum(Y1);
subplot(2,1,1),plot(edges,[Y1_cumulative 80],'LineWidth',1.2,'Color',colours(2,:));
hold on

[Y2,edges] = histcounts(delta_2(:,1),80);
Y2_cumulative  = cumsum(Y2);
plot(edges,[Y2_cumulative 80],'LineWidth',1.2,'Color',colours(1,:));

hold off
box off

xlabel('?rcD×C (Original - Modified)'), ylabel('Cumulative counts');
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

top = 1:3;
mid_1 = 47:50;
mid_2 = 97:100;
bottom = 147:150;

% display (axis definition [left bottom width height])

%congruent
for i = 1:1
for p = 1:5
    current_cong= imresize(imread(fullfile(folder1,theFiles1(delta_1_x((i-1).*5+p,2)).name)),[150 150]);
    if delta_1_x((i-1).*5+p,4) == 1
        current_cong([top mid_1],[top mid_1],1) = uint8(255);
        current_cong([top mid_1],[top mid_1],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 2
        current_cong([top mid_1],[mid_1 mid_2],1) = uint8(255);
        current_cong([top mid_1],[mid_1 mid_2],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 3
        current_cong([top mid_1],[mid_2 bottom],1) = uint8(255);
        current_cong([top mid_1],[mid_1 bottom],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 4
        current_cong([mid_1 mid_2],[top mid_1],1) = uint8(255);
        current_cong([mid_1 mid_2],[top mid_1],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 5
        current_cong([mid_1 mid_2],[mid_1 mid_2],1) = uint8(255);
        current_cong([mid_1 mid_2],[mid_1 mid_2],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 6
        current_cong([mid_1 mid_2],[mid_2 bottom],1) = uint8(255);
        current_cong([mid_1 mid_2],[mid_2 bottom],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 7
        current_cong([mid_2 bottom],[top mid_1],1) = uint8(255);
        current_cong([mid_2 bottom],[top mid_1],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 8
        current_cong([mid_2 bottom],[mid_1 mid_2],1) = uint8(255);
        current_cong([mid_2 bottom],[mid_1 mid_2],2:3) = uint8(0);
    else
        current_cong([mid_2 bottom],[mid_2 bottom],1) = uint8(255);
        current_cong([mid_2 bottom],[mid_2 bottom],2:3) = uint8(0);
    end
    
    ax1= axes('Position',[0.015+(p-1).*0.19 0.52 0.2 0.2]);
    image(ax1,current_cong);
    axis square
    title({['\delta = ', num2str(delta_1_x((i-1).*5 + p,1)), ', ', num2str(delta_1_x((i-1).*5 + p,3)), ' dva']});
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);
    
    current_incong= imresize(imread(fullfile(folder2,theFiles2(delta_1_x((i-1).*5+p,2)).name)),[150 150]);
    
    if delta_1_x((i-1).*5+p,4) == 1
        current_cong([top mid_1],[top mid_1],1) = uint8(255);
        current_cong([top mid_1],[top mid_1],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 2
        current_cong([top mid_1],[mid_1 mid_2],1) = uint8(255);
        current_cong([top mid_1],[mid_1 mid_2],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 3
        current_cong([top mid_1],[mid_2 bottom],1) = uint8(255);
        current_cong([top mid_1],[mid_1 bottom],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 4
        current_cong([mid_1 mid_2],[top mid_1],1) = uint8(255);
        current_cong([mid_1 mid_2],[top mid_1],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 5
        current_cong([mid_1 mid_2],[mid_1 mid_2],1) = uint8(255);
        current_cong([mid_1 mid_2],[mid_1 mid_2],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 6
        current_cong([mid_1 mid_2],[mid_2 bottom],1) = uint8(255);
        current_cong([mid_1 mid_2],[mid_2 bottom],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 7
        current_cong([mid_2 bottom],[top mid_1],1) = uint8(255);
        current_cong([mid_2 bottom],[top mid_1],2:3) = uint8(0);
    elseif delta_1_x((i-1).*5+p,4) == 8
        current_cong([mid_2 bottom],[mid_1 mid_2],1) = uint8(255);
        current_cong([mid_2 bottom],[mid_1 mid_2],2:3) = uint8(0);
    else
        current_cong([mid_2 bottom],[mid_2 bottom],1) = uint8(255);
        current_cong([mid_2 bottom],[mid_2 bottom],2:3) = uint8(0);
    end
    ax2= axes('Position',[0.015+(p-1).*0.19 0.3 0.2 0.2]);
    image(ax2,current_incong);
    axis square
    set(gca,'FontName','Arial','FontSize',8,'FontWeight','normal','Box','off','XColor','none','YColor','none');
    xticks([]),yticks([]);

end

%     if i < 10
%         filename = ['cong_0', num2str(i),'.jpg'];
%     else
%         filename = ['cong_', num2str(i), '.jpg'];
%     end
%     saveas(gcf,filename);
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
    title({['delta = ', num2str(delta_2_x((i-1).*5 + p,1)), ', ', num2str(delta_2_x((i-1).*5 + p,3)), ' dva']});
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
[A,img_num1] = sort(delta_1_x(:,2));
[B,img_num2] = sort(delta_2_x(:,2));

scatter(delta_1_x(img_num1,1),delta_2_x(img_num2,1),'filled');
xlabel(['\delta' 'rcDxC in congruent condition']), ylabel(['\delta' 'rcDxC in incongruent condition']);
set(gca,'FontName','Arial','FontSize',12);




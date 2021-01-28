%%%% This code finds and cuts the smallest rectangle area for the critical
%%%% objects in the sample congruent and incongruent image pairs.

% load image
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\congruent'; 
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incongruent';
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);

quantile_ranges = [0.7 1; 0.8 1; 0.9 1]; % set possible quantile ranges

for i = 1:20
    % read sample image
    cong = imread(fullfile(folder1,theFiles1(i).name)); 
    incong = imread(fullfile(folder2,theFiles2(i).name));

     % compute difference matrix between congruent and incongruent images
    difference = sum(abs(cong-incong),3);
    temp = difference(difference>0);

    for q = 1:3
        current_quantile = quantile(temp,quantile_ranges(q,:));
        [row,column] = find(difference > current_quantile(:,1) & difference < current_quantile(:,2)); % find the area of difference, by returning the row and column indices of pixels within that area 

        % mark critical object on congruent
        ax1= axes('Position',[0.015+(q-1).*0.3 0.52 0.3 0.3]); 
        image(ax1,cong);
        axis square

        hold on
        plot(ax1,[min(column) max(column) max(column) min(column) min(column)],...
            [max(row) max(row) min(row) min(row) max(row)],'r-','LineWidth',1.4);
        hold off

        set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
        xticks([]),yticks([]);
        title(['[' num2str(quantile_ranges(q,1)) ' ' num2str(quantile_ranges(q,2)) ']']);

        % mark critical object on incongruent
        ax2= axes('Position',[0.015+(q-1).*0.3 0.2 0.3 0.3]);
        image(ax2,incong);
        axis square

        hold on
        plot(ax2,[min(column) max(column) max(column) min(column) min(column)],...
            [max(row) max(row) min(row) min(row) max(row)],'r-','LineWidth',1.4);
        hold off
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal','Box','off','XColor','none','YColor','none');
        xticks([]),yticks([]);

    end

% save image

if i < 10
    filename = ['img_0', num2str(i),'.jpg'];
else
    filename = ['img_', num2str(i), '.jpg'];
end
saveas(gcf,filename);

if i<20
 clf;
end
clear temp
end

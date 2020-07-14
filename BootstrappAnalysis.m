% this script runs the bootstrap analysis of null patch d x c vs modified
% patch dxc, based on RGB, RGB, and RGB differences

%% load data
load('All absent decisions.mat');
%% bootstrapping
location = [0 6.50 9.20];

btp_means = zeros(200,3);

for e = 1:3

% create histogram of modified patch RGB
h = histogram(data_allab.RGB(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)),20);

% bin up the null patch responses, based on the modified patch histogram
Null_bin = cell(1,h.NumBins);

for i = 1:h.NumBins
Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.RGB >= h.BinEdges(i)& data_allab.RGB <= h.BinEdges(i+1)...
    & data_allab.Eccentricity == location(e),:);
end

% bootstrapping

btp_trials = cell(200,20);


for n = 1:200
    for c = 1:h.NumBins
        column_mean = zeros(1,h.NumBins);
        if h.Values(c) == 0
            column_mean(c)= 0;
        else
            rand_index = randperm(height(Null_bin{:,c})); 
            old_sq = Null_bin{:,c};
            new_sq = old_sq(rand_index,:);
            selected = new_sq(1:h.Values(:,c),:);
            btp_data{n,c} = selected;
            column_mean(c) = mean(selected.Decision);
        end
    end
    btp_means_RGB(n,e) = mean(column_mean(column_mean~=0));
end    
clear h;
end
clear e;

%% t tests

modified_original = data_allab.Decision(data_allab.PatchType~=0);
modified_median = median(modified_original);

all_means(:,:,1) = btp_means_Luminance;
all_means(:,:,2) = btp_means_Contrast;
all_means(:,:,3) = btp_means_RGB;

all_means = mean(all_means,3);

% not separating eccentricity
across_ecc = mean(all_means,2);
[h1,p1] = ttest(across_ecc,modified_median);

clear modified_median
clear modified_original
% separating eccentricity;

for e = 1:3
    modified_original = data_allab.Decision(data_allab.PatchType~=0& data_allab.Eccentricity == location(e));
    modified_median = median(modified_original);
    [h2,p2] = ttest(all_means(:,e),modified_median,'alpha',0.017);
    ecc_t(e) = p2;
    clear p2
end



%% plot histograms
for e = 1:3

    
subplot(3,2,2.*e-1), histogram(data_allab.RGB(data_allab.PatchType == 0 & data_allab.Eccentricity == location(e)));
hold on
histogram(data_allab.RGB(data_allab.PatchType ~= 0 & data_allab.Eccentricity == location(e)));
hold off
set(gca,'FontName','Arial','FontSize',12);

ylabel(['Eccentricity = ' num2str(location(e))]);
if e == 1
    title('Raw distribution'),legend({'Null','Modified'});
elseif e == 3
    xlabel('RGB Difference');
end


% q1 = quantile(btp_means,[0.025 0.25 0.75 0.975]);
modified_original = data_allab.Decision(data_allab.PatchType~=0& data_allab.Eccentricity == location(e));
% btp_means_interval = btp_means(btp_means<= q1(4)& btp_means>=q1(1));
modified_median = median(modified_original);

subplot(3,2,2.*e), histogram(btp_means(:,e));
hold on
d = gca;
limit = d.YLim;
plot([modified_median modified_median],[0 limit(2)],'r-','LineWidth',1);
hold off

set(gca,'FontName','Arial','FontSize',12,'XLim',[-4.5 4.5]);
xticks([-4:1:4]);

if e == 1
    title('Bootstrap results');
elseif e == 3
    xlabel('Decision x confidence');
end

clear modified_original
clear modified_median
end

%% not separating eccentricity

btp_means = zeros(200,1);

% create histogram of modified patch RGB
h = histogram(data_allab.RGB(data_allab.PatchType ~= 0),20);

% bin up the null patch responses, based on the modified patch histogram
Null_bin = cell(1,h.NumBins);

for i = 1:h.NumBins
Null_bin{:,i} = data_allab(data_allab.PatchType == 0 & data_allab.RGB >= h.BinEdges(i)& data_allab.RGB <= h.BinEdges(i+1),:);
end

% bootstrapping

btp_trials = cell(200,20);


for n = 1:200
    for c = 1:h.NumBins
        column_mean = zeros(1,h.NumBins);
        if h.Values(c) == 0
            column_mean(c)= 0;
        else
            rand_index = randperm(height(Null_bin{:,c})); 
            old_sq = Null_bin{:,c};
            new_sq = old_sq(rand_index,:);
            selected = new_sq(1:h.Values(:,c),:);
            btp_data{n,c} = selected;
            column_mean(c) = mean(selected.Decision);
        end
    end
    btp_means(n) = mean(column_mean(column_mean~=0));
end    
clear h;


% plot histograms


    
subplot(1,2,1), histogram(data_allab.RGB(data_allab.PatchType == 0));
hold on
histogram(data_allab.RGB(data_allab.PatchType ~= 0));
hold off
set(gca,'FontName','Arial','FontSize',12);
legend({'Null','Modified'});
title('Raw distribution');
ylabel('Count');
xlabel('RGB difference');
ylim([0 1100]);

% q1 = quantile(btp_means,[0.025 0.25 0.75 0.975]);
modified_original = data_allab.Decision(data_allab.PatchType~=0);
% btp_means_interval = btp_means(btp_means<= q1(4)& btp_means>=q1(1));
modified_median = median(modified_original);

subplot(1,2,2), histogram(btp_means);
hold on
d = gca;
limit = d.YLim;
plot([modified_median modified_median],[0 limit(2)],'r-','LineWidth',1);
hold off
set(gca,'FontName','Arial','FontSize',12,'XLim',[-4.5 4.5]);
xticks([-4:1:4]);
title('Bootstrap results');
ylabel('Count');
xlabel('Decision x confidence');






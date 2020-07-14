% plotting individual results ------------------------------
colours = [cbrewer('qual', 'Set1', 9); cbrewer('qual', 'Dark2', 8)];
for sub = 1:num_sub
  subplot(1,2,1),plot([0 6.48 9.16],matrix1(sub,:,1),'d-','Color',colours(sub,:),'MarkerSize',8,'LineWidth',0.8),title('Hit');
    xticks([0 6.48 9.16]),set(gca,'XTickLabel',{'0','6.50','9.20'},'FontSize',12);
    xlabel('Eccentricity'),ylabel('%Accurate judgments'),ylim([0 1]);
  hold on
  if sub == num_sub
      hold off
  end
  
  subplot(1,2,2),plot([0 6.48 9.16],matrix1(sub,:,2),'d-','Color',colours(sub,:),'MarkerSize',8,'LineWidth',0.8),title('CR');
    xticks([0 6.48 9.16]),set(gca,'XTickLabel',{'0','6.50','9.20'},'FontSize',12);
    xlabel('Eccentricity'),ylabel('%Accurate judgments'),ylim([0 1]);
    
    hold on
    if sub == num_sub
        hold off
    end
%     subplot(2,2,2),title('Congruent CR'),plot([0 6.48 9.16],matrix1(sub,:,2),'d-','MarkerSize',6,'LineWidth',0.5);
%     subplot(2,2,3),title('Incongruent Hit'),plot([0 6.48 9.16],matrix2(sub,:,1),'d-','MarkerSize',6,'LineWidth',0.5);
%     subplot(2,2,4),title('Incongruent CR'),plot([0 6.48 9.16],matrix2(sub,:,2),'d-','MarkerSize',6,'LineWidth',0.5);
end

patch = categorical({'present','absent'});
patch = reordercats(patch,{'present','absent'});
out = figure;
subplot(1,2,1), bar(patch,grandmatrix(1,:),'BarWidth',0.7), ylabel('%Accurate judgments');
hold on
errorbar(patch,grandmatrix(1,:),grandmatrix(2,:),'k.','LineWidth',1);
hold off

% ----------------------------------------------------------

% plotting individual results ------------------------------

colours = [cbrewer('qual', 'Set1', 9); cbrewer('qual', 'Dark2', 8)];
for sub = 1:num_sub
  subplot(2,2,1),plot([0 6.48 9.16],matrix1(sub,:,1),'d-','Color',colours(sub,:),'MarkerSize',8,'LineWidth',0.8),title('Congruent Hit');
    xticks([0 6.48 9.16]),set(gca,'XTickLabel',{'0','6.50','9.20'},'FontSize',12);
    xlabel('Eccentricity'),ylabel('%Accurate judgments'),ylim([0 1]);
  hold on
  if sub == num_sub
      hold off
  end
  
  subplot(2,2,2),plot([0 6.48 9.16],matrix1(sub,:,2),'d-','Color',colours(sub,:),'MarkerSize',8,'LineWidth',0.8),title('Congruent CR');
    xticks([0 6.48 9.16]),set(gca,'XTickLabel',{'0','6.50','9.20'},'FontSize',12);
    xlabel('Eccentricity'),ylabel('%Accurate judgments'),ylim([0 1]);
    
    hold on
    if sub == num_sub
        hold off
    end

    subplot(2,2,3),plot([0 6.48 9.16],matrix2(sub,:,1),'d-','Color',colours(sub,:),'MarkerSize',8,'LineWidth',0.8),title('Incongruent Hit');
        xticks([0 6.48 9.16]),set(gca,'XTickLabel',{'0','6.50','9.20'},'FontSize',12);
        xlabel('Eccentricity'),ylabel('%Accurate judgments'),ylim([0 1]);
      hold on
      if sub == num_sub
          hold off
      end
  
  subplot(2,2,4),plot([0 6.48 9.16],matrix2(sub,:,2),'d-','Color',colours(sub,:),'MarkerSize',8,'LineWidth',0.8),title('Incongruent CR');
    xticks([0 6.48 9.16]),set(gca,'XTickLabel',{'0','6.50','9.20'},'FontSize',12);
    xlabel('Eccentricity'),ylabel('%Accurate judgments'),ylim([0 1]);
    
    hold on
    if sub == num_sub
        hold off
    end

end


% ----------------------------------------------------------


%confidence plot for regan

grandmatrix = zeros(2,2);
matrix1 = zeros(num_sub,2);

for sub = 1:num_sub
    R_indv = Results(Results(:,1)==sub &(Find_IAP|Find_CAP|Find_N),:);
    R_indv(:,9) = abs(R_indv(:,9));
for condition = 1:2
    if condition == 1
        Find_patch = R_indv(:,5)~=1;
        Results_P = R_indv(Find_patch,:);
        confidence = sum(Results_P(:,9))/size(Results_P,1);
    else
        Find_patch = R_indv(:,5)==1;
        Results_P = R_indv(Find_patch,:);
        confidence = sum(Results_P(:,9))/size(Results_P,1);
    end
    
    matrix1(sub,condition) = confidence;
    
    clear Results_P
    clear condition_mean 
    clear confidence
    clear Find_patch
end
end

% calculate standard errors
grandmatrix(1,:) = mean(matrix1);


for condition = 1:2
grandmatrix(2,condition) = std(matrix1(:,condition))/sqrt(num_sub);
end

patch = categorical({'present','absent'});
patch = reordercats(patch,{'present','absent'});
subplot(1,2,2), bar(patch,abs(grandmatrix(1,:)),'BarWidth',0.7), ylabel('%Accurate judgments');
hold on
errorbar(patch,abs(grandmatrix(1,:)),grandmatrix(2,:),'k.','LineWidth',1);
hold off
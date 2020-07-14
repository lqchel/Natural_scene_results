criterion = {'No exclusion','Mean +/- 3SD', 'Mean +/- 2.5SD', 'Mean +/- 2SD'};
y = [0.62 0.65 0.84 0.87];
d = categorical(y,y,criterion,'Ordinal',true);
plot(d,y,'bo-','LineWidth',1.5,'MarkerSize',5), ylim([0 1]),...
   xlabel('Exclusion criterion, from lenient to strict'), ylabel('AUC for hypothesis 1');

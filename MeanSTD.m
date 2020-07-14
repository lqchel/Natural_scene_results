clear all;

data1 = importdata('OnlinePilotData.mat');
data2 = importdata('Pooled Results 2.mat');
[m1,sd1] = mean_std(data1, 10);
[m2,sd2] = mean_std(data2,15);

plot(m1,sd1,'bo','MarkerFaceColor','blue');
xlabel('Mean'),ylabel('SD');
hold on
plot(m2,sd2,'ro','MarkerFaceColor','red');
hold off
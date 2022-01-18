% ori = [0 0+0.125i 0+0.25i 0+0.375i 0+0.5i 0+1i 0+1.5i 0+2i 0+2.5i 0+3.24i ...
%     3.24+3.24i 3.24-3.24i 0-3.24i 0-2.5i 0-2i 0-1.5i 0-1i 0-0.5i 0-0.375i ...
%     0-0.25i 0-0.125i 0];
% ori = 10.47.*ori;
% cor_ori = log(1+ori);
% ori_2 = [0+9.72i 3.24+9.72i 3.24+3.24i 0+3.24i 0+9.72i];
% ori_2 = ori_2.*10.47;
% cor_ori_2 = log(1+ori_2);
% ori_3 = [3.24+9.72i 9.72+9.72i 9.72+3.24i 3.24+3.24i 3.24+9.72i];
% ori_3 = 10.47.*ori_3;
% cor_ori_3 = log(1+ori_3);
% 
% fill(cor_ori,'b-','LineWidth',1.2);
% hold on
% fill(cor_ori_2,'r-','LineWidth',1.2);
% fill(cor_ori_3,'m-','LineWidth',1.2);
% hold off
% legend({'Eccentricity = 0','Eccentricity = 6.48','Eccentricity = 9.16'});
% axis square

% polyfill(real(cor_ori),imag(cor_ori))
% 
% ans =
% 
%     9.5848
% 
% polyfill(real(cor_ori_2),imag(cor_ori_2))
% 
% ans =
% 
%     0.5668
% 
% polyfill(real(cor_ori_3),imag(cor_ori_3))
% 
% ans =
% 
%     0.5016
% 
%% central
x = [zeros(1,length(0:0.01:3.24)) 0:0.02:3.24 3.24+zeros(1,length(0:0.02:6.48)) -1.*[-3.24:0.02:0] zeros(1,length(0:0.02:3.24)) 0];
y = [0:0.01:3.24 3.24+zeros(1,length(0:0.02:3.24)) -1.*[-3.24:0.02:3.24] -3.24+zeros(1,length(0:0.02:3.24)) -3.24:0.02:0 0];
ori = 10.47.*complex(x,y);
cor_ori = log(ori+1);

%%%


%% para_peripheral
y1 = y+6.48;
ori_2 = 10.47 .*complex(x,y1);
cor_ori_2 = log(ori_2+1);


%% peripheral
x1 = [3.24+zeros(1,length(0:0.02:3.24)) 3.24:0.02:9.72 9.72+zeros(1,length(0:0.02:6.48)) -1.*[-9.72:0.02:-3.24] 3.24+zeros(1,length(0:0.02:3.24))];
y2 = [6.48:0.02:9.72 9.72+zeros(1,length(0:0.02:6.48)) -1.*[-9.72:0.02:-3.24] 3.24+zeros(1,length(0:0.02:6.48)) 3.24:0.02:6.48];
 ori_3 = 10.47 .*complex(x1,y2);

cor_ori_3 = log(ori_3+1);

%% para_peripheral 2 (full patch on one side of visual field)

y3 = y2-6.48;
ori_4 = 10.47 .*complex(x1,y3);

cor_ori_4 = log(ori_4+1);

%% fill projection
colours = cbrewer('qual','Pastel2',8);
fill(real(cor_ori),imag(cor_ori),colours(3,:));
hold on
fill(real(cor_ori_2),imag(cor_ori_2),colours(2,:));
fill(real(cor_ori_3),imag(cor_ori_3),colours(6,:));
fill(real(cor_ori_4),imag(cor_ori_4),colours(2,:));
fill(-1.*real(cor_ori),imag(cor_ori),colours(3,:));
fill(-1.*real(cor_ori_2),imag(cor_ori_2),colours(2,:));
fill(-1.*real(cor_ori_3),imag(cor_ori_3),colours(6,:));
fill(-1.*real(cor_ori_2),-1.*imag(cor_ori_2),colours(2,:));
fill(-1.*real(cor_ori_3),-1.*imag(cor_ori_3),colours(6,:));
fill(real(cor_ori_2),-1.*imag(cor_ori_2),colours(2,:));
fill(real(cor_ori_3),-1.*imag(cor_ori_3),colours(6,:));
fill(-1.*real(cor_ori_4),imag(cor_ori_4),colours(2,:));
hold off
legend({'Eccentricity = 0','Eccentricity = 6.48','Eccentricity = 9.16'});
xlim([-5 5]),set(gca,'FontSize',12)

%% calculate fill of projection

fill_1 = polyarea(real(cor_ori),imag(cor_ori));
fill_2 = 2.*polyarea(real(cor_ori_2),imag(cor_ori_2));
fill_3 = polyarea(real(cor_ori_3),imag(cor_ori_3));


%% new experiment

% calculate half visual field size

half_field = fill_1 + 2.*fill_2 + 2.*fill_3; 

% dividing visual field

dv_x = real(cor_ori);
dv_x = dv_x(1:length(cor_ori)/2);
dv_y = imag(cor_ori);
dv_y = dv_y(1:length(cor_ori)/2);

num = reshape([1:1:length(dv_x)], 1, length(dv_x));

% index = find(polyarea([dv_x(1:num) dv_x(num)],[dv_y(1:num) 0]) == half_field/2);


window_min = [2.5 3 3.5];
window_max = [3 3.5 4];

for window = 1:3
    
dv_x1 = dv_x(dv_x>= window_min(window) & dv_x <= window_max(window));
% dv_y1 = dv_y(dv_x>= window_min(window) & dv_x <= window_max(window));

    for i = 1:length(dv_x1)
        if window == 1
            diff_1(i) = abs(polyarea([dv_x(dv_x<=dv_x1(i)) dv_x1(i)],[dv_y(1:sum(dv_x<=dv_x1(i))) 0]) - half_field/4);
            x_index1(i) = dv_x1(i);   
        elseif window == 2
            diff_2(i) = abs(polyarea([dv_x(dv_x<=dv_x1(i)) dv_x1(i)],[dv_y(1:sum(dv_x<=dv_x1(i))) 0]) - half_field/4);
            x_index2(i) = dv_x1(i); 
        elseif window == 3
            diff_3(i) = abs(polyarea([dv_x(dv_x<=dv_x1(i)) dv_x1(i)],[dv_y(1:sum(dv_x<=dv_x1(i))) 0]) - half_field/4);
            x_index3(i) = dv_x1(i); 
            
        end
    end
clear dv_x1
clear dv_y1
end

min(diff_1),x_index1(diff_1 == min(diff_1))

min(diff_2),x_index2(diff_2 == min(diff_2)) % line of dividing the shape sits at x = x_index1

min(diff_3),x_index3(diff_3 == min(diff_3))

smallest = x_index1(diff_1 == min(diff_1));
smallest_y = dv_y(dv_x == smallest);

%% convert it back to visual field space

x5 = [dv_x(dv_x <= smallest) flip(dv_x(dv_x <= smallest))];
y5 = [dv_y(1:sum(dv_x<=smallest)) -1.*flip(dv_y(1:sum(dv_x<=smallest)))];

y6(1,:) = -smallest_y:0.01:smallest_y;
x6 = smallest + zeros(1,length(y6));


converted = exp(complex(x6,y6));
converted = reshape(converted,length(x6),1);

% wholefield_x = 4.28.*[real(converted-1); -1.*(real(converted-1))];
% wholefield_y = [imag(converted-1); flip(imag(converted-1))].*4.28;

wholefield_x_right = 4.28.*real(converted-1); 
wholefield_x_left = -4.28.*(real(converted-1));
wholefield_y_right = 4.28.*imag(converted-1); 
wholefield_y_left = flip(imag(converted-1)).*4.28;


colours = cbrewer('qual','Pastel2',8);
subplot(1,2,1),fill(wholefield_x_left,wholefield_y_left,colours(2,:));
axis square
hold on
fill(wholefield_x_right,wholefield_y_right,colours(3,:));
plot([-440 -440 440 440 -440], [-440 440 440 -440 -440],'k--','LineWidth',1);
legend({'Left foveal area','Right foveal area','Image boundary'});
xlabel('Pixels'), ylabel('Pixels');
title('Centre and periphery with equal cortical area representation');
hold off
set(gca,'FontSize','12','FontName','Arial');

%% also how this is presented in the cortical space
% x7 = [-9.72 + zeros(1,length(-9.72:0.02:9.72)) -9.72:0.02:0 0+zeros(1,length(-9.72:0.02:9.72),1) -1.*[0:0.02:9.72]];
% y7 = [-9.72:0.02:9.72 9.72+zeros(1,length(-9.72:0.02:0)) -1.*[-9.72:0.02:9.72] -9.72+zeros(1,length(0:0.02:9.72))];
% 
% edge_ori = complex(x7,y7).*10.48;
% edge_cor = log(edge_ori +1);
% x7 = [real(cor_ori) real(cor_ori_2) real(cor_ori_3) real(cor_ori_4) real(cor_ori_2) real(cor_ori_3)];
% y7 = [imag(cor_ori) imag(cor_ori_2) imag(cor_ori_3) imag(cor_ori_4) -1.*imag(cor_ori_2) -1.*imag(cor_ori_3)];
% 
% k = boundary(reshape(x7,length(x7),1),reshape(y7,length(y7),1),1);

subplot(1,2,2),plot(real(cor_ori),imag(cor_ori),'k--','LineWidth',1);
hold on
plot(real(cor_ori_2),imag(cor_ori_2),'k--','LineWidth',1);
plot(real(cor_ori_3),imag(cor_ori_3),'k--','LineWidth',1);
plot(real(cor_ori_4),imag(cor_ori_4),'k--','LineWidth',1);
plot(real(cor_ori_2),-1.*imag(cor_ori_2),'k--','LineWidth',1);
plot(real(cor_ori_3),-1.*imag(cor_ori_3),'k--','LineWidth',1);
plot(-1.*real(cor_ori_4),imag(cor_ori_4),'k--','LineWidth',1);
plot(-1.*real(cor_ori),imag(cor_ori),'k--','LineWidth',1);
plot(-1.*real(cor_ori_2),imag(cor_ori_2),'k--','LineWidth',1);
plot(-1.*real(cor_ori_3),imag(cor_ori_3),'k--','LineWidth',1);
plot(-1.*real(cor_ori_2),-1.*imag(cor_ori_2),'k--','LineWidth',1);
plot(-1.*real(cor_ori_3),-1.*imag(cor_ori_3),'k--','LineWi dth',1);
fill(x5,y5,colours(2,:));
fill(-1.*x5,y5,colours(3,:));
axis square
hold off
xlabel('mm of cortical area'), ylabel('mm of cortical area');
title('Centre and periphery projected to V1');
set(gca,'FontSize','12','FontName','Arial');


%% draw the mask

centre_x = [0; -440; -440; 0; 0;wholefield_x_left; wholefield_x_right;0;0;440;440;0];
centre_y = [-440; -440; 440; 440; 62.30; wholefield_y_left; wholefield_y_right;62.30;440;440;-440;-440];

mask_centre = patch(centre_x,centre_y,[.5 .5 .5],'EdgeColor','none');
axis square;
set(gcf,'Color','none');
set(gca, 'XColor','none','YColor','none')

mask_peripheral = plot([-440.5 -440.5 440.5 440.5 -440.5], [-440 440 440 -440 -440],'-','Color',[.5 .5 .5],'LineWidth',0.5);
hold on;
fill([0;wholefield_x_left; wholefield_x_right;0],[62.30;wholefield_y_left; wholefield_y_right;62.30],[.5 .5 .5],'EdgeColor','none');
hold off
axis square
set(gcf,'Color','none');
set(gca, 'XColor','none','YColor','none');

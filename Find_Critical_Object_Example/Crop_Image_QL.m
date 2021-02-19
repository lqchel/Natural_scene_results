%this code cuts images into 12 patches, among which 4 locate in the center
%of the square image

%load the files and file names
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incongruent'; 
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\congruent';
folder3 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\photo difference'; 
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
filePattern3 = fullfile(folder3, '*.jpg');

theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);
theFiles3 = dir(filePattern3);

% create save locations

mkdir('cong_periphery')
mkdir('incong_periphery')
mkdir('cong_center')
mkdir('incong_center')

 
%% Crop images
%crop patches in location 1-8
edge = 440/3;
edge_s = 440/6;
for a = 1:2

%%%Reading congruent image

filename1 = theFiles1(a).name;
fullname1 = fullfile(folder1,filename1);
fprintf(1,'Now reading %s\n',fullname1');
congruent = imread(fullname1);

%%%Reading incongruent image
filename2 = theFiles2(a).name;%see listing-file attributes in dir
fullname2 = fullfile(folder2,filename2);
fprintf(1,'Now reading %s\n',fullname2);
incongruent = imread(fullname2);

%reading difference image. these images shows the difference between
%congruent and incongruent image, by deducting incongruent image (in RGB matrix) from
%congruent image (in RGB matrix)
%so that in difference images, only the critical objects part are not black
%(therefore sum of this area is not 0)
% for cutting other images e.g. scrambled masks, this part can be ignored

filename3 = theFiles3(a).name;
fullname3 = fullfile(folder3,filename3);
fprintf(1,'Now reading %s\n',fullname3');
difference = imread(fullname3);

% calculate average difference on every pixel location, and determine which
% patch contains the largest area of the critical object
image_index_l = zeros(3,3);

for i = 1:3
    for b = 1:3 
    dividedImage_d = imcrop(difference, [(b-1).*edge (i-1).*edge edge edge]);
        if i == 2 && b == 2
            
            image_index_s = zeros(2,2);
            for s = 1:2
                for d = 1:2
                 patch_temp_d = imcrop(dividedImage_d,[(d-1).*edge_s (s-1).*edge_s edge_s edge_s]);
                 image_index_s(d,s) = sum(sum(sum(patch_temp_d)))/edge_s^2;
                end
                clear patch_temp_d
            end
            
        else
            image_index_l(b,i) = sum(sum(sum(dividedImage_d)))/edge^2;
        end
    clear dividedImage_d    
    end
end
s_max = max(max(image_index_s));
l_max = max(max(image_index_l));
P_index = max([s_max; l_max]);

% then, crop the images into small patches
patch_number_s = [9 10; 11 12];
patch_number_l = [1 2 3; 4 0 5; 6 7 8];

for i = 1:3
    for b = 1:3 
        dividedImage_I = imcrop(incongruent,[(b-1).*edge (i-1).*edge edge edge]);
        dividedImage_C = imcrop(congruent, [(b-1).*edge (i-1).*edge edge edge]);
 
        if i == 2 && b == 2           
            for s = 1:2
                for d = 1:2
                 dividedImage_Is = imcrop(dividedImage_I,[(d-1).*edge_s (s-1).*edge_s edge_s edge_s]);
                 dividedImage_Cs = imcrop(dividedImage_C,[(d-1).*edge_s (s-1).*edge_s edge_s edge_s]);
                 
                 if image_index_s(d,s)== P_index
                     filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_center',...
                         ['cong_',num2str(a),'_', num2str(patch_number_s(d,s)),'_p', '.jpg']);
                     filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_center',...
                         ['incong_',num2str(a),'_', num2str(patch_number_s(d,s)),'_p', '.jpg']);
                 else
                     filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_center',...
                         ['cong_',num2str(a),'_', num2str(patch_number_s(d,s)),'_a', '.jpg']);
                     filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_center',...
                         ['incong_',num2str(a),'_', num2str(patch_number_s(d,s)),'_a', '.jpg']);

                 end
                 imwrite(dividedImage_Is,filename_I);
                 imwrite(dividedImage_Cs,filename_C);
                clear dividedImage_Is
                clear dividedImage_Cs
                clear filename_C
                clear filename_I
                end
            end 
            
        else
              if image_index_l(b,i)== P_index
                 filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_periphery',...
                     ['cong_',num2str(a),'_', num2str(patch_number_l(b,i)),'_p','.jpg']);
                 filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_periphery',...
                     ['incong_',num2str(a),'_', num2str(patch_number_l(b,i)),'_p', '.jpg']);
              else
                 filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_periphery',...
                     ['cong_',num2str(a),'_', num2str(patch_number_l(b,i)),'_a','.jpg']);
                 filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_periphery',...
                     ['incong_',num2str(a),'_', num2str(patch_number_l(b,i)),'_a', '.jpg']);
              end
             imwrite(dividedImage_I,filename_I);
             imwrite(dividedImage_C,filename_C);
        end
        clear dividedImage_I    
        clear dividedImage_C
        clear filename_C
        clear filename_I
    end
end


end
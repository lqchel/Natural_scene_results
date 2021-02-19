%this code cuts images into 12 patches, among which 4 locate in the center
%of the square image

%load the files and file names
folder2 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incongruent'; 
folder1 = 'C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\congruent';

%%% full resource of the difference image, go to
%%% "Qianchen_Liad_Natural_Scene\experiment codes\squarephotodifferences"
%%% on google share drive
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

for a = 1:length(theFiles1)

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
center_or_periphery = zeros(length(theFiles1),1);

for i = 1:3
    for b = 1:3 
    dividedImage_d = imcrop(difference, [(b-1).*edge (i-1).*edge edge edge]);
    image_index_l(i,b) = sum(sum(sum(dividedImage_d)));
    end
end
disp(image_index_l)

%%% if the central location (consist of location 9, 10, 11, 12 altogether),
%%% takes the largest critical object area, then find the largest object 
%%% area among the four small patches in the center; if not, find the
%%% largest object area in peripheral patch

[row, col] = find(image_index_l== max(max(image_index_l)));
if row == 2 && col == 2
    center_or_periphery(a,1) = 1;
    image_index_s = zeros(2,2);
    middle_patch = imcrop(differece,[edge edge edge edge]);
    for s = 1:2
        for d = 1:2
         patch_temp_d = imcrop(middle_patch,[(d-1).*edge_s (s-1).*edge_s edge_s edge_s]);
         image_index_s(s,d) = sum(sum(sum(patch_temp_d))); % finds the largest object area among the four small patches in the center
        end
    end
    P_index = max(max(image_index_s));            
else
    P_index = max(max(image_index_l));
    center_or_periphery(a,1) = 0;
end
   
    
% then, crop the images into small patches
patch_number_s = [9 10; 11 12];
patch_number_l = [1 2 3; 4 0 5; 6 7 8];

for i = 1:3
    for b = 1:3 
        if a == 1
        disp([(b-1).*edge (i-1).*edge edge edge])
        end
        dividedImage_I = imcrop(incongruent,[(b-1).*edge (i-1).*edge edge edge]);
        dividedImage_C = imcrop(congruent, [(b-1).*edge (i-1).*edge edge edge]);
 
        if i == 2 && b == 2           
            for s = 1:2
                for d = 1:2
                 dividedImage_Is = imcrop(dividedImage_I,[(d-1).*edge_s (s-1).*edge_s edge_s edge_s]);
                 dividedImage_Cs = imcrop(dividedImage_C,[(d-1).*edge_s (s-1).*edge_s edge_s edge_s]);
                 
                 % filename of patch with largest critical object area ends
                 % with "p", and filenames of other patches end with "a" 
                 if image_index_s(s,d)== P_index && center_or_periphery(a,1) == 1
                     filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_center',...
                         ['cong_',num2str(a),'_', num2str(patch_number_s(s,d)),'_p', '.jpg']);
                     filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_center',...
                         ['incong_',num2str(a),'_', num2str(patch_number_s(s,d)),'_p', '.jpg']);
                 else
                     filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_center',...
                         ['cong_',num2str(a),'_', num2str(patch_number_s(s,d)),'_a', '.jpg']);
                     filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_center',...
                         ['incong_',num2str(a),'_', num2str(patch_number_s(s,d)),'_a', '.jpg']);

                 end
                 imwrite(dividedImage_Is,filename_I);
                 imwrite(dividedImage_Cs,filename_C);
                end
            end 
            
        else
                 % filename of patch with largest critical object area ends
                 % with "p", and filenames of other patches end with "a" 
              if image_index_l(i,b)== P_index && center_or_periphery(a,1) == 0
                 filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_periphery',...
                     ['cong_',num2str(a),'_', num2str(patch_number_l(i,b)),'_p','.jpg']);
                 filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_periphery',...
                     ['incong_',num2str(a),'_', num2str(patch_number_l(i,b)),'_p', '.jpg']);
              else
                 filename_C = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\cong_periphery',...
                     ['cong_',num2str(a),'_', num2str(patch_number_l(i,b)),'_a','.jpg']);
                 filename_I = fullfile('C:\Users\liang\OneDrive\Documents\honours\research project\Natural_scene_results\Find_Critical_Object_Example\incong_periphery',...
                     ['incong_',num2str(a),'_', num2str(patch_number_l(i,b)),'_a', '.jpg']);
              end
             imwrite(dividedImage_I,filename_I);
             imwrite(dividedImage_C,filename_C);
        end
    end
end


end
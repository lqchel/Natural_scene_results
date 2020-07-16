%this code cuts the congruent/incongruent images into patches, find out
%the patches that contains the largest proportion of the critical object,
%thereby determining the presence/absence of incongruent/congruent object and name them
%accordingly


clear all;
%%% load the files and file names
folder2 = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\squareimage\incongruent cropped'; 
folder1 = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\squareimage\congruent cropped';
folder3 = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\squarephotodifference'; 
filePattern1 = fullfile(folder1,'*.jpg');
filePattern2 = fullfile(folder2,'*.jpg');
filePattern3 = fullfile(folder3, '*.jpg');
theFiles1 = dir(filePattern1);
theFiles2 = dir(filePattern2);
theFiles3 = dir(filePattern3);

%%% access and crop each of the congruent/incongruent image
for a = 1:length(theFiles1)
    
%Reading incongruent image
filename2 = theFiles2(a).name;
fullname2 = fullfile(folder2,filename2);
fprintf(1,'Now reading %s\n',fullname2);
incongruent = imread(fullname2);

%Reading congruent image

filename1 = theFiles1(a).name;
fullname1 = fullfile(folder1,filename1);
fprintf(1,'Now reading %s\n',fullname1');
congruent = imread(fullname1);

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
edge = 440/3; %defining the edge and therefore size of the patch

%cut each difference image into 9 equal
%patches, calculate the sum of the difference within each patch, and mark down the number of the patch that contains the largest difference 

    for i = 1:3

        for b = 1:3 
        dividedImage_d = imcrop(difference, [(b-1).*edge (i-1).*edge edge edge]);
        image_index(b,i) = sum(sum(sum(dividedImage_d)));
        end
    end
    P_index = max(max(image_index));
 %cut congruent and incongruent images into patches; name the patches corresponding to the maximal difference patch as "..._p",
 %meaning that critical objecct is present; and name others as "..._a"
 %I'm having trouble here, as many patches that are detected as
 %object-present in fact do not contains the object at all
 
    for i = 1:3
        for b = 1:3
            dividedImage_I = imcrop(incongruent,[(b-1).*edge (i-1).*edge edge edge]);
            dividedImage_C = imcrop(congruent, [(b-1).*edge (i-1).*edge edge edge]);
            
      if image_index(b,i)== P_index
         imwrite(dividedImage_I,strcat('incong_',num2str(a),'_', num2str((i-1).*3 + b),'_p', '.jpg'))
         imwrite(dividedImage_C,strcat('cong_',num2str(a),'_', num2str((i-1).*3 + b),'_p','.jpg'))
         
     else
         imwrite(dividedImage_I,strcat('incong_',num2str(a),'_', num2str((i-1).*3 + b),'_a', '.jpg'))
         imwrite(dividedImage_C,strcat('cong_',num2str(a),'_', num2str((i-1).*3 + b),'_a','.jpg'))
     end
        end
    end
end

    

  
 
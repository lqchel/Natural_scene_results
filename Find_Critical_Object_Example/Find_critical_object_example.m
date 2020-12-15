% load image
cong = imread('SquareCongruent_002.jpg');
incong = imread('SquareIncongruent_002.jpg');

difference = abs(cong-incong); % compute difference matrix between congruent and incongruent images
[row,column] = find(sum(difference,3)>0); % find the area of difference, by returning the row and column indices of pixels within that area 
object_row_coordinate = min(row):max(row); % find min and max of the row and column indices, which defines the size and location of the rectangle cue
object_col_coordinate = min(column):max(column);

% now, we can find and cut out the area for critical objects in both
% congruent and incongruent images; you can utilize the coordinates for the
% cues
imwrite(cong(object_row_coordinate,object_col_coordinate,:),'cong_object_002.jpg');
imwrite(incong(object_row_coordinate,object_col_coordinate,:),'incong_object_002.jpg');
% this is the instruction of use for my data

% there are 13 columns in the data matrix:
% 1st column = participant number, 2nd = image/trial number (starting from 3, as the first two images are for pratice trials)
% 3rd = patch number within a trial, 4th --> congruent image = 0, incongruent = 1,
% 5th = patch source (1 = Nishimoto's image set, 2 = congruent image, 3 = incongruent image)
% 6th = critical object presence (1 = randomly selected from Nishimoto's image, so no critical object,
% 2 = no critical object, 3 = critical object present)
% 7th = presentation location, 8th = yes(1)/no(-1) judgment
% 9th = confidence (1-4), 10th = reaction time in seconds
% 13th = eccentricity in dva

%import data file
Results = importdata('Pooled Results.mat');
Results(:,9) = Results(:,8).*Results(:,9); % combined yes/no with confidence levels; if necessary, Results(:,9)-0.5 before combining


%% patch indexes -- this helps you to find specific patch responses / rows 

%Find Null patches(patches randomly selected from pther images) -- signal absent for hypo 1
Find_N = Results(:,5) ==1;

% present patch trials -- signal present for hypo 1
Find_CAP = Results(:,4) == 0 & Results(:,5) == 2; % congruent image presented, patch present
Find_IAP = Results(:,4) == 1 & Results(:,5) == 3;  % incongruent image present, patch present

%Congruent trial (congruent image presented) with Congruent object, and incongruent trial with
%incongruent object -- signal present for hypo 2
Find_Congruent_CP = Results(:,4) == 0 & Results(:,5) == 2 & Results(:,6) == 3;
Find_Incongruent_IP = Results(:,4) == 1 & Results(:,5) == 3 & Results(:,6) == 3;

%Incongruent trial with congruent object, congruent trial with incongruent
%object -- signal absent for hypo 2

Find_Congruent_IP = Results(:,4) == 0 & Results(:,5) == 3 & Results(:,6) == 3;
Find_Incongruent_CP = Results(:,4) == 1 & Results(:,5) == 2 & Results(:,6) == 3;



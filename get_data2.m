%%% this code loads the Exp 2 data and trial indexes

function [Results] = get_data2(Exp)
if Exp == 1
    Results = importdata('Pooled Results.mat');
    Results(:,6) = Results(:,6) == 3;
    Results(:,9) = Results(:,8).*Results(:,9);
    location1 = Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8;
    location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 2;
    eccentricity_level = zeros(length(Results),1)+ location1 + location2;
    real_location1 = (Results(:,7)==2 | Results(:,7)==4| Results(:,7)==6|Results(:,7)== 8) .* 6.5;
    real_location2 = (Results(:,7)==1 | Results(:,7)==3| Results(:,7)==7|Results(:,7)== 9) .* 9.2;
    actual_eccentrcity = zeros(length(Results),1) + real_location1 + real_location2;
    Results = [Results eccentricity_level actual_eccentrcity];
    Results(:,11) = Results(:,2);
elseif Exp == 2
    Results = importdata('Exp2_data.mat');
    Results = Results(Results(:,11)~=0,:);
    Results(:,9) = Results(:,8).*Results(:,9);
end
end


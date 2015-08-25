function [fitness,threshold] = calcVolumetricMatchEdgeExclusion(targetStructure,exclusionStructure,edgeExclusionStructure, exposedStruct)
%Calculates the amount of overlap between the target structure and
%simulated structure
%Note: target structure and simulated structure are inverted: 1=void, 0=SU8
%Sets fill factor of simulated structure equal to that of the target

% if (max(max(max(intensityDist))) - min(min(min(intensityDist))) ) / min(min(min(intensityDist))) < 0.01 %if approx no variation in intensity
%     fitness = 0;
%     threshold =  0;
%     return 
% end

%nbins = 256; %number of bins for intensity histogram

%Need fill of target structure to normalize
fill = 1 - sum(sum(sum(targetStructure))) / (size(targetStructure,1)*size(targetStructure,1)*size(targetStructure,3)); %gives fill fraction of SU8 of target
%threshold = fixfill(reshape(intensityDist,1,[]),nbins,fill); %Calculates the threshold value that will yield desired fill fraction

%simStruct = intensityDist<threshold;  %1 if intensity below threshold (void), 0 if above (SU8)
simStruct = ~exposedStruct; 

%f_void = sum(sum(sum( targetStructure & simStruct )))/(size(targetStructure,1)*size(targetStructure,1)*size(targetStructure,3))/fill; %Fraction of points that are contained by both
%f_filled = sum(sum(sum( ~targetStructure & ~simStruct )))/(size(targetStructure,1)*size(targetStructure,1)*size(targetStructure,3))/fill; %Fraction of poitns contained by neither
%Add -1 of 1@1 for target structure ( 1 = in helix)  %Negative fitness is good
f1 = -1*sum(sum(sum( targetStructure & simStruct )))/(size(targetStructure,1)*size(targetStructure,1)*size(targetStructure,3))/fill; %Fraction of poitns contained by neither
%Add 0@1 for exclusion structure (1 = in helix)  %If there are points outside the big helix, add to the fitness (bad)
f2 = sum(sum(sum( ~exclusionStructure & simStruct )))/(size(targetStructure,1)*size(targetStructure,1)*size(targetStructure,3))/fill; %Fraction of poitns contained by neither
%Add -1 of 1@1 (1 is at edge) %If there are points at the edge, add to the
%fitness % normalized to ~total number of edge points
f3 = sum(sum(sum( edgeExclusionStructure & simStruct )))/(size(targetStructure,1)*size(targetStructure,3)*2+size(targetStructure,2)*size(targetStructure,3)*2);

fitness = f1*4 + f2 + f3;

% %TEST
% figure
% patched = patch(isosurface(padarray(simStruct,[1,1,1]),0.5));
% set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
% view(3);
%     camlight
% 
%     lighting gouraud
%     xlim([1 size(intensityDist,1)])
%     ylim([1 size(intensityDist,2)])
%     zlim([1 size(intensityDist,3)])


end
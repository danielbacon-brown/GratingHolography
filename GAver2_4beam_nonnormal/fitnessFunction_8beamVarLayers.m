function fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome,doPlots)
%Designed to be a fairly general fitness function
%Use in a separate file allows the chromosome to be measured directly

%Get layer chromosome length
chromNlayer=0;
for i=1:length(GAoptions.S4interfaceOptions.layers)
   chromNlayer=chromNlayer+GAoptions.S4interfaceOptions.layers(i).getChromosomeSize();
end

%Splits chromosome according to each module
[gratingChromosome ,incidentLightChromosome, offsetChromosome,layerChromosome ] = splitChromosome(chromosome,[ ...
    GAoptions.gratingFunction.getChromosomeSize(), ...
    GAoptions.incidentLightFunction.getChromosomeSize(),  ...
    ...GAoptions.calcStructureFunction.getChromosomeSize() ...
    GAoptions.offsetConductor.getChromosomeSize(), ...
    chromNlayer ...
    ...sum(GAoptions.S4interfaceOptions.layers(:).getChromosomeSize()) ...
    ]);

%Splits layer chromosome into 1 for each layer
layerChromosomes = cell(length(GAoptions.S4interfaceOptions.layers),1);
u=0;
for i_l = 1:length(GAoptions.S4interfaceOptions.layers)  %IF the layer has constant thickness, don't make a chromosome for it
    if GAoptions.S4interfaceOptions.layers(i_l).constThick >= 0
        layerChromosomes{i_l} = [];
    else
        layerChromosomes{i_l} = layerChromosome( (u+1):(u+GAoptions.S4interfaceOptions.layers(i_l).getChromosomeSize()) );
        u = u + GAoptions.S4interfaceOptions.layers(i_l).getChromosomeSize();
    end
end


%Generate Grating
grating = GAoptions.gratingFunction.generateGrating(gratingChromosome);

%TEST!!!!!!!!!!!!!!!!!
%GAoptions.fill = 0.85
%grating = createTestGrating();


if exist('doPlots','var')
    GAoptions.gratingFunction.plotGrating(grating);  %plot grating
end



%Incident Field
incidentFieldParams = GAoptions.incidentLightFunction.generateField(incidentLightChromosome);

if exist('doPlots','var')
   GAoptions.incidentLightFunction.plotPolarizationEsp(incidentFieldParams.Esp); 
end

%TEST!!!!!!!!!
%incidentEsp = [0;1];


% %Do RCWA analysis with GD-Calc:
% %[~, scat_field, ~] = gdc(grating,incidentField,GAoptions.gratingOptions.order);
% 
% %Calculate propagating diffracted beams
% %[E231, k231] = field_convert(incidentEsp,scat_field);  %DO OWN REWRITE
% 
% %Calculate intensity distribution:
% %intensityDist = GAoptions.calcStructureFunction.calcIntensity( E231,k231,GAoptions.dimensions,GAoptions.cells,offsetChromosome);

%Do RCWA analysis with S4:
intensityDist = GAoptions.S4interface.doRCWA(GAoptions,grating,incidentFieldParams,layerChromosomes);  %Note: this intensityDist is based on a hexagonal parallelogram
if length(intensityDist)<1  %If the RCWA failed, return bad fitness and move on
   fitness = 10;
   return;
end



%Do offset of interference pattern
intensityDist = GAoptions.offsetConductor.doOffset(intensityDist, offsetChromosome); %W/(unitarea)


%Calculate Fitness
[fitness,threshold] = calcVolumetricMatchEdgeExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure,GAoptions.edgeExclusionStructure, intensityDist,GAoptions.fill);
%threshold = fixfill(reshape(intensityDist,1,[]),256,GAoptions.fill); %Calculates the threshold value that will yield desired fill fraction
%fitness = rand();



if exist('doPlots','var') %Plot simulated structure
    figure
    patched = patch(isosurface(padarray(intensityDist,[1,1,1],100),threshold));
    set(patched,'FaceColor', [30 255 30]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([1 size(intensityDist,1)+2])
    ylim([1 size(intensityDist,2)+2])
    zlim([1 size(intensityDist,3)+2])
    
    
%     %Do Lumerical simulation
%     writeLumericalRunFileSquare(GAoptions, intensityDist>threshold);
%     system(['fdtd-solutions -run ', GAoptions.dir, GAoptions.LumRunScript]);
%     while(~exist([GAoptions.dir,GAoptions.currentLumResultsFile,'.mat'],'file'))
%        pause(0.1)
%     end
%     LumResults = load([GAoptions.dir,GAoptions.currentLumResultsFile]);
%     transmissionRight = LumResults.transmission_right/sqrt(2);
%     transmissionLeft = LumResults.transmission_left/sqrt(2);
%     reflectionRight = LumResults.reflection_right*-1/sqrt(2);
%     reflectionLeft = LumResults.reflection_left*-1/sqrt(2);
%     
%     %Plot transmission and reflection of RCP and LCP waves
%     figure
%     frequencies = linspace(2.99e8/GAoptions.fdtd.maxMeasWL, 2.99e8/GAoptions.fdtd.minMeasWL, GAoptions.fdtd.numMeasWL); %Linear in frequency space
%     %wavelengths = linspace(GAoptions.fdtd.minMeasWL,GAoptions.fdtd.maxMeasWL,GAoptions.fdtd.numMeasWL);
%     wavelengths = 2.99e8./frequencies; %Convert to wavelength
%     plot(wavelengths,transmissionRight,'r',wavelengths,transmissionLeft,'b');
%     figure
%     plot(wavelengths,reflectionRight,'r',wavelengths,reflectionLeft,'b');




    
    
    %Do sensitizer sim
    acidCount = excitePAG(intensityDist,GAoptions.dimensions,GAoptions.sensSim.sensDens,GAoptions.sensSim.absCrossSection,GAoptions.sensSim.Texposure);
    acidMax = max(max(max(acidCount)))
    acidMin = min(min(min(acidCount)))
    threshold = fixfill(acidCount,256,GAoptions.fill); %Calculates the threshold value that will yield desired fill fraction
    figure
    patched = patch(isosurface(padarray(acidCount,[1,1,1],100),threshold));
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([1 size(acidCount,1)+2])
    ylim([1 size(acidCount,2)+2])
    zlim([1 size(acidCount,3)+2])
    
    
    %DO SMOOTHING//DIFFUSION (WITH SUMMING)
    diffRate = 5;
    smoothed = zeros([size(acidCount),8]);
    sizesmoothedout= size(smoothed)
    smoothed(:,:,:,1) = smooth3(acidCount,'gaussian',diffRate);
    for smooth_i = 2:size(smoothed,4)
        smoothed(:,:,:,smooth_i) = smooth3(smoothed(:,:,:,smooth_i-1),'gaussian',diffRate);
        disp(smooth_i)
    end
    crosslinkCount = sum(smoothed,4);
    sizesmoothedout= size(crosslinkCount)
    crosslinkThreshold = fixfill(crosslinkCount, 256, GAoptions.fill)
    figure
    patched = patch(isosurface(padarray(crosslinkCount,[1,1,1],100),crosslinkThreshold));
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([1 size(acidCount,1)+2])
    ylim([1 size(acidCount,2)+2])
    zlim([1 size(acidCount,3)+2])
    
    
    %CALCULATE REPEAT CELLS
    %GAoptions.calcStructureFunction.cacheSet = false;
    %GAoptions.calcStructureFunction.prepareCommon(k231,GAoptions.dimensions*2,GAoptions.cells*2);
    %intensityDistRepeat = GAoptions.calcStructureFunction.calcIntensity( E231,k231,GAoptions.dimensions*2,GAoptions.cells,offsetChromosome);
    
%     figure
%     patched = patch(isosurface(padarray(intensityDistRepeat,[1,1,1],100),threshold));
%     set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
%     view(3);
%     camlight
%     axis equal
%     lighting gouraud
%     xlim([1 size(intensityDistRepeat,1)+2])
%     ylim([1 size(intensityDistRepeat,2)+2])
%     zlim([1 size(intensityDistRepeat,3)+2])
    
    
    %Save grating and incident field for S4 simulation
%    save('lastGratingandField','grating','incidentEsp')
    
%    doS4(grating,incidentEsp)
    
end
        
end
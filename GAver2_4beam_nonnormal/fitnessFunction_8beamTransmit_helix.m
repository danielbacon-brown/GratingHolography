function fitness = fitnessFunction_8beamTransmit_helix(GAoptions,chromosome,doPlots)
%Designed to be a fairly general fitness function
%Use in a separate file allows the chromosome to be measured directly



%Splits chromosome according to each module
[gratingChromosome ,incidentLightChromosome, offsetChromosome ] = splitChromosome(chromosome,[ ...
    GAoptions.gratingFunction.getChromosomeSize(), ...
    GAoptions.incidentLightFunction.getChromosomeSize(),  ...
    ...GAoptions.calcStructureFunction.getChromosomeSize() ...
    GAoptions.offsetConductor.getChromosomeSize() ...
    ]);

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
intensityDist = GAoptions.S4interface.doRCWA(GAoptions,grating,incidentFieldParams);  %Note: this intensityDist is based on a hexagonal parallelogram




%Do offset of interference pattern
intensityDist = GAoptions.offsetConductor.doOffset(intensityDist, offsetChromosome);


%Calculate Fitness
[fitness,threshold] = calcVolumetricMatchExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure, intensityDist,GAoptions.fill);
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
    
    
    writeLumericalRunFileSquare(GAoptions, intensityDist>threshold);
    system(['fdtd-solutions -run ', GAoptions.dir, GAoptions.LumRunScript]);
    
    
    while(~exist([GAoptions.dir,GAoptions.currentLumResultsFile,'.mat'],'file'))
        pause(0.1)
    end
        
    LumResults = load([GAoptions.dir,GAoptions.currentLumResultsFile]);
    transmissionRight = LumResults.transmission_right/sqrt(2);
    transmissionLeft = LumResults.transmission_left/sqrt(2);
    reflectionRight = LumResults.reflection_right*-1/sqrt(2);
    reflectionLeft = LumResults.reflection_left*-1/sqrt(2);
    
    
    figure
    frequencies = linspace(2.99e8/GAoptions.fdtd.maxMeasWL, 2.99e8/GAoptions.fdtd.minMeasWL, GAoptions.fdtd.numMeasWL); %Linear in frequency space
    %wavelengths = linspace(GAoptions.fdtd.minMeasWL,GAoptions.fdtd.maxMeasWL,GAoptions.fdtd.numMeasWL);
    wavelengths = 2.99e8./frequencies; %Convert to wavelength
    plot(wavelengths,transmissionRight,'r',wavelengths,transmissionLeft,'b');
    
    figure
    plot(wavelengths,reflectionRight,'r',wavelengths,reflectionLeft,'b');
    
    
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
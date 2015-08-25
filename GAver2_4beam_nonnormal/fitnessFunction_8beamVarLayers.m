function fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome,doPlots)
%Designed to be a fairly general fitness function
%Use in a separate file allows the chromosome to be measured directly

%Get layer chromosome length
chromNlayer=0;
for i=1:length(GAoptions.S4interfaceOptions.layers)
   chromNlayer=chromNlayer+GAoptions.S4interfaceOptions.layers(i).getChromosomeSize();
end
%Get material chromosome length
chromNmaterial=0;
for i=1:length(GAoptions.S4interfaceOptions.materials)
   chromNmaterial=chromNmaterial+GAoptions.S4interfaceOptions.materials(i).getChromosomeSize();
end

%Splits chromosome according to each module
[gratingChromosome ,incidentLightChromosome, offsetChromosome,layerChromosome,materialChromosome ] = splitChromosome(chromosome,[ ...
    GAoptions.gratingFunction.getChromosomeSize(), ...
    GAoptions.incidentLightFunction.getChromosomeSize(),  ...
    GAoptions.offsetConductor.getChromosomeSize(), ...
    chromNlayer ...
    chromNmaterial ...
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

%Splits layer chromosome into 1 for each layer
materialChromosomes = cell(length(GAoptions.S4interfaceOptions.materials),1);
u=0;
for i_m = 1:length(GAoptions.S4interfaceOptions.materials)  %IF the layer has constant thickness, don't make a chromosome for it
    if GAoptions.S4interfaceOptions.materials(i_m).constRI >= 0
        materialChromosomes{i_m} = [];
    else
        materialChromosomes{i_m} = materialChromosome( (u+1):(u+GAoptions.S4interfaceOptions.materials(i_m).getChromosomeSize()) )
        u = u + GAoptions.S4interfaceOptions.materials(i_m).getChromosomeSize();
    end
end


%Generate Grating
grating = GAoptions.gratingFunction.generateGrating(gratingChromosome);

if exist('doPlots','var')
    GAoptions.gratingFunction.plotGrating(grating);  %plot grating
end



%Incident Field
incidentFieldParams = GAoptions.incidentLightFunction.generateField(incidentLightChromosome);

if exist('doPlots','var')
   GAoptions.incidentLightFunction.plotPolarizationEsp(incidentFieldParams.Esp); 
end


%Do RCWA analysis with S4:
intensityDist = GAoptions.S4interface.doRCWA(GAoptions,grating,incidentFieldParams,layerChromosomes,materialChromosomes);  %Note: this intensityDist is based on a hexagonal parallelogram
if length(intensityDist)<1  %If the RCWA failed, return bad fitness and move on
   fitness = 10;
   return;
end



%Do offset of interference pattern
intensityDist = GAoptions.offsetConductor.doOffset(intensityDist, offsetChromosome); %W/(unitarea)


%Calculate Fitness
%[fitness,threshold] = calcVolumetricMatchEdgeExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure,GAoptions.edgeExclusionStructure, intensityDist,GAoptions.fill);
[fitness,threshold] = calcVolumetricMatchExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure, intensityDist,GAoptions.fill);
%threshold = fixfill(reshape(intensityDist,1,[]),256,GAoptions.fill); %Calculates the threshold value that will yield desired fill fraction
%fitness = rand();



if exist('doPlots','var') %Plot simulated structure
    figure
    patched = patch(isosurface(padarray(intensityDist,[1,1,1],1e20),threshold));
    set(patched,'FaceColor', [30 255 30]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([2 size(intensityDist,2)+1])
    ylim([2 size(intensityDist,1)+1])
    zlim([2 size(intensityDist,3)+1])
    
    
    if strcmp(GAoptions.lattice, 'square') && ~strcmp(strtrim(GAoptions.hostname),'Daniel-netbook')
        
        
        %Do Lumerical simulation
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
        
        %Plot transmission and reflection of RCP and LCP waves
        figure
        frequencies = linspace(2.99e8/GAoptions.fdtd.maxMeasWL, 2.99e8/GAoptions.fdtd.minMeasWL, GAoptions.fdtd.numMeasWL); %Linear in frequency space
        %wavelengths = linspace(GAoptions.fdtd.minMeasWL,GAoptions.fdtd.maxMeasWL,GAoptions.fdtd.numMeasWL);
        wavelengths = 2.99e8./frequencies; %Convert to wavelength
        plot(wavelengths,transmissionRight,'r',wavelengths,transmissionLeft,'b');
        figure
        plot(wavelengths,reflectionRight,'r',wavelengths,reflectionLeft,'b');
        

    elseif strcmp(GAoptions.lattice, 'hexagonal')
        %If hexagonal, convert to cartesian system for plotting and lumerical:
        intensityDist = hex2cart(intensityDist, GAoptions.cellsCart);
        
        %Plot cartesian
        figure
        patched = patch(isosurface(padarray(intensityDist,[1,1,1],1e20),threshold));
        set(patched,'FaceColor', [30 255 30]/256, 'EdgeColor', 'none');
        view(3);
        camlight
        axis equal
        lighting gouraud
        xlim([2 size(intensityDist,2)+1])
        ylim([2 size(intensityDist,1)+1])
        zlim([2 size(intensityDist,3)+1])
        
        if ~strcmp(strtrim(GAoptions.hostname),'Daniel-netbook')
            %Do Lumerical simulation
            
            %Make hexagonal repeat of structure:
            lx = size(intensityDist,1);
            left = intensityDist(1:floor(lx/2),:,:);
            right = intensityDist((floor(lx/2)+1):end,:,:);
            top = cat(1,right,left);
            intensityDistComb = cat(2,intensityDist,top);
            
            writeLumericalRunFileSquare(GAoptions, intensityDistComb>threshold);
            system(['fdtd-solutions -run ', GAoptions.dir, GAoptions.LumRunScript]);
            while(~exist([GAoptions.dir,GAoptions.currentLumResultsFile,'.mat'],'file'))
                pause(0.1)
            end
            LumResults = load([GAoptions.dir,GAoptions.currentLumResultsFile]);
            transmissionRight = LumResults.transmission_right/sqrt(2);
            transmissionLeft = LumResults.transmission_left/sqrt(2);
            reflectionRight = LumResults.reflection_right*-1/sqrt(2);
            reflectionLeft = LumResults.reflection_left*-1/sqrt(2);
            
            %Plot transmission and reflection of RCP and LCP waves
            figure
            frequencies = linspace(2.99e8/GAoptions.fdtd.maxMeasWL, 2.99e8/GAoptions.fdtd.minMeasWL, GAoptions.fdtd.numMeasWL); %Linear in frequency space
            %wavelengths = linspace(GAoptions.fdtd.minMeasWL,GAoptions.fdtd.maxMeasWL,GAoptions.fdtd.numMeasWL);
            wavelengths = 2.99e8./frequencies; %Convert to wavelength
            plot(wavelengths,transmissionRight,'r',wavelengths,transmissionLeft,'b');
            figure
            plot(wavelengths,reflectionRight,'r',wavelengths,reflectionLeft,'b');
        
        end
        
    end
    
    
%     
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
% 



    
    
    %Do sensitizer sim
    acidCount = excitePAG(intensityDist,GAoptions.dimensions,GAoptions.sensSim.sensDens,GAoptions.sensSim.absCrossSection,GAoptions.sensSim.Texposure);
    acidMax = max(max(max(acidCount)))
    acidMin = min(min(min(acidCount)))
    threshold = fixfill(acidCount,256,GAoptions.fill); %Calculates the threshold value that will yield desired fill fraction
    figure
    patched = patch(isosurface(padarray(acidCount,[1,1,1],1e20),threshold));
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([2 size(acidCount,2)+1])
    ylim([2 size(acidCount,1)+1])
    zlim([2 size(acidCount,3)+1])
    
    
    %DO SMOOTHING//DIFFUSION (WITH SUMMING)
    diffRate = 5;
    smoothed = zeros([size(acidCount),8]);
    smoothed(:,:,:,1) = smooth3(acidCount,'gaussian',diffRate);
    for smooth_i = 2:size(smoothed,4)
        smoothed(:,:,:,smooth_i) = smooth3(smoothed(:,:,:,smooth_i-1),'gaussian',diffRate);
        disp(smooth_i)
    end
    crosslinkCount = sum(smoothed,4);
    crosslinkThreshold = fixfill(crosslinkCount, 256, GAoptions.fill)
    figure
    patched = patch(isosurface(padarray(crosslinkCount,[1,1,1],1e20),crosslinkThreshold));
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([2 size(acidCount,2)+1])
    ylim([2 size(acidCount,1)+1])
    zlim([2 size(acidCount,3)+1])
    
    
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
    
    
    
end
        
end
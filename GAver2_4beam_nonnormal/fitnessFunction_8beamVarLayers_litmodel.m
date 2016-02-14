function fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome)
%Designed to be a fairly general fitness function
%Use in a separate file allows the chromosome to be measured directly

plotHeatMaps = 0;
twoPhoton = 1;
colormap(hot)


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
[gratingChromosome ,incidentLightChromosome, offsetChromosome,layerChromosome,materialChromosome, fillChromosome ] = splitChromosome(chromosome,[ ...
    GAoptions.gratingFunction.getChromosomeSize(), ...
    GAoptions.incidentLightFunction.getChromosomeSize(),  ...
    GAoptions.offsetConductor.getChromosomeSize(), ...
    chromNlayer ...
    chromNmaterial ...
    GAoptions.fillHandler.getChromosomeSize() ...
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


    %plotS4grating(grating,GAoptions.lattice)
if GAoptions.runSingle == 1
    %GAoptions.gratingFunction.plotGrating(grating);  %plot grating
%    plotS4grating(grating,GAoptions.lattice)
end



%Incident Field
incidentFieldParams = GAoptions.incidentLightFunction.generateField(incidentLightChromosome);

if GAoptions.runSingle == 1
   GAoptions.incidentLightFunction.plotPolarizationEsp(incidentFieldParams.Esp); 
end


%Do RCWA analysis with S4:
intensityDist = GAoptions.S4interface.doRCWA(GAoptions,grating,incidentFieldParams,layerChromosomes,materialChromosomes);  %Note: this intensityDist is based on a hexagonal parallelogram
%Intensitydist is in units of (W/m^2)
if length(intensityDist)<1  %If the RCWA failed, return bad fitness and move on
   fitness = 10;
   return;
end



%Do offset of interference pattern
intensityDist = GAoptions.offsetConductor.doOffset(intensityDist, offsetChromosome); %W/(unitarea)
%Check for about 0 intensity difference
if (max(max(max(intensityDist))) - min(min(min(intensityDist))) ) / min(min(min(intensityDist))) < 0.01 %if approx no variation in intensity
    fitness = 0;
    return 
end


%Calculate Fitness
if strcmp(GAoptions.fitnessType, 'structure')
    
    %Calc structure
    [exposedStruct,fillfrac,threshold] = GAoptions.fillHandler.applyFill(fillChromosome,intensityDist);
    
    thresholdout = threshold
    
    if GAoptions.useExclusion == 1
        %[fitness,threshold] = calcVolumetricMatchEdgeExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure,GAoptions.edgeExclusionStructure, intensityDist,GAoptions.fill);
        if GAoptions.useEdgeExclusion
            fitness = calcVolumetricMatchEdgeExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure,GAoptions.edgeExclusionStructure, exposedStruct);
        else
            fitness = calcVolumetricMatchExclusion(GAoptions.targetStructure, GAoptions.exclusionStructure, exposedStruct);
        end
    else
       fitness = calcVolumetricMatch(GAoptions.targetStructure, exposedStruct); 
    end
    
elseif strcmp(GAoptions.fitnessType, 'fdtd')
    
    %If hexagonal, need to make it repeating hexagonal for it to fit into
    %Lumerical periodicity.
    if strcmp(GAoptions.lattice, 'hexagonal')
        
        %PLOT HEATMAP:
        if plotHeatMaps == 1
            figure
            colormap(hot)
            imagesc( intensityDist(:,:,20), [0,max(max(intensityDist(:,:,20)))] )
            colorbar
        end

    
        %If hexagonal, convert to cartesian system for plotting and lumerical:
        intensityDistInt = hex2cart(intensityDist, GAoptions.cellsCart);
        
        
        %PLOT HEATMAP:
        if plotHeatMaps == 1
            figure
            colormap(hot)
            imagesc(intensityDistInt(:,:,20), [0,max(max(intensityDistInt(:,:,20)))])
            colorbar
        end

        
        
%         
%         %TEST OF HEX 2 CART
%         [exposedStruct,fillfrac,threshold] = GAoptions.fillHandler.applyFill(fillChromosome,intensityDist);
%         figure
%         patched = patch(isosurface(padarray(intensityDist,[1,1,1],1e20),threshold));
%         set(patched,'FaceColor', [30 255 30]/256, 'EdgeColor', 'none');
%         view(3);
%         camlight
%         axis equal
%         lighting gouraud
%         xlim([2 size(exposedStruct,2)+1])
%         ylim([2 size(exposedStruct,1)+1])
%         zlim([2 size(exposedStruct,3)+1])
%         
%         
        
        
        
        
%         %Make hexagonal repeat of structure:
%         lx = size(intensityDist,1);
%         left = intensityDist(1:floor(lx/2),:,:);
%         right = intensityDist((floor(lx/2)+1):end,:,:);
%         top = cat(1,right,left);
%         intensityDist = cat(2,intensityDist,top);

        %Make hexagonal repeat of structure:
        ly = size(intensityDistInt,2);
        left = intensityDistInt(:,1:floor(ly/2),:);
        right = intensityDistInt(:,(floor(ly/2)+1):end,:);
        top = cat(2,right,left);
        intensityDistInt = cat(1,intensityDistInt,top);

%         
%         
%         %TEST OF HEXAGONAL REPEAT
%         [exposedStruct,fillfrac,threshold] = GAoptions.fillHandler.applyFill(fillChromosome,intensityDist);
%         figure
%         patched = patch(isosurface(padarray(intensityDist,[1,1,1],1e20),threshold));
%         set(patched,'FaceColor', [30 255 30]/256, 'EdgeColor', 'none');
%         view(3);
%         camlight
%         axis equal
%         lighting gouraud
%         xlim([2 size(exposedStruct,2)+1])
%         ylim([2 size(exposedStruct,1)+1])
%         zlim([2 size(exposedStruct,3)+1])
        
    else
        intensityDistInt = intensityDist;  %Do no interpolation if already cartesian
        
        
    end
    
    
    
    
    %Calc structure
    [exposedStruct,fillfrac,threshold] = GAoptions.fillHandler.applyFill(fillChromosome,intensityDistInt);
    
    thresholdOut = threshold
    
    %Do repeating structure:
    %exposedStructRepeat = repmat(exposedStruct,[1,1,GAoptions.num]
    
    %Export structure/runfile
    %writeLumericalRunFileSquare(GAoptions, simStruct);
    fitness = GAoptions.lumInterfaceFitness.calcCDfitness(exposedStruct);
    
    
    %transmissionRight = LumResults.transmission_right/sqrt(2);
    %transmissionLeft = LumResults.transmission_left/sqrt(2);
    %fitness = abs( transmissionRight - transmissionLeft );
end

%threshold = fixfill(reshape(intensityDist,1,[]),256,GAoptions.fill); %Calculates the threshold value that will yield desired fill fraction
%fitness = rand();



if GAoptions.runSingle == 1 %Plot simulated structure
    
   if twoPhoton == 1
      intensityDist = intensityDist.^2 / 10000; 
       
   end
    
    
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
    
    
%     max(max(intensityDist(10,:,:))) - min(min(intensityDist(10,:,:)))
%     %Plot intensity cross-section
%     figure
%     HeatMap(squeeze(intensityDist(10,:,:)) - min(min(squeeze(intensityDist(10,:,:)))) )
%     
    
    if strcmp(GAoptions.lattice, 'square') %&& ~strcmp(strtrim(GAoptions.hostname),'Daniel-netbook')
        
        %Apply skin
        if GAoptions.useSkinSingle == 1
            skin = calcSkin(exposedStruct,GAoptions.skinIterations);
            writeLumericalRunFileSkin(GAoptions, skin);
        else
            writeLumericalRunFileSquare(GAoptions, exposedStruct);
        end
        
        %Do Lumerical simulation
        %writeLumericalRunFileSquare(GAoptions, intensityDist>threshold);
        
        
        system(['fdtd-solutions -run ', GAoptions.dir, GAoptions.LumRunScript]);
        while(~exist([GAoptions.dir,GAoptions.currentLumResultsFile],'file'))
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
        
        %PLOT HEATMAP:
        if plotHeatMaps == 1
            figure
            colormap(hot)
            imagesc( intensityDist(:,:,20),[0,max(max(intensityDist(:,:,20)))])
            colorbar
            distMax  = max(max(intensityDist(:,:,20)))
            distMin  = min(min(intensityDist(:,:,20)))
            figure
            hist(reshape(intensityDist(:,:,20),1,[]),100)
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
        end
        
        %If hexagonal, convert to cartesian system for plotting and lumerical:
        intensityDist = hex2cart(intensityDist, GAoptions.cellsCart);
        
        
        %PLOT HEATMAP:
        if plotHeatMaps == 1
            figure
            colormap(hot)
            imagesc(intensityDist(:,:,20),[0,max(max(intensityDist(:,:,20)))])
            colorbar
            distMax  = max(max(intensityDist(:,:,20)))
            distMin  = min(min(intensityDist(:,:,20)))
            figure
            hist(reshape(intensityDist(:,:,20),1,[]),100)
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
        end
        
        
        
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
            while(~exist([GAoptions.dir,GAoptions.currentLumResultsFile],'file'))
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
    %acidCount = excitePAG(intensityDist,GAoptions.dimensions,GAoptions.sensSim.sensDens,GAoptions.sensSim.absCrossSection,GAoptions.sensSim.Texposure);
    acidCount = excitePAG(intensityDist,GAoptions.dimensions,GAoptions.sensSim.sensDens,GAoptions.sensSim.QYtimesAbsCrossSection,GAoptions.sensSim.Texposure);
    %intensitydist is in units of W/m^2
    acidMax = max(max(max(acidCount)))
    acidMin = min(min(min(acidCount)))
    threshold = fixfill(acidCount,256,fillfrac); %Calculates the threshold value that will yield desired fill fraction
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
    
    
    %DIFFUSION MODEL BASED ON PROFILE SIMULATION OF SU8 THICK FILM RESIST
    R = 8.617e-5; %eV/K
    T = 335; %K
    
%     Ksci0 = 2.9; %1/s   %prefactor for polymerization
%     Easci = 0.062; %eV  %activation energy of crosslinking
%     Kloss0 = 3.5; %1/s  %prefactor for polymerization termination
%     Eloss = 0.051; %1/s

%Based on newer version of article (jstage)
    Ksci0 = 5.6812; %1/s   %prefactor for polymerization
    Easci = 0.1687; %eV  %activation energy of crosslinking
    Kloss0 = 9.539; %1/s  %prefactor for polymerization termination
    Eloss = 0.3135; %eV
%     EaD = 0.1476; %eV
%     D0 = 5.0767; %1/s %but what does this mean?
    
    Ksci = Ksci0*exp(-Easci/(R*T))  %constants given in paper
    Kloss = Kloss0*exp(-Eloss/(R*T))
    
    %sadly, a value for activation energy of diffusion constant is not
    %given, assume D=D0?
    D0 = 6.922e-7; %um^2/s
    D0 = 0.3; %um^2/s  %valid for diffusion of water through SU8
    D0 = 3;
    D0 = 30;
    D0 = 300;
    D0 = 900;
    %D0 = 90000;
    %D=D0*exp(0.5/(R*T)); 
    D=D0*exp(-0.005/(R*T)); 
    l = size(acidCount,1)/GAoptions.dimensions(1); %length of cell
    probWalk = 1;
    %D=l^2/(2*tau)
    tau = l^2/(2*D)*probWalk  %choose timestep based on desired diffusivity %s
    
    %FOR RANDOM WALK DIFFUSION:
    %numRandWalkSteps = 90    %for 1x
    numRandWalkSteps = 9  %for 10x
    
    
    

    %probWalk = 0.3 %for all directions
    %probWalk = 1 %for all directions  %since timestep is based on matching diffusion set
    prob_north = 1/6; %probability of moving north given moving
    prob_south = 1/5; %probability of moving south given moving but not moving north
    prob_east = 1/4;
    prob_west = 1/3;
    prob_up = 1/2;
    prob_down = 1;
    
    
    %probCrosslink = 1 %probability a particle makes a crosslink each time  (Can make into a Bernouli test)
    
    % n_moving = zeros(Nx,Ny,Nz);
    % n_north = zeros(Nx,Ny,Nz);  %Holds the number of particles moving north from each spot
    % n_south = zeros(Nx,Ny,Nz);
    % n_east = zeros(Nx,Ny,Nz);
    % n_west = zeros(Nx,Ny,Nz);
    % n_up = zeros(Nx,Ny,Nz);
    % n_down = zeros(Nx,Ny,Nz);
    
    
    crosslinkCount = zeros(size(acidCount));  %accumulated number of crosslinks  %should be normalized to max theoretical # crosslinks
%concentration of unreacted epoxide:
    epoxideConcM = 4.96 * 1e-15 * 6.022e23; %mol/L * 1L/1000cm^3 * (1cm/10000um)^3 * 6.022e23molecules/mol -> molecules/um^3
    epoxidePerCell = epoxideConcM*(GAoptions.dimensions(1)*GAoptions.dimensions(2)*GAoptions.dimensions(3))/(GAoptions.cells(1)*GAoptions.cells(2)*GAoptions.cells(3)) %Number of molecules/cell
    epoxideCount = epoxidePerCell * ones(size(acidCount));
    
    %Acid count in model is normalized, but it is not explicit about what
    %it is normalized to
    %Assuming acid count is normalized to concentration of triaryl
    %sulfonium salt PAG concentration (0.0349M)
    baseTASSconcentration = 6.3e7; %    =0.1047M
    nTASScell = baseTASSconcentration/( size(acidCount,1)*size(acidCount,2)*size(acidCount,3) );
    %normalizedAcidConc = acidCount / baseTASSconcentration;
    
    
    for i_step = 1:numRandWalkSteps
        
        epoxideCountMaxt = max(max(max(epoxideCount)))
        epoxideCountMint = min(min(min(epoxideCount)))

        acidCountMaxt = max(max(max(acidCount)))
        acidCountMint = min(min(min(acidCount)))
        
        
        n_moving_tot = binornd( acidCount, probWalk);  %number moving
        n_moving = n_moving_tot; %number left to assign to different directions
        n_north = binornd(n_moving, prob_north);  %number moving north
        n_moving = n_moving - n_north;
        n_south = binornd(n_moving, prob_south);
        n_moving = n_moving - n_south;
        n_east = binornd(n_moving, prob_east);
        n_moving = n_moving - n_east;
        n_west = binornd(n_moving, prob_west);
        n_moving = n_moving - n_west;
        n_up = binornd(n_moving, prob_up);
        n_moving = n_moving - n_up;
        n_down = n_moving;
        
        %remove particles that moved away
        acidCount = acidCount - n_moving_tot;
        
        %Add particles that arrived
        %order of dimensions doesn't matter
        
        %east
        acidCount(1,:,:) = acidCount(1,:,:) + n_east(end,:,:);
        acidCount(2:end,:,:) = acidCount(2:end,:,:) + n_east(1:(end-1),:,:);
        %west
        acidCount(end,:,:) = acidCount(end,:,:) + n_west(1,:,:);
        acidCount(1:(end-1),:,:) = acidCount(1:(end-1),:,:) + n_west(2:end,:,:);
        %north
        acidCount(:,1,:) = acidCount(:,1,:) + n_north(:,end,:);
        acidCount(:,2:end,:) = acidCount(:,2:end,:) + n_north(:,1:(end-1),:);
        %south
        acidCount(:,end,:) = acidCount(:,end,:) + n_south(:,1,:);
        acidCount(:,1:(end-1),:) = acidCount(:,1:(end-1),:) + n_south(:,2:end,:);
        %up
        acidCount(:,:,1) = acidCount(:,:,1) + n_up(:,:,end);
        acidCount(:,:,2:end) = acidCount(:,:,2:end) + n_up(:,:,1:(end-1));
        %down
        acidCount(:,:,end) = acidCount(:,:,end) + n_down(:,:,1);
        acidCount(:,:,1:(end-1)) = acidCount(:,:,1:(end-1)) + n_south(:,:,2:end);
        
        
        %Calculate num crosslinks
        %crosslinkCount = crosslinkCount + acidCount;  %Assuming each acid produces exactly 1 crosslink per step
        %crosslinkCount = crosslinkCount + Ksci.*acidCount.*(1-crosslinkCount).*tau;
        %crosslinkCount = crosslinkCount + Ksci.*(acidCount/nTASScell).*(1-crosslinkCount).*tau;  %Doing a normalization of acid concentration relative to TASS, since the model uses values for normalized acid concentration
 %exponential decay of unreacted epoxide:
        epoxideCount = epoxideCount.*exp(-Ksci.*(acidCount/nTASScell).*tau); 
        
        %Calculate 'acid' loss (polymerization termination)
        %        acidCount = acidCount - Kloss.*(acidCount/nTASScell).*tau;
        %Should instead be an exponential loss:
%        acidCount = acidCount -  binornd( acidCount, Kloss*tau );
        

        

        
    end
    
    
    crosslinkCount = epoxidePerCell - epoxideCount;
    crosslinkMax = max(max(max(crosslinkCount)))
    crosslinkMin = min(min(min(crosslinkCount)))
    
    
%     %DO SMOOTHING//DIFFUSION (WITH SUMMING)
%     diffRate = 5;
%     smoothed = zeros([size(acidCount),8]);
%     smoothed(:,:,:,1) = smooth3(acidCount,'gaussian',diffRate);
%     for smooth_i = 2:size(smoothed,4)
%         smoothed(:,:,:,smooth_i) = smooth3(smoothed(:,:,:,smooth_i-1),'gaussian',diffRate);
%         disp(smooth_i)
%     end
%     crosslinkCount = sum(smoothed,4);
    crosslinkThreshold = fixfill(crosslinkCount, 256, fillfrac)
    figure
    patched = patch(isosurface(padarray(crosslinkCount,[1,1,1],1e99),crosslinkThreshold));
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([2 size(acidCount,2)+1])
    ylim([2 size(acidCount,1)+1])
    zlim([2 size(acidCount,3)+1])

    
    %PLOT HEATMAP:
    if plotHeatMaps == 1
        figure
        colormap(hot)
        imagesc(crosslinkCount(:,:,20),[0,max(max(crosslinkCount(:,:,20)))])
        colorbar
    end
    
    
    
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

    function outindex = applyPeriod(inindex,maxindex)
        if inindex < 0
            outindex = maxindex;
        elseif inindex > maxindex
            outindex = 0;
        end
    end

        
end
function runGA_8beamVarLayer

%GA version 2: modular code.
%This runfile assumes a 4-beam symmetric configuration (no central mode/ non-normal incidence)
%Each module (e.g. the grating generator, fitness measurer, each has own
%piece of GAoptions data

%Different code for different computers
[idum,hostname] = system('hostname')
GAoptions.hostname = strtrim(hostname);

%%%%% Genetic Algorithm Options %%%%%
    GAoptions.popSize = 300;
    GAoptions.numGen = 25;
    GAoptions.elite = 1;
    GAoptions.numRepetitions = 1; %Number of times to repeat the GA
    %Options for built-in GA algorithm:
    if strcmp(strtrim(hostname),'lotus-bud')
        UseParallelVar = true;
    elseif strcmp(strtrim(hostname),'berzerk')
        UseParallelVar = 'always'
    elseif strcmp(strtrim(hostname),'Daniel-netbook')
        UseParallelVar = 'never';
    end
    GAoptions.optimset = gaoptimset('PopulationSize', GAoptions.popSize,...
        'PopulationType', 'bitstring',...
        'EliteCount', GAoptions.elite,...
        'CrossoverFraction', 0.95,...
        'CrossoverFcn', @crossoverscattered,...
        'Generations', GAoptions.numGen,...
        'SelectionFcn', @selectionroulette,...
        'TolFun', 0,...
        ... %'MutationFcn', @mutationfcn,...  use default gaussian mutation distribution
        'FitnessScalingFcn', @fitscalingshiftlinear,...
        'StallTimeLimit', 120*3600, ... %five days
        'Display', 'diagnose', ...
        'PlotFcns',@gaplotbestf, ...
        'UseParallel',UseParallelVar, ...  
        'StallGenLimit', 20);
    %Get and record random number generator
    GAoptions.randomStream = RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(GAoptions.randomStream);
    
    
    
    %%%%% Filenames %%%%%
    if strcmp(strtrim(hostname),'lotus-bud')
        GAoptions.dir = ['/home/daniel/GA_8beamReflect_helix_data/','GA_',datestr(now,'mm-dd-yy_HH:MM'),'_N',int2str(GAoptions.popSize), '/']  %directory that all this goes into
    elseif strcmp(strtrim(hostname),'berzerk')
% figure
% patched = patch(isosurface(padarray(intensityDist,[1,1,1],100),threshold));
% set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
% view(3);
% camlight
% lighting gouraud
% xlim([1 size(intensityDist,1)])
        GAoptions.dir = ['/home/danielbacon-brown/GA_8beamReflect_helix_data/','GA_',datestr(now,'mm-dd-yy_HH:MM'),'_N',int2str(GAoptions.popSize), '/']  %directory that all this goes into
    elseif strcmp(strtrim(hostname),'Daniel-netbook')
        GAoptions.dir = ['C:/Users/daniel/GA_8beamReflect_helix_data/','GA_ver2_',datestr(now,'mm-dd-yy_HH;MM'),'_N',int2str(GAoptions.popSize), '/']  %directory that all this goes into
    end
    mkdir(GAoptions.dir);%Create folder if it doesn't exist yet
    GAoptions.GARecordFileBase = 'GAResults';  %Stores results
    GAoptions.currentLumSave = 'LumSave.fsp';
    GAoptions.LumRunScript = 'LumRun.lsf';
    GAoptions.currentLumResultsFile = 'LumResults'
    
    
    
%     %%%%% Grating Generation %%%%%
%     %gratingOptions.periodicity = 0.3572;  %um  %should give c=a
%     gratingOptions.periodicity = 0.3928; %um %should give c=a
%     gratingOptions.n_filled = 1.58 %1.58;  %refractive index of grating material
%     gratingOptions.n_void = 1;  %refractive index of void space
%     gratingOptions.NpixelX = 30; %Number of blocks
%     gratingOptions.NpixelY = 30; %Number of lines
%     gratingOptions.chromNthickness = 8;  %Number of chromosome positions for each variable
%     gratingOptions.chromNpointX = 8;
%     gratingOptions.chromNpointY = 8;
%     gratingOptions.chromThicknessMax = 0.2; %Maximum thickness of chromosome %um
%     gratingOptions.penRadius = 0.050; %pen radius (um)
%     gratingOptions.numLines = 2; %Number of lines to draw
%     gratingOptions.order = make_order(4, 'hexagonal')  %Diffracted orders
% % temporder =  make_order(4, 'hexagonal')
% % tempm1 = temporder(:).m1
% % tempm2 = temporder(:).m2
% % disp(num2str(tempm1))
% % disp(num2str(tempm2))
% % save('temporderout','temporder')
% %gratingOptions.order = 
%     GAoptions.gratingOptions = gratingOptions;
%     GAoptions.gratingFunction = GratingDrawLineHexagonal(gratingOptions);
%     disp(GAoptions.gratingFunction)
    

    %%%%%% Lattice Dimensions %%%%%
    GAoptions.laserWavelength = 0.532; %um
    GAoptions.C_over_A = 0.75;    %Max C/A for air gap is 0.578 %for PDMS prism, max C/A = 1.396
    GAoptions.lattice = 'square';
    GAoptions.n_PR = 1.58; %refractive index of the photoresist (SU8)
    GAoptions.n_substrate = 1.5; %Glass slide as substrate
    GAoptions.n_prism = 1.5; %PDMS prism
    %GAoptions.n_prism = 1; %no prism
    GAoptions.n_gratingVoid = 1; %assuming vacuum-SU8 grating
    %This assumes the 4-beam symmetric configuration
    GAoptions.period = GAoptions.laserWavelength/(2*GAoptions.n_PR) * sqrt(2+1/(GAoptions.C_over_A^2))
    GAoptions.cells = floor(25*[1,1,GAoptions.C_over_A]); %number of 'ticks' in each dimension
    GAoptions.dimensions = [GAoptions.period, GAoptions.period, GAoptions.period*GAoptions.C_over_A]; %dimensions of unit cell
    %Vectors describing periodicity:
    GAoptions.u = [GAoptions.period,0,0];
    GAoptions.v = [0,GAoptions.period,0];
    GAoptions.w = [0,0,GAoptions.period*GAoptions.C_over_A];
    GAoptions.isAirGap = 0;
    


    
    
    
    
    
    %%%%% Incident light %%%%%
    incidentLightOptions.wavelength = GAoptions.laserWavelength ;  %um
    incidentLightOptions.chromNpsi = 8;
    incidentLightOptions.chromNchi = 8;
    incidentLightOptions.n_interference = GAoptions.n_PR;
    incidentLightOptions.n_incidence = GAoptions.n_prism; %The refractive index of the material that the plane wave is in
    incidentLightOptions.period = GAoptions.period;
    incidentLightOptions.C_over_A = GAoptions.C_over_A;
    incidentLightOptions.beamPowerDens = 20935  %W/m^2      %=20935W/m^2
    GAoptions.incidentLightOptions = incidentLightOptions;
    GAoptions.incidentLightFunction = IncidentLightAngled4BeamSymmetric(incidentLightOptions);


%For square GratingGrid
    gratingOptions.NblockX = 3;
    gratingOptions.NblockY = 3;
    gratingOptions.chromNspacingX = 8;
    gratingOptions.chromNspacingY = 8;
    gratingOptions.chromNSU8thickness = 8; 
    gratingOptions.n_filled = GAoptions.n_PR; %1.58;  %refractive index of grating material
    gratingOptions.n_void = GAoptions.n_gratingVoid;  %refractive index of void space of grating
    %gratingOptions.periodicity = 0.3928; %um %should give c=a
    %gratingOptions.periodicity = 0.5334; %um should give c=a*2
    gratingOptions.C_over_A = GAoptions.C_over_A;
    gratingOptions.periodicity = GAoptions.period; %getA(gratingOptions.c_over_a, gratingOptions.n_filled, incidentLightOptions.wavelength) * gratingOptions.n_filled/gratingOptions.n_void;  %Need to scale by refractive index because of critical angle at air-SU8 interface
    gratingOptions.chromNthickness = 8;  %Number of chromosome positions for each variable
    gratingOptions.thicknessMax = 0.3; %Maximum thickness of grating %um
    gratingOptions.spacingMin = 0.2; %Minimum block width 
    gratingOptions.order = make_order(4, 'hexagonal');
    gratingOptions.SU8thicknessMax = 10; %um
    gratingOptions.SU8thicknessMin = 5; %um
    GAoptions.gratingOptions = gratingOptions;
    GAoptions.gratingFunction = GratingGridSquare(gratingOptions);
    disp(GAoptions.gratingFunction)




    
    
       
    %Note: for S4, the cells for u and v vectors need to be equal
    GAoptions.cells = floor(25*[1,1,GAoptions.C_over_A]);

    GAoptions.dimensions = [GAoptions.period, GAoptions.period, GAoptions.period*GAoptions.C_over_A];

    %%%%% Target Structure %%%%%
    %targetStructureOptions.cells = [40,floor(40*sqrt(3)/2), 40];  %c=a
    %targetStructureOptions.cells = [40,floor(40*sqrt(3)/2), 40*1.4]; c=a*1.4
    %Calc c from a:
    %c = gratingOptions.periodicity/ ...
    %    ( gratingOptions.periodicity*gratingOptions.n_filled/incidentLightOptions.wavelength - sqrt( (gratingOptions.periodicity*gratingOptions.n_filled/incidentLightOptions.wavelength)^2 - 4/3))
    
    %targetStructureOptions.cells = floor([40,40*sqrt(3)/2, 40*gratingOptions.c_over_a]);
    targetStructureOptions.cells = GAoptions.cells;
    %targetStructureOptions.dimensions = [GAoptions.gratingOptions.periodicity, GAoptions.gratingOptions.periodicity*sqrt(3)/2, GAoptions.gratingOptions.periodicity];
    %targetStructureOptions.dimensions = [GAoptions.gratingOptions.periodicity, GAoptions.gratingOptions.periodicity*sqrt(3)/2, gratingOptions.periodicity*gratingOptions.c_over_a];
    targetStructureOptions.u = GAoptions.u;
    targetStructureOptions.v = GAoptions.v;
    targetStructureOptions.w = GAoptions.w;
    targetStructureOptions.radb = GAoptions.period/4; %radius / distance of middle of helix to axis (um)
    targetStructureOptions.radl = GAoptions.period/8; %radius of the wire (um)
    targetStructureOptions.relativeZ = GAoptions.C_over_A; %ratio of height to width of ellipsoidal 'pen'
    GAoptions.targetStructureOptions = targetStructureOptions;
    GAoptions.targetStructure = generateHelixStructureParallelogram(targetStructureOptions)  %Create helix
    targetfill = 1 - sum(sum(sum(GAoptions.targetStructure))) / (size(GAoptions.targetStructure,1)*size(GAoptions.targetStructure,1)*size(GAoptions.targetStructure,3));

    %This is a larger helix that counts as negative for filled spots outside of it
    %exclusionStructureOptions.cells = floor([40,40*sqrt(3)/2, 40*gratingOptions.c_over_a]);
    exclusionStructureOptions.cells = GAoptions.cells;    
    %exclusionStructureOptions.dimensions = [GAoptions.gratingOptions.periodicity, GAoptions.gratingOptions.periodicity*sqrt(3)/2, gratingOptions.periodicity*gratingOptions.c_over_a];
    exclusionStructureOptions.u = GAoptions.u;
    exclusionStructureOptions.v = GAoptions.v;
    exclusionStructureOptions.w = GAoptions.w;
    exclusionStructureOptions.radb = GAoptions.period/4;
    exclusionStructureOptions.radl = GAoptions.period/4;
    exclusionStructureOptions.relativeZ = GAoptions.C_over_A;  %ratio of height to width of ellipsoidal 'pen'
    GAoptions.exclusionStructureOptions = exclusionStructureOptions;
    GAoptions.exclusionStructure = generateHelixStructureParallelogram(exclusionStructureOptions)  %Helix
    exclusionfill = 1 - sum(sum(sum(GAoptions.exclusionStructure))) / (size(GAoptions.exclusionStructure,1)*size(GAoptions.exclusionStructure,1)*size(GAoptions.exclusionStructure,3));
    
    GAoptions.fill = (targetfill+exclusionfill)/2;  %Matches interference pattern fill to average of the target and exclusion structures. Consider adding fill factor to chromosome

    %1 at the edges of the unit cell to reduce fitness of x-y continuous
    %structures
    edgeExclusionStructure = zeros(GAoptions.cells);
    edgeExclusionStructure(1,:,:) = 1;
    edgeExclusionStructure(end,:,:) = 1;
    edgeExclusionStructure(:,1,:) = 1;
    edgeExclusionStructure(:,end,:) = 1;
    GAoptions.edgeExclusionStructure = edgeExclusionStructure;
    
    
    % Not used for S4-based sims
%     %%%%% Calc structure %%%%%  
%     %GAoptions.dimensions = [GAoptions.gratingOptions.periodicity, GAoptions.gratingOptions.periodicity*sqrt(3)/2, GAoptions.gratingOptions.periodicity];    %MAKE SURE DIMENSIONS ARE ACCURATE
%     GAoptions.dimensions = targetStructureOptions.dimensions;
%     calcStructureOptions.chromNoffsetX = 8;%ratio of height to width of ellipsoidal 'pen'
%     calcStructureOptions.chromNoffsetY = 8;
%     calcStructureOptions.chromNoffsetZ = 8;
%     GAoptions.calcStructureOptions = calcStructureOptions;
%     %GAoptions.calcStructureFunction = CalculateIntensitySimple(calcStructureOptions); %Does basic interference calculation, monochromatic, no PAG, no reflection
    
    %%%%% Offset %%%%%
    offsetOptions.offsetType = 'parallogram';
    offsetOptions.chromNoffsetX = 8;
    offsetOptions.chromNoffsetY = 8;
    offsetOptions.chromNoffsetZ = 8;
    offsetOptions.cells = GAoptions.cells;
    GAoptions.offsetConductor = OffsetConductor(offsetOptions);
    
    
    
    %%%% Calc Fitness %%%%%
    GAoptions.calcFitness = @calcVolumetricMatchExclusion;
    %GAoptions.fitnessFunction.targetstructure = 'helix_1to1';
    
    
%     %S4 layer data:
%     S4interfaceOptions.materialNames =   ['Vacuum',  'SU8',  'Prism',    'Glass',        'ITO'];
%     S4interfaceOptions.materialRI =      [1,         1.58,   n_prism,    n_substrate,    1.94-0.046i];      
%     
%     S4interfaceOptions.layerNames =        ['Front',   'TCO',  'PR',   'Grating',  'Back'];
%     S4interfaceOptions.layerMaterials =    ['Prism',   'ITO',  'SU8',  'Vacuum',   'Vacuum'];
%     S4interfaceOptions.layerThicknesses =  [0,         0.1,    0,      0,          0];      
%     

    %%%%% Materials %%%%%
    S4interfaceOptions.materialNames = { ...
        'Vacuum', ...
        'Glass', ...
        'ITO', ...
        'SU8'};
    lengthout = length(S4interfaceOptions.materialNames)
    S4interfaceOptions.materialRI = [ ...
        1, ...
        1.5, ...
        1.94 + 0.046i, ...
        1.58 ];
    
    %%%%% Layers: %%%%%
    chromNlayer = 8
    S4interfaceOptions.layers(1) = Layer('Front','Glass',0);
    S4interfaceOptions.layers(2) = Layer('TCO','ITO',-1,0.03,0.2,chromNlayer);
    S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
    S4interfaceOptions.layers(4) = Layer('Grating','Vacuum', -1, 0,gratingOptions.thicknessMax,chromNlayer); 
    S4interfaceOptions.layers(5) = Layer('Back','Vacuum', 0);


    %Writer for S4 runfile
    S4interfaceOptions.dimensions = GAoptions.dimensions;
    S4interfaceOptions.cells = GAoptions.cells;
    S4interfaceOptions.isAirGap = GAoptions.isAirGap;
    S4interfaceOptions.gratingCoatingMetal = 'gold';
    S4interfaceOptions.gratingCoatingThickness  = 0.03; %um
    S4interfaceOptions.n_glass = GAoptions.n_substrate;
    S4interfaceOptions.n_prism = GAoptions.n_prism;
    GAoptions.S4interface = S4interfaceSquareGeneral(S4interfaceOptions)
    GAoptions.S4interfaceOptions = S4interfaceOptions;
    %GAoptions.S4interface = S4interfaceSquareReflect2(S4interfaceOptions)
    %GAoptions.S4interface = S4interfaceSquareTransmit(S4interfaceOptions)
    
    %%%%% Do FDTD %%%%%
    fdtd.simulationTime = 500; %fs
    fdtd.shrinkFactor = 0.9; %Thickness after development / thickness during interference
    fdtd.repeatUnits = 3; %Number of unit cells to do FDTD simulation
    fdtd.meshAccuracy = 1;
    fdtd.minSourceWL = 0.5e-6; %m
    fdtd.maxSourceWL = 4e-6;  %m
    fdtd.minMeasWL = 0.5e-6; %m
    fdtd.maxMeasWL = 4e-6;  %m
    fdtd.numMeasWL = 500;
    fdtd.addCubesDirectly = 1;
    fdtd.SU8matrix = 1;
    GAoptions.fdtd = fdtd;
    
    
    %%%%% Do sensitizer simulation %%%%%
    sensSim.sensDens = 8290000 *4; %Density of sensitizer molecules per um^3
    %abs_prob_per_photon = 2.91e-7; %nm^2/molecule %BCPI %based on 5-9-15 UV-vis measurements 
    sensSim.absCrossSection = 2.91e-7; %nm^2/molecule %BCPI  %MAY NEED TO UPDATE
    sensSim.Texposure = 200; %s
    GAoptions.sensSim = sensSim;
    
    
    
    %Get layer chromosome length
    chromNlayer=0;
    for i=1:length(GAoptions.S4interfaceOptions.layers)
        chromNlayer=chromNlayer+GAoptions.S4interfaceOptions.layers(i).getChromosomeSize();
    end
    %Creation of 'GAproblem' for built-in Matlab GA
    GAoptions.GAproblem.options = GAoptions.optimset;  %options set above
    GAoptions.GAproblem.nvars = ...
        + GAoptions.gratingFunction.getChromosomeSize()  ...
        + GAoptions.incidentLightFunction.getChromosomeSize()  ...
        ... + GAoptions.calcStructureFunction.getChromosomeSize() ...
        + GAoptions.offsetConductor.getChromosomeSize() ...
        + chromNlayer ...
        ;
    GAoptions.GAproblem.fitnessfcn = @gfit;  %fitness function
    GAoptions.GAproblem.Aineq = [];  %A matrix for linear inequality constraints
    GAoptions.GAproblem.Bineq = [];  %B vector for linear inequality constraints
    GAoptions.GAproblem.Aeq = [];  %A matrix for linear equality constraints
    GAoptions.GAproblem.Beq = [];  %B vector for linear inequality constraints
    GAoptions.GAproblem.LB = [];  %lower bound on 'x'
    GAoptions.GAproblem.VB = [];  %upper bound on 'x'
    GAoptions.GAproblem.non1con = [];  %nonlinear constraint function
    
    

    cd(GAoptions.dir)
    
tic
[chromosome, fitness, reason, output, pop, scores] = ga(GAoptions.GAproblem);
save([GAoptions.dir,GAoptions.GARecordFileBase,datestr(now,'mm-dd-yy_HH;MM')],'chromosome','fitness','reason','output','pop','scores','GAoptions');
toc

% %DO BEST CHROMOSOME:
% bestchrom = chromosome
doPlots=1;
fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome,doPlots) %Should also work for reflections





% [gratingChromosome, incidentLightChromosome, offsetChromosome] = splitChromosome(chromosome,[ ...
%     GAoptions.gratingFunction.getChromosomeSize(), ...
%     GAoptions.incidentLightFunction.getChromosomeSize(),  ...
%     GAoptions.calcStructureFunction.getChromosomeSize()]);
% %Generate Grating
% grating = GAoptions.gratingFunction.generateGrating(gratingChromosome);
% %GAoptions.gratingFunction.plotGrating(grating)
% 
% %Incident Field
% [incidentField,incidentEsp] = GAoptions.incidentLightFunction.generateField(incidentLightChromosome);
% 
% %Do RCWA analysis
% [~, scat_field, ~] = gdc(grating,incidentField,GAoptions.gratingOptions.order);
% %Calculate propagating diffracted beams
% [E231, k231] = field_convert(incidentEsp,scat_field);  %DO OWN REWRITE
% 
% %Calculate intensity distribution:
% intensityDist = GAoptions.calcStructureFunction.calcIntensity( E231,k231,GAoptions.dimensions,GAoptions.targetStructureOptions.cells,offsetChromosome);
% 
% 
% [fitness,threshold] = calcVolumetricMatch(GAoptions.targetStructure, intensityDist)
% 
% figure
% patched = patch(isosurface(padarray(intensityDist,[1,1,1],100),threshold));
% set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
% view(3);
% camlight
% lighting gouraud
% xlim([1 size(intensityDist,1)])
% ylim([1 size(intensityDist,2)])
% zlim([1 size(intensityDist,3)])




        
        
        
%%%%% Declaration of Fitness Function %%%%%
    function fitness = gfit(chromosome)  %This is done so that it can pass GAoptions to the fitness function.
        tic
        fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome)
        toc 
%         %Splits chromosome according to each module
%         [gratingChromosome, incidentLightChromosome, offsetChromosome] = splitChromosome(chromosome,[ ...
%             GAoptions.gratingFunction.getChromosomeSize(), ...
%             GAoptions.incidentLightFunction.getChromosomeSize(),  ...
%             GAoptions.calcStructureFunction.getChromosomeSize()]);
%         
%         %Generate Grating
%         grating = GAoptions.gratingFunction.generateGrating(gratingChromosome);
%         %GAoptions.gratingFunction.plotGrating(grating)
%         
%         %Incident Field
%         [incidentField,incidentEsp] = GAoptions.incidentLightFunction.generateField(incidentLightChromosome);
%         
%         
%         %Do RCWA analysis
%         [~, scat_field, ~] = gdc(grating,incidentField,GAoptions.gratingOptions.order);
%         
%         %Calculate propagating diffracted beams
%         [E231, k231] = field_convert(incidentEsp,scat_field);  %DO OWN REWRITE
%         
%         %Calculate intensity distribution:
%         intensityDist = GAoptions.calcStructureFunction.calcIntensity( E231,k231,GAoptions.dimensions,GAoptions.targetStructureOptions.cells,offsetChromosome);
%         
%         
%         %Calculate Fitness
%         [fitness,threshold] = calcVolumetricMatch(GAoptions.targetStructure, intensityDist);
%         
        
        
    end




end







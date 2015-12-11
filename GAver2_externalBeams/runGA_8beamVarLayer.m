function runGA_8beamVarLayer

%GA version 2: modular code.
%This runfile assumes a 4-beam symmetric configuration (no central mode/ non-normal incidence)
%Each module (e.g. the grating generator, fitness measurer, each has own
%piece of GAoptions data

%Different code for different computers
[idum,hostname] = system('hostname')
GAoptions.hostname = strtrim(hostname);


    %%%%% Genetic Algorithm Options %%%%%
    GAoptions.popSize = 6;
    GAoptions.numGen = 1;
    GAoptions.elite = 1;
    GAoptions.numRepetitions = 1; %Number of times to repeat the GA
    
    
    %%%%% FITNESS TYPE %%%%%
    GAoptions.fitnessType = 'structure';
    %GAoptions.fitnessType = 'fdtd';
    GAoptions.useExclusion = 0; %for 'structure'
    GAoptions.useEdgeExclusion = 0; %for 'structure'
    
    
    %%%%% Filenames %%%%%
    GAoptions.LumRunScript = 'LumRunScript.lsf'; %Runs for single
    GAoptions.currentLumSave = 'currentLumSave.lsf';
    GAoptions.currentLumResultsFile = 'currentLumResultsFile.mat';
    if strcmp(strtrim(hostname),'lotus-bud')
        GAoptions.dir = ['/home/daniel/GA_8beamReflect_helix_data/','GA_',datestr(now,'mm-dd-yy_HH:MM'),'_N',int2str(GAoptions.popSize), '/']  %directory that all this goes into
    elseif strcmp(strtrim(hostname),'berzerk')
        GAoptions.dir = ['/home/danielbacon-brown/GA_8beamReflect_helix_data/','GA_',datestr(now,'mm-dd-yy_HH:MM'),'_N',int2str(GAoptions.popSize), '/']  %directory that all this goes into
    elseif strcmp(strtrim(hostname),'Daniel-netbook')
        GAoptions.dir = ['C:/Users/daniel/GA_8beamReflect_helix_data/','GA_ver2_',datestr(now,'mm-dd-yy_HH;MM'),'_N',int2str(GAoptions.popSize), '/']  %directory that all this goes into
    end
    mkdir(GAoptions.dir);%Create folder if it doesn't exist yet
    GAoptions.GARecordFileBase = 'GAResults';  %Stores results

    
    

    %%%%%% Lattice Dimensions %%%%%
    GAoptions.normalIncidence = 0; 
    GAoptions.laserWavelength = 0.532; %um
    GAoptions.C_over_A = 1;    %Max C/A for air gap is 0.578 %for PDMS prism, max C/A = 1.396 %Will be overwritten if normal incidence
    GAoptions.lattice = 'square';
    %GAoptions.lattice = 'hexagonal';
    GAoptions.n_PR = 1.58; %refractive index of the photoresist (SU8)
    GAoptions.n_PDMS = 1.43; 
    GAoptions.n_gratingVoid = 1; %assuming vacuum-SU8 grating
    GAoptions.n_glassSubstrate = 1.48;   %Glass slide as substrate
    GAoptions.n_front = GAoptions.n_glassSubstrate; %GAoptions.n_PDMS;  %refractive index for incident light



    %%%%% Materials %%%%%
    chromNmaterial = 8;
    S4interfaceOptions.materials(1) = Material('Vacuum',1);
    %S4interfaceOptions.materials(2) = Material('Glass',-1,1.47,1.52,chromNmaterial);  %Variable refractive index glass
    S4interfaceOptions.materials(2) = Material('Glass',GAoptions.n_glassSubstrate);
    S4interfaceOptions.materials(3) = Material('ITO',1.94+0.046i);
    S4interfaceOptions.materials(4) = Material('SU8',GAoptions.n_PR);
    S4interfaceOptions.materials(5) = Material('PDMS',GAoptions.n_PDMS);
    S4interfaceOptions.materials(6) = Material('Graphene',2.6793+1.2227i);  %For single layer
    
            
    
    %%%%% Layers: %%%%%    
    chromNlayer = 8;
    
%     %Prism-coupled, glass->ITO->SU8->grating->vacuum
%     S4interfaceOptions.layers(1) = Layer('Front','Glass',0);
%     S4interfaceOptions.layers(2) = Layer('TCO','ITO',-1,0.015,0.15,chromNlayer);
%     S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
%     S4interfaceOptions.layers(4) = Layer('Grating','Vacuum', -1, 0,0.3,chromNlayer); 
%     S4interfaceOptions.layers(5) = Layer('Back','Vacuum', 0);
    
    %Prism-coupled, glass->SU8->grating->vacuum    --requires some
    %electroless deposition
    S4interfaceOptions.layers(1) = Layer('Front','Glass',0);
    S4interfaceOptions.layers(2) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
    S4interfaceOptions.layers(3) = Layer('Grating','Vacuum', -1, 0,0.3,chromNlayer); 
    S4interfaceOptions.layers(4) = Layer('Back','Vacuum', 0);
    
%     %Prism-coupled, glass->graphene->SU8->grating->vacuum
%     S4interfaceOptions.layers(1) = Layer('Front','Glass',0);
%     S4interfaceOptions.layers(2) = Layer('Graphene','Graphene',0.00034);  %Single layer
%     S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
%     S4interfaceOptions.layers(4) = Layer('Grating','Vacuum', -1, 0,0.3,chromNlayer); 
%     S4interfaceOptions.layers(5) = Layer('Back','Vacuum', 0);
    
%     %Incident on air.  Air->grating->SU8->graphene->Glass
%     S4interfaceOptions.layers(1) = Layer('Front','Vacuum',0);
%     S4interfaceOptions.layers(2) = Layer('Grating','Vacuum', -1, 0,0.3,chromNlayer);
%     S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
%     S4interfaceOptions.layers(4) = Layer('Graphene','Graphene',0.0034);
%     S4interfaceOptions.layers(5) = Layer('Back','Glass', 0);
    
    

%     %Incident on PDMS.  PDMS->grating->SU8->ITO->Glass
%     S4interfaceOptions.layers(1) = Layer('Front','PDMS',0);
%     S4interfaceOptions.layers(2) = Layer('Grating','Vacuum', -1, 0,0.2,chromNlayer);
%     S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
%     S4interfaceOptions.layers(4) = Layer('TCO','ITO',-1,0.015,0.15,chromNlayer);
%     S4interfaceOptions.layers(5) = Layer('Back','Glass', 0);
    
    
%     %Incident on air.  Air->grating->SU8->ITO->Glass
%     S4interfaceOptions.layers(1) = Layer('Front','Vacuum',0);
%     S4interfaceOptions.layers(2) = Layer('Grating','Vacuum', -1, 0,0.3,chromNlayer);
%     S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
%     S4interfaceOptions.layers(4) = Layer('TCO','ITO',-1,0.015,0.15,chromNlayer);
%     S4interfaceOptions.layers(5) = Layer('Back','Glass', 0);
    
%     %Incident on air.  No TCO.  Air->grating->SU8->Glass
%     S4interfaceOptions.layers(1) = Layer('Front','Vacuum',0);
%     S4interfaceOptions.layers(2) = Layer('Grating','Vacuum', -1, 0,0.3,chromNlayer);
%     S4interfaceOptions.layers(3) = Layer('PrInterference','SU8',-1,5,15,chromNlayer);
%     S4interfaceOptions.layers(4) = Layer('Back','Glass', 0);
    %013
    
    
    
    %Periodicity
    if GAoptions.normalIncidence %Normally incident light
        
        if strcmp(GAoptions.lattice, 'square')
            GAoptions.period = 3 * GAoptions.laserWavelength/(sqrt(2)*2*GAoptions.n_PR);
            GAoptions.C_over_A = sqrt(2);
        elseif strcmp(GAoptions.lattice, 'hexagonal')
            GAoptions.period = sqrt(3)*GAoptions.laserWavelength/ (sqrt(2)*GAoptions.n_PR);
            GAoptions.C_over_A = sqrt(3)/sqrt(2);
        end
        
        
    else   %Non-normal incidence
        
        if strcmp(GAoptions.lattice, 'square')
            GAoptions.period = GAoptions.laserWavelength/(2*GAoptions.n_PR) * sqrt(2+1/(GAoptions.C_over_A^2));  %for c/a<0.578, light can be coupled in without prism, but then can get transmitted modes in grating
        elseif strcmp(GAoptions.lattice, 'hexagonal')
            GAoptions.period = GAoptions.laserWavelength/GAoptions.n_PR/2 * sqrt( 1/GAoptions.C_over_A^2 + (4/3)^2); %for c/a<0.613, light can be coupled in without prism
        end

    end
    
    %Cell # and dimensions
    GAoptions.repeatingUnits = 1;  %Used as a test for periodicity, should be 1 for actual optimization (except for fdtd test)
    GAoptions.cells = floor(25*[1,1,GAoptions.C_over_A*GAoptions.repeatingUnits]); %number of 'ticks' in each dimension
    if strcmp(GAoptions.lattice, 'square')
        GAoptions.dimensions = [GAoptions.period, GAoptions.period, GAoptions.period*GAoptions.C_over_A]; %dimensions of unit cell
    elseif strcmp(GAoptions.lattice, 'hexagonal')
        GAoptions.dimensions = [GAoptions.period, GAoptions.period*sqrt(3)/2, GAoptions.period*GAoptions.C_over_A]; %dimensions of unit cell
        GAoptions.cellsCart = floor(20*[1,sqrt(3)/2,GAoptions.C_over_A*GAoptions.repeatingUnits]) %Cells for when converting to cartesian coordinates (plotting and Lumerical)
    end
    
    %Vectors describing periodicity (for target structure)
    if strcmp(GAoptions.lattice, 'square')
        GAoptions.u = [GAoptions.period,0,0];
        GAoptions.v = [0,GAoptions.period,0];
        GAoptions.w = [0,0,GAoptions.period*GAoptions.C_over_A];
    elseif strcmp(GAoptions.lattice, 'hexagonal')
        GAoptions.u = [GAoptions.period,0,0];
        GAoptions.v = [GAoptions.period/2,GAoptions.period*sqrt(3)/2,0];
        GAoptions.w = [0,0,GAoptions.period*GAoptions.C_over_A];
    end


    
    
    
    
    
    %%%%% Incident light %%%%%
    incidentLightOptions.normalIncidence = GAoptions.normalIncidence;
    incidentLightOptions.wavelength = GAoptions.laserWavelength ;  %um
    incidentLightOptions.lattice = GAoptions.lattice;
    incidentLightOptions.chromNpsi = 8;
    incidentLightOptions.chromNchi = 8;
    incidentLightOptions.n_interference = GAoptions.n_PR;
    incidentLightOptions.n_incidence = GAoptions.n_front; %The refractive index of the material that the plane wave is in
    incidentLightOptions.period = GAoptions.period;
    incidentLightOptions.C_over_A = GAoptions.C_over_A;
    incidentLightOptions.beamPowerDens = 20935  %W/m^2      %=20935W/m^2
    GAoptions.incidentLightOptions = incidentLightOptions;
    %GAoptions.incidentLightFunction = IncidentLightAngled(incidentLightOptions);
    GAoptions.incidentLightFunction = IncidentLightGeneral(incidentLightOptions);

%For  GratingGrid
    gratingOptions.NblockX = 4;
    gratingOptions.NblockY = 4;
    gratingOptions.chromNspacingX = 8;
    gratingOptions.chromNspacingY = 8;
    gratingOptions.chromNSU8thickness = 8; 
    gratingOptions.n_filled = GAoptions.n_PR; %1.58;  %refractive index of grating material
    gratingOptions.n_void = GAoptions.n_gratingVoid;  %refractive index of void space of grating
    gratingOptions.C_over_A = GAoptions.C_over_A;
    gratingOptions.periodicity = GAoptions.period; 
    gratingOptions.spacingMin = 0.15; %Minimum block width %relative to periodicity
    GAoptions.gratingOptions = gratingOptions;
    if strcmp(GAoptions.lattice, 'square')
        GAoptions.gratingFunction = GratingGridSquare(gratingOptions);
    elseif strcmp(GAoptions.lattice, 'hexagonal')
        GAoptions.gratingFunction = GratingGridHexagonal(gratingOptions);
    end
    disp(GAoptions.gratingFunction)


       

    if strcmp(GAoptions.lattice,'square')
        GAoptions.dimensions = [GAoptions.period, GAoptions.period, GAoptions.period*GAoptions.C_over_A];
    elseif strcmp(GAoptions.lattice,'hexagonal')
        GAoptions.dimensions = [GAoptions.period, GAoptions.period*sqrt(3)/2, GAoptions.period*GAoptions.C_over_A];
    end

    
    %%%%% Target Structure %%%%%
    
    %targetStructureOptions.structureType = 'doublehelix';
    targetStructureOptions.structureType = 'helix';

    targetStructureOptions.cells = GAoptions.cells;
    targetStructureOptions.u = GAoptions.u;
    targetStructureOptions.v = GAoptions.v;
    targetStructureOptions.w = GAoptions.w;
    targetStructureOptions.radb = GAoptions.period/5; %radius / distance of middle of helix to axis (um)
    targetStructureOptions.radl = GAoptions.period/6; %radius of the wire (um)
    targetStructureOptions.relativeZ = GAoptions.C_over_A; %ratio of height to width of ellipsoidal 'pen'
    GAoptions.targetStructureOptions = targetStructureOptions;
    if strcmp(targetStructureOptions.structureType, 'doublehelix')
        GAoptions.targetStructure = generateDoubleHelixStructureParallelogram(targetStructureOptions)  %Create double helix
    elseif strcmp(targetStructureOptions.structureType, 'helix')
        GAoptions.targetStructure = generateHelixStructureParallelogram(targetStructureOptions)  %Create helix
    end
    targetfill = 1 - sum(sum(sum(GAoptions.targetStructure))) / (size(GAoptions.targetStructure,1)*size(GAoptions.targetStructure,1)*size(GAoptions.targetStructure,3));

    %This is a larger helix that counts as negative for filled spots outside of it
    exclusionStructureOptions.structureType = targetStructureOptions.structureType;
    exclusionStructureOptions.cells = GAoptions.cells;    
    exclusionStructureOptions.u = GAoptions.u;
    exclusionStructureOptions.v = GAoptions.v;
    exclusionStructureOptions.w = GAoptions.w;
    exclusionStructureOptions.radb = GAoptions.period/5;
    exclusionStructureOptions.radl = GAoptions.period/5;
    exclusionStructureOptions.relativeZ = GAoptions.C_over_A;  %ratio of height to width of ellipsoidal 'pen'
    GAoptions.exclusionStructureOptions = exclusionStructureOptions;
    %GAoptions.exclusionStructure = generateDoubleHelixStructureParallelogram(exclusionStructureOptions)  %Create double helix
    %GAoptions.exclusionStructure = generateHelixStructureParallelogram(exclusionStructureOptions)  %Helix
    if strcmp(targetStructureOptions.structureType, 'doublehelix')
        GAoptions.exclusionStructure = generateDoubleHelixStructureParallelogram(exclusionStructureOptions)  %Create double helix
    elseif strcmp(targetStructureOptions.structureType, 'helix')
        GAoptions.exclusionStructure = generateHelixStructureParallelogram(exclusionStructureOptions)  %Create helix
    end
    exclusionfill = 1 - sum(sum(sum(GAoptions.exclusionStructure))) / (size(GAoptions.exclusionStructure,1)*size(GAoptions.exclusionStructure,1)*size(GAoptions.exclusionStructure,3));
    
    %GAoptions.fill = (3*targetfill+exclusionfill)/4;  %Matches interference pattern fill to average of the target and exclusion structures. Consider adding fill factor to chromosome

    %1 at the edges of the unit cell to reduce fitness of x-y continuous structures
    edgeExclusionStructure = zeros(GAoptions.cells);
    edgeExclusionStructure(1,:,:) = 1;
    edgeExclusionStructure(end,:,:) = 1;
    edgeExclusionStructure(:,1,:) = 1;
    edgeExclusionStructure(:,end,:) = 1;
    GAoptions.edgeExclusionStructure = edgeExclusionStructure;
    
    
    
    %%%%% FILL %%%%%
    fillOptions.constFill = -1;
    fillOptions.fillMax = 0.85;
    fillOptions.fillMin = 0.8;
    fillOptions.chromNfill = 8;
    GAoptions.fill = fillOptions;
    GAoptions.fillHandler = FillFactorHandler(fillOptions)
    
    
    %%%%% SKIN %%%%%
    GAoptions.useSkinSingle = 1;  %Whether or not to do the skin calculation when doing a "single" calc
    GAoptions.skinIterations = 3; 
    
    
    
 
    %%%%% Offset %%%%%
    offsetOptions.offsetType = 'parallogram';
    offsetOptions.chromNoffsetX = 8;
    offsetOptions.chromNoffsetY = 8;
    offsetOptions.chromNoffsetZ = 8;
    offsetOptions.cells = GAoptions.cells;
    GAoptions.offsetConductor = OffsetConductor(offsetOptions);
    
    
    
    %%%% Calc Fitness %%%%%
    %GAoptions.calcFitness = @calcVolumetricMatchExclusion;
    


    
    
    
    
    
    
    %Writer for S4 runfile
    S4interfaceOptions.dimensions = GAoptions.dimensions;
    S4interfaceOptions.cells = GAoptions.cells;
    S4interfaceOptions.repeatingUnits = GAoptions.repeatingUnits;
    S4interfaceOptions.lattice = GAoptions.lattice;
    GAoptions.S4interface = S4interfaceSquareGeneral(S4interfaceOptions)
    GAoptions.S4interfaceOptions = S4interfaceOptions;
    
    %%%%% Do FDTD %%%%%
    fdtd.simulationTime = 500; %fs
    fdtd.shrinkFactor = 0.9; %Thickness after development / thickness during interference
    fdtd.repeatUnits = 3; %Number of unit cells to do FDTD simulation
    fdtd.meshAccuracy = 3;
    fdtd.minSourceWL = 0.5e-6; %m
    fdtd.maxSourceWL = 4e-6;  %m
    fdtd.minMeasWL = 0.5e-6; %m
    fdtd.maxMeasWL = 4e-6;  %m
    fdtd.numMeasWL = 500;
    fdtd.addCubesDirectly = 1;
    fdtd.SU8matrix = 1;
    fdtd.measureReflection = 1;
    GAoptions.fdtd = fdtd;
    
    
    if strcmp(GAoptions.fitnessType, 'fdtd')
        %%%%% Do FDTD as a test %%%%%
        fdtdfitness.dir = GAoptions.dir;
        fdtdfitness.dimensions = GAoptions.dimensions;
        fdtdfitness.cells = GAoptions.cells;
        fdtdfitness.lattice = GAoptions.lattice;
        fdtdfitness.baseLumSave = 'baseLum.fsp'; %base Lumerical file with common objects
        fdtdfitness.baseScriptFile = 'baseScript.lsf'; %Sets up the base file
        fdtdfitness.modLumSave = 'LumSave.fsp';  %File containing setup+structure
        fdtdfitness.runScriptFile = 'LumRun.lsf';  %Executed each sim
        fdtdfitness.resultsFile = 'LumResults';  %Location of output by Lumerical
        fdtdfitness.NKfile = 'LumNK.txt'  %Location of n,k structure file
        fdtdfitness.simulationTime = 250; %fs
        fdtdfitness.shrinkFactor = 0.9; %Thickness after development / thickness during interference
        fdtdfitness.repeatUnits = 3; %Number of unit cells to do FDTD simulation
        fdtdfitness.meshAccuracy = 1;
        fdtdfitness.minSourceWL = 0.7e-6; %m
        fdtdfitness.maxSourceWL = 2e-6;  %m
        fdtdfitness.minMeasWL = 0.7e-6; %m
        fdtdfitness.maxMeasWL = 2e-6;  %m
        fdtdfitness.numMeasWL = 100;
        fdtdfitness.addCubesDirectly = 0;
        fdtdfitness.SU8matrix = 1;
        fdtdfitness.measureReflection = 0;
        fdtdfitness.n_exposed = 1.6;    %for n,k fdtdstructure
        fdtdfitness.k_exposed = 0;
        fdtdfitness.n_inverted = 0.10858;  %Ag @ 1.33um  Babar and Weaver
        fdtdfitness.k_inverted = 9.6590;    % "
        GAoptions.lumInterfaceFitness = LumericalInterfaceFitness(fdtdfitness);
        GAoptions.fdtdfitness = fdtdfitness;
    end
    
    %%%%% Do sensitizer simulation %%%%%
    sensSim.sensDens = 8290000 *4; %Density of sensitizer molecules per um^3
    %abs_prob_per_photon = 2.91e-7; %nm^2/molecule %BCPI %based on 5-9-15 UV-vis measurements 
    sensSim.absCrossSection = 2.91e-7; %nm^2/molecule %BCPI  %MAY NEED TO UPDATE
    sensSim.Texposure = 200; %s
    GAoptions.sensSim = sensSim;
    
    
    
    
    %Options for built-in GA algorithm:
    if strcmp(strtrim(hostname),'lotus-bud')
        UseParallelVar = true;
    elseif strcmp(strtrim(hostname),'berzerk')
        if strcmp(GAoptions.fitnessType, 'fdtd')
            UseParallelVar = 'never'
        else
            UseParallelVar = 'always'
        end
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
    
    %Tells whether or not to do output of data/files etc.
    GAoptions.runSingle = 0;
    
    
    %Get layer chromosome length
    chromNlayer=0;
    for i=1:length(GAoptions.S4interfaceOptions.layers)
        chromNlayer=chromNlayer+GAoptions.S4interfaceOptions.layers(i).getChromosomeSize();
    end
    %Get material chromosome length
    chromNmaterial=0;
    for i=1:length(GAoptions.S4interfaceOptions.materials)
        chromNmaterial=chromNmaterial+GAoptions.S4interfaceOptions.materials(i).getChromosomeSize()
    end
    %Creation of 'GAproblem' for built-in Matlab GA
    GAoptions.GAproblem.options = GAoptions.optimset;  %options set above
    GAoptions.GAproblem.nvars = ...
        + GAoptions.gratingFunction.getChromosomeSize()  ...
        + GAoptions.incidentLightFunction.getChromosomeSize()  ...
        + GAoptions.offsetConductor.getChromosomeSize() ...
        + chromNlayer ...
        + chromNmaterial ...
        + GAoptions.fillHandler.getChromosomeSize() ...
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
GAoptions.runSingle = 1;
fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome) %Should also work for reflections




        
        
%%%%% Declaration of Fitness Function %%%%%
    function fitness = gfit(chromosome)  %This is done so that it can pass GAoptions to the fitness function.
        %tic
        fitness = fitnessFunction_8beamVarLayers(GAoptions,chromosome);
        %toc 
    end




end







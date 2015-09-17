%%MAKE SURE THIS IS UP TO DATE!!!


function runSingleChromosome(GAoptions,chromosome)

    messagestring = 'You are using the generic GA_writing version of runSingleChromosome';
    messagestring

    GAoptions.useRealisticInversion = false;
    GAoptions.addCubesDirectly = false;
    GAoptions.runSim = true;
    
    chromosome %output
    
    %change file names so data isn't overwritten
    GAoptions.GACreatorRecordFile = 'GACreatorRecord_single.mat';
    GAoptions.GARecordFileBase = 'GARecord_single';
    GAoptions.currentNKfile = 'currentNK_single.txt';
    GAoptions.baseLumSave = 'baseLumSave_single.fsp';
    GAoptions.currentLumSave = 'currentLumSave_single.fsp';
    GAoptions.currentLumResultsFile = 'currentLumResults_single.mat';
    GAoptions.LumSetupScript = 'LumSetup_single.lsf';
    GAoptions.LumRunScript = 'LumRun_single.lsf';
    GAoptions.fitnessFunctionExport = 'FitFuncExport_single.m';
    GAoptions.GACreatorExport = 'GACreatorExport_single.m';

    
    %change # of cells in simulation
    Ncellx = GAoptions.Ncellx*2
    Ncelly = GAoptions.Ncelly*2
    Ncellz = GAoptions.Ncellz*2  %number of cells used in the intensity map
    N = Ncellx*Ncelly*Ncellz;  %just want a single structure
    GAoptions.Ncellx = Ncellx;  GAoptions.Ncelly = Ncelly;  GAoptions.Ncellz = Ncellz; GAoptions.Ncell = N;
    
   
    
    %Change wavelengths
    GAoptions.minSourceWL = 500*1e-9;
    GAoptions.maxSourceWL = 4000*1e-9;
    GAoptions.minMeasWL = 500*1e-9; %1400*1e-9;
    GAoptions.maxMeasWL = 4000*1e-9; %1700*1e-9;
    GAoptions.numMeasWL = 300; %Number of wavelengths to get data from
    GAoptions.allMeasWL = GAoptions.minMeasWL:(GAoptions.minMeasWL-GAoptions.maxMeasWL)/GAoptions.numMeasWL:GAoptions.maxMeasWL;
    GAoptions.meshAccuracy = 2;
    
    writeLumericalSetupFile(GAoptions)
    system(['/opt/lumerical/fdtd/bin/fdtd-solutions -run ',GAoptions.dir,GAoptions.LumSetupScript]);
    %Write Lumerical Run File
    writeLumericalRunFile(GAoptions)
    
    
    chromosomalData = readChromosomalData(chromosome,GAoptions);
    
    
    %MAKE GRATING
    grating = makeGrating_addCircles(chromosomalData, GAoptions);
    
    intensityDist = calcIntensities(grating,chromosomalData,GAoptions);
    
    
    %Normalizes the threshold
    maxIntensity = max(real(reshape(intensityDist,GAoptions.Ncell,1)));
    minIntensity = min(real(reshape(intensityDist,GAoptions.Ncell,1)));
    threshold = chromosomalData.thresholdfraction*(maxIntensity-minIntensity)+minIntensity;
    
    %Convert the intensity distribution into filling data of SU8
    structureSU8 = intensityDist > threshold;  
    
    %Determine whether the structure is bicontinuous (returns 1) and if not, returns approx. how bicontinuous it is    
    [bicontinuousFactor, finalF, finalE] = ContinuityTest(structureSU8, GAoptions.useRealisticInversion);  
      
    %Calculate the bicontinuous factors of structures with different thresholds
    thresholdLow = (2*threshold+minIntensity)/3  %between actual and minimum-possible thresholds, but weighted towards actual
    thresholdHigh = (2*threshold+maxIntensity)/3  
    bicontinuousFactorLowThresh = ContinuityTest(intensityDist > thresholdLow, GAoptions.useRealisticInversion);
    bicontinuousFactorHighThresh = ContinuityTest(intensityDist > thresholdHigh, GAoptions.useRealisticInversion);
        
    
    %Save the structure to a file
    if GAoptions.hexagonalGrating
        if GAoptions.useRealisticInversion
            outputLumericalFile(finalE,[GAoptions.period,2*GAoptions.yperiod,GAoptions.Z_T_shrunk],GAoptions.n_Au, GAoptions.k_Au, 0, [GAoptions.dir,GAoptions.currentNKfile]); %0=don't invert
        else
            outputLumericalFile(structureSU8,[GAoptions.period,2*GAoptions.yperiod,GAoptions.Z_T_shrunk], GAoptions.n_Ag, GAoptions.k_Ag, 1, [GAoptions.dir,GAoptions.currentNKfile]); %1=invert
        end
    else %square grating
        if GAoptions.useRealisticInversion
            outputLumericalFile(finalE,[GAoptions.period,GAoptions.period,GAoptions.Z_T_shrunk],GAoptions.n_Au, GAoptions.k_Au, 0, [GAoptions.dir,GAoptions.currentNKfile]); %0=don't invert
        else
            outputLumericalFile(structureSU8,[GAoptions.period,GAoptions.period,GAoptions.Z_T_shrunk], GAoptions.n_Au, GAoptions.k_Au, 1, [GAoptions.dir,GAoptions.currentNKfile]); %1=invert
        end
    end
    
    
    %Remove old Lumerical output file (if it exists)
    if exist([GAoptions.dir,GAoptions.currentLumResultsFile], 'file')
        delete([GAoptions.dir, GAoptions.currentLumResultsFile]);
    end
    
    %Tell Lumerical to run the script
    system(['/opt/lumerical/fdtd/bin/fdtd-solutions -run ',GAoptions.dir,GAoptions.LumRunScript]);
    %Wait for the appearance of the Lumerical output file
    while (~exist([GAoptions.dir,GAoptions.currentLumResultsFile], 'file'))
        pause(0.1);
    end
    
    %Analyze the output data
    LumResults = load([GAoptions.dir,GAoptions.currentLumResultsFile]);
    transmission_right = LumResults.transmission_right;
    transmission_left = LumResults.transmission_left;
    trans_right = sum(transmission_right(1:end));
    trans_left = sum(transmission_left(1:end));
    trans_delta = abs(trans_right-trans_left)
    
    bicontinuousFitness = (3*bicontinuousFactor + bicontinuousFactorLowThresh + bicontinuousFactorHighThresh)^2
    fitness = -1 * trans_delta * bicontinuousFitness %More negative fitness = better fit

    
    
    %PLOT STRUCTURE;
    figure %creates a figure on the display
    %threshold=fixfill(IN,128, fillfraction);  %determines intensity threshold that will give the fill-fraction <- but threshold is given by chromosome
    patched = patch(isosurface(shave(structureSU8),threshold));  %creates surface at threshold
    %patched = patch(isosurface(shave(IN),threshold));  %creates surface at threshold
    %patched = patch(isosurface(IN,threshold));
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');  %display stuff
    view(3);
    axis tight
    axis equal 
    camlight
    lighting gouraud
    xlim([1 Ncelly])
    ylim([1 Ncellx])
    zlim([1 Ncellz])
    %xlim('auto')
    %ylim('auto')
    %zlim('auto')
    box on
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    
    %PLOT INVERSE STRUCTURE;
    figure %creates a figure on the display
    %threshold=fixfill(IN,128, fillfraction);  %determines intensity threshold that will give the fill-fraction <- but threshold is given by chromosome
    patched = patch(isosurface(shave_inverse(structureSU8),threshold));  %creates surface at threshold
    %patched = patch(isosurface(IN,threshold));
    set(patched,'FaceColor', [212 175 55]/256, 'EdgeColor', 'none');  %display stuff
    view(3);
    axis tight
    axis equal
    camlight
    lighting gouraud
    xlim([1 Ncelly])
    ylim([1 Ncellx])
    zlim([1 Ncellz])
    %xlim('auto')
    %ylim('auto')
    %zlim('auto')
    box on
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    
    
    
    frequencyStart = 2.99e8/GAoptions.maxMeasWL; %7.494e13; %0.5um
    frequencyStop = 2.99e8/GAoptions.minMeasWL;%5.996e14; %4um
    frequencies = frequencyStart:(frequencyStop-frequencyStart)/(GAoptions.numMeasWL-1):frequencyStop;
    wavelengths = 2.99e8./frequencies *1e6; %microns
    figure
    plot(wavelengths,transmission_left,wavelengths,transmission_right);
    
    
    
    
%     %CHANGE LUMERICAL RUN FILE SO IT DOESN'T CLOSE
%     clear oldLumRun;
%     oldLumRun = importdata([GAoptions.dir,GAoptions.LumRunScript], '\t', 400); 
%     for i=1:size(oldLumRun,1)  %replaces line that says 'exit' if it exists
%         oldLumRun(i,:);
%         if strcmp(oldLumRun(i,:),'exit; ')
%             oldLumRun(i,:) = cellstr(' ');
%         end
%         %if strcmp(oldLumRun(i,:),'runparallel; ')
%         %    oldLumRun(i,:) = cellstr(' ');
%         %end
%     end
%     newLumRun=char(oldLumRun);
%     LumRunOut = fopen([GAoptions.dir,LumRunScript], 'wt');
%     for l=1:size(newLumRun,1)
%         t=newLumRun(l,:);
%         fprintf(LumRunOut,' %s \n', t);
%     end
%     fclose (LumRunOut);
    
    %FIGURING OUT TALBOT REPEAT LENGTH  - should be done in GACreator
    %Talbot repeat length
     %d = GAoptions.period;
     %n = GAoptions.n_SU8;
     %lamda0 = GAoptions.incidentField.wavelength;
     %sqrt2 = sqrt(2);
     %GAoptions.verticalShrinkage=0;  %currently, assume no shrinkage
     %GAoptions.Z_T = (d/sqrt2)/(d*n/(lamda0*sqrt2)-sqrt((d*n/(lamda0*sqrt2))^2-1/2))  %SEEMS TO WORK - yeah
     %GAoptions.Z_T_shrunk = GAoptions.Z_T*(1-GAoptions.verticalShrinkage);
     %GAoptions.Ncellx = 30; GAoptions.Ncelly = 30; GAoptions.Ncellz = floor(GAoptions.Ncellx * GAoptions.Z_T/GAoptions.period);
   
    
    %%FROM FITNESS FUNCTION - WITH PLOTTING ADDED%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %Read chromosomal data
%     ucell = zeros(GAoptions.Ngratingx,GAoptions.Ngratingy);  %creates a 2D matrix for the mask
%     uc1 = reshape(chromosome(1:GAoptions.Ngrating),GAoptions.Ngratingx,GAoptions.Ngratingy,GAoptions.Nlev);  %reshapes first part of chromosome into a 3D matrix
%     for i = 1:GAoptions.Nlev  %for each level
%         ucell = ucell+uc1(:,:,i);   %increase height of 2D matrix by one if 
%     end  %basically the height of each square in ucell is equal to the sum of the filled blocks in that column
%     ucell %display ucell
%     u=GAoptions.Ngrating;  %u is the position in the chromosome that has already been read
%     thick = zeros(GAoptions.Nlev);
%     for i= 1:GAoptions.Nlev
%         thick(i) = convert_binary(chromosome(u+1:u+GAoptions.Nthick))/(2^GAoptions.Nthick-1);  %converts the binary strings signifying thickness into an integer array %fractional thickness
%         u=u+GAoptions.Nthick;
%     end
%     if GAoptions.useDielectric == 1  %only determine dielectric thickness if there is a dielectric
%         thick(GAoptions.Nlev+1) = convert_binary(chromosome(u+1:u+GAoptions.Ndl))/(2^GAoptions.Ndl-1);  %thickness of the dielectric
%         u=u+GAoptions.Ndl;
%     end
%     psi0 = convert_binary(chromosome(u+1:u+GAoptions.Npsi0))*2*pi/(2^GAoptions.Npsi0-1)  %convert psi
%     u=u+GAoptions.Npsi0;
%     chi0 = convert_binary(chromosome(u+1:u+GAoptions.Nchi0))*2*pi/(2^GAoptions.Nchi0-1) %convert chi
%     u=u+GAoptions.Nchi0;
%     thresholdfraction = convert_binary(chromosome(u+1:u+GAoptions.Nthresholdfraction))/(2^GAoptions.Nthresholdfraction-1)
%     u=u+GAoptions.Nthresholdfraction;
    
    
%     %FROM SID'S THESIS - FOR TESTING ONLY
%     %psi0 = 5.03; %from thesis
%     psi0 = 5.03 - 3/4*pi;
%     %psi0 = 3.00;
%     %psi0 = 3.20;
%     %psi0 = 6.15;
%     chi0 = 1.52;
%     %psi0 = 3.4456; %from "Best/spiral (ar=1.0)" <- actually a multilevelgrating
%     %chi0 = 3.4456;
%     %psi0 = 0.2027;
%     %chi0 = 0;
%     ucell =[1,1,1;0,1,0;0,0,0];
%     %thick(1) = 0.1613 %fraction
%     thick(1) = 0.199   %199nm
%     thresholdfraction = 0.45 %complete guess - not in thesis
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %MAKE GRATING
%     clear stratum stripe stratum1 stratum2 stripe1 stripe2 block1 block2
%     tmax = GAoptions.thickmax;%thickness maximum in terms of d/Nx
%     tmin = GAoptions.thickmin;
%     if GAoptions.useDielectric == 1
%         tmxdl = GAoptions.tmxdl;%thickness maximum in terms of d/Nx [ must be less than or equal to than tmin]
%         tmndl = GAoptions.tmndl;
%     end
%     d = GAoptions.period; %grating period (unitless)
%     Nlev = GAoptions.Nlev;
%     [Nx,Ny] = size(ucell);
%     
%     for q=1:Nlev
%         thick(q) = ((tmax-tmin)*thick(q)+tmin)*d/Nx  %converts the thickness fraction to a physical value 
%     end
%     
%     if GAoptions.useDielectric == 1
%         thick(Nlev+1) = ((tmxdl-tmndl)*thick(Nlev+1)+tmndl)*d/Nx;   %converts the thickness fraction of the dielectric to a physical value
%     end
%     
%     %makes it so that substrate is displayed with different color
%     GAoptions.n_SU8_substrate = GAoptions.n_SU8;
%     GAoptions.constantGrating.pmt={GAoptions.n_air^2,GAoptions.n_SU8^2,GAoptions.n_SU8_substrate};  %only so that substrate can be plotted with a different color
%     GAoptions.pmt_index_air = 1;
%     GAoptions.pmt_index_SU8 = 2;
%     GAoptions.pmt_index_SU8_substrate = 3;
%     GAoptions.constantGrating.pmt_sub_index=GAoptions.pmt_index_SU8_substrate;  %Index of the substrate in the above array
%     
%     grating = GAoptions.constantGrating;  %gets grating with parameters common to all gratings
    
%     if GAoptions.useDielectric == 1   %%%If there is a dielectric%%%
% %         qc = 1;%current structure-layer
% %         %SETTING UP STRATA
% %         stratum{2*Nlev+1}.type = 2;  %preallocates array
% %         for q=1:(2*Nlev)+1 %need two levels within each level because of dielectric
% %             stratum{q}.type=2;%biperiodic **
% %             stratum{q}.h11=1;%stratum periodicity is same as grating periodicity
% %             stratum{q}.h12=0;  %%CHANGE???
% %             stratum{q}.h21=0;
% %             stratum{q}.h22=1;    
% %             %stripe{q}.type = 1;%inhomogeneous  %%huh??  stripe{q} is never used
% %             if (rem(q,2))  %if q is odd
% %                 stratum{q}.thick = thick(Nlev+1);  %thickness of lay is equal to thickness of dielectric
% %             else  %if q is even
% %                 stratum{q}.thick = thick(qc) - thick(Nlev+1);  %thickness =  layer thickness - dielectric thickness
% %                 qc = qc + 1;  %number of structure-layers completed
% %             end      
% %         end
% % 
% %         %SETTING RELIEF PROFILE
% %         %1 = air; 2 = SU8; 3 = TiO2
% %         %lowest layer
% %         q = 1;
% %         ut = (ucell >= q);  %ut is equal to 1 where the corresponding ucell is atleast 1
% %         for i = 1:Nx
% %             stratum{q}.stripe{i}.type = 1;
% %             stratum{q}.stripe{i}.c1 = i/Nx;  %thickness of stripe
% %             for j = 1:Ny
% %                 stratum{q}.stripe{i}.block{j}.c2 = j/Ny;  %length of block
% %                 stratum{q}.stripe{i}.block{j}.pmt_index =  3 - 1*(ut(i,j));  %2 if filled, 3 if empty
% %             end
% %         end 
% % 
% % 
% %         for q = 1:Nlev;
% %             %upper half of layer q and lower half of q+1
% %             %qt = ones(Nx,Ny)*q;  %equal to layer number
% %             ut = (ucell >= q);  %equal to 1 if the height is at least height of current layer
% %             utm = (ucell >= (q-1));  %equal to 1 if the height is at least height of (current layer - 1)  %not used
% %             utp = (ucell >= (q+1));  %equal to 1 if the height is at least height of (current layer + 1)
% %             %Upper half of layer q
% %             for i = 1:Nx
% %                 stratum{(2*q)}.stripe{i}.type = 1;  %inhomogeneous stripe
% %                 stratum{(2*q)}.stripe{i}.c1 = i/Nx;  %width of stripe
% %                 for j = 1:Ny
% %                     stratum{(2*q)}.stripe{i}.block{j}.c2 = j/Ny;  %length of block
% %                     stratum{(2*q)}.stripe{i}.block{j}.pmt_index = 1 +(ut(i,j));  %2 if filled and 1 if empty
% %                 end
% %             end  
% %             %Lower half of (q+1)
% %             for i = 1:Nx
% %                 stratum{(2*q) +1}.stripe{i}.type = 1;
% %                 stratum{(2*q) +1}.stripe{i}.c1 = i/Nx;
% %                 for j = 1:Ny
% %                     stratum{(2*q) +1}.stripe{i}.block{j}.c2 = j/Ny;
% %                     stratum{(2*q) +1}.stripe{i}.block{j}.pmt_index =  0*(utm(i,j)) + 2*(ut(i,j)) -1*(utp(i,j)) + 1 ;  %3 if filled and above is empty;  2 if filled and above is filled; 1 if empty
% %                 end
% %             end  
% %         end
% % 
% %         for q = 1:(2*Nlev)+1  %copy the strata data to grating
% %             grating.stratum{q} = stratum{q};
% %         end
%         
%     else %%%%if there is no dielectric%%%
%         
%         stratum{Nlev}.type=2; %preallocates space for strata data
%         for q=1:(Nlev) %need two levels within each level because of dielectric
%             stratum{q}.type=2;%biperiodic **
%             stratum{q}.h11=0;%stratum periodicity is same as grating periodicity
%             stratum{q}.h12=1; 
%             stratum{q}.h21=1;
%             stratum{q}.h22=0;    
%             stratum{q}.thick = thick(q);  %thickness =  layer thickness
%         end   
%         for q = 1:Nlev;
%             for i = 1:Nx
%                 stratum{q}.stripe{i}.type = 1;  %inhomogeneous stripe
%                 stratum{q}.stripe{i}.c1 = i/Nx;  %width of stripe
%                 for j = 1:Ny
%                     stratum{q}.stripe{i}.block{j}.c2 = j/Ny;  %length of block
%                     %A block is filled at the nth level and below if there
%                     %are n '1's in that column of the chromosome.
%                     if ucell(i,j)>=q
%                         stratum{q}.stripe{i}.block{j}.pmt_index = GAoptions.pmt_index_SU8;
%                     else
%                         stratum{q}.stripe{i}.block{j}.pmt_index = GAoptions.pmt_index_air;
%                     end
%                 end
%             end
%             grating.stratum{q} = stratum{q};
%         end
%     end         
%     
%     
%     %PLOT GRATING:
%     %draw the grating
%     clear pmt_display;
%     pmt_display(1).name='Air';
%     pmt_display(1).color=[255/255,255/255,240/255];
%     pmt_display(1).alpha=0;
%     pmt_display(2).name = 'su-8';
%     pmt_display(2).color=[0/255,0/255,128/255];
%     pmt_display(2).alpha=1;
%     pmt_display(3).name = 'su-8_substrate';    %only for display
%     pmt_display(3).color=[255/255,255/255,255/255];
%     pmt_display(3).alpha=1;
%     if GAoptions.useDielectric
%         pmt_display(3).name = 'TiO2';
%         pmt_display(3).color=[127/255,127/255,127/255];
%         pmt_display(3).alpha=1;
%     end
% 
%     x_limit=[-0.5,0,0;0.5,2*GAoptions.period,2*GAoptions.yperiod];
%     h_plot=gdc_plot(grating,1,pmt_display,x_limit);
%     %if ~isempty(h_plot);
%     %    view(37.5,20);
%     %end
%     
%     
%     %CALCULATE TRANSMISSION/REFLECTION MATRICES
%     [param_size, scat_field, inc_field] = gdc(grating,GAoptions.incidentField,GAoptions.order);
% 
%     
%     %INCIDENT POLARIZATION
%     Esp = [cos(psi0)*cos(chi0)-1i*sin(psi0)*sin(chi0);...
%            sin(psi0)*cos(chi0)+1i*cos(psi0)*sin(chi0)];
%     
%     %CALCULATE PROPAGATING DIFFRACTED BEAMS
%     [E231, k231] = field_convert(Esp,scat_field);
%     
%     
%     %APPLY COORDINATE TRANSFORMATION - rotate the intensity matrix
%     %E = g.transform*E231;  %rotates it a set angle (currently zero)
%     %k = g.transform*k231;
% 
%     
%     offset = [0.1,0,0.3];
%     %N = icalc(E231,k231,GAoptions.period,Ncellx,Ncelly,Ncellz,offset,2,1); %offset is in microns
%     if GAoptions.hexagonalGrating
%         IN = icalc_mod(E231,k231,GAoptions.period,GAoptions.yperiod*2,GAoptions.Z_T,Ncellx,Ncelly,Ncellz,offset,1,1);  %offset is in microns  %don't want it repeated yet
%     else
%         IN = icalc_mod(E231,k231,GAoptions.period,GAoptions.period,GAoptions.Z_T,Ncellx,Ncelly,Ncellz,offset,1,1);  %offset is in microns  %don't want it repeated yet
%     end
%     %Normalizes the threshold
%     maxIN = max(real(reshape(IN,N,1)));
%     minIN = min(real(reshape(IN,N,1)));
%     maxThreshold = GAoptions.maxThreshold;
%     minThreshold = GAoptions.minThreshold;
%     thresholdfraction = thresholdfraction*(maxThreshold-minThreshold)+minThreshold;
%     threshold = thresholdfraction*(maxIN-minIN)+minIN;
%     
%     %Convert the intensity distribution into filling data of SU8
%     structureSU8 = IN > threshold; 
%     
    
    
    
    
    %determines whether the structure is bicontinuous (returns 1) and if not, returns approx. how bicontinuous it is    
%    [bicontinuousFactor, finalF, finalE] = ContinuityTest(structureSU8,useRealisticInversion);  
    
%     %PLOT STRUCTURE;
%     figure %creates a figure on the display
%     %threshold=fixfill(IN,128, fillfraction);  %determines intensity threshold that will give the fill-fraction <- but threshold is given by chromosome
%     patched = patch(isosurface(shave(IN),threshold));  %creates surface at threshold
%     %patched = patch(isosurface(IN,threshold));
%     set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');  %display stuff
%     view(3);
%     axis tight
%     axis equal
%     camlight
%     lighting gouraud
%     xlim([1 Ncelly])
%     ylim([1 Ncellx])
%     zlim([1 Ncellz])
%     %xlim('auto')
%     %ylim('auto')
%     %zlim('auto')
%     box on
%     set(gca,'XTick',[],'YTick',[],'ZTick',[])
%     
%     
%     %PLOT INVERSE STRUCTURE;
%     figure %creates a figure on the display
%     %threshold=fixfill(IN,128, fillfraction);  %determines intensity threshold that will give the fill-fraction <- but threshold is given by chromosome
%     patched = patch(isosurface(shave_inverse(IN),threshold));  %creates surface at threshold
%     %patched = patch(isosurface(IN,threshold));
%     set(patched,'FaceColor', [212 175 55]/256, 'EdgeColor', 'none');  %display stuff
%     view(3);
%     axis tight
%     axis equal
%     camlight
%     lighting gouraud
%     xlim([1 Ncelly])
%     ylim([1 Ncellx])
%     zlim([1 Ncellz])
%     %xlim('auto')
%     %ylim('auto')
%     %zlim('auto')
%     box on
%     set(gca,'XTick',[],'YTick',[],'ZTick',[])
%     
    
    
    
%     %NEED TO SAVE THE STRUCTURE TO A FILE
%     %Save the structure to a file
%     %%%%%%%Lumerical Run File%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%     runout = fopen([GAoptions.dir,GAoptions.LumRunScript], 'w');
%     fprintf(runout,'load("%s"); \n',[GAoptions.dir,GAoptions.baseLumSave]);
%     fprintf(runout,'switchtolayout; \n');
%     if addCubesDirectly
%         
%         %ASSUMES HEXAGONAL
%         dx = GAoptions.period/Ncellx*1e-6;
%         dy = 2*GAoptions.yperiod/Ncelly*1e-6;
%         dz = GAoptions.Z_T_shrunk/Ncellz*1e-6;
%         %initial rect that will be copied
%         fprintf(runout,'addrect; set("name","base"); set("x span",%e); set("y span",%e); set("z span",%e); \n',...
%             dx, dy, dz);
%         %fprintf(runout,'setmaterial("Au(Gold)-CRC");');
%         fprintf(runout,'set("material","Au (Gold) - CRC"); \n');
%         fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -GAoptions.period/2*1e-6 - dx/2, -2*GAoptions.yperiod/2*1e-6 - dy/2, -GAoptions.Z_T_shrunk/2*1e-6 - dz/2 - GAoptions.Z_T_shrunk*1e-6*(GAoptions.Nstructrepeat-1)/2);
%         for i = 1:Ncellx
%             for j = 1:Ncelly
%                 for k = 1:Ncellz
%                     if ~structureSU8(i,j,k)
%                         fprintf(runout,'select("base"); copy(%e,%e,%e); set("name","%s"); \n',i*dx,j*dy,k*dz, 'C');
%                     end
%                 end
%             end
%         end
%         fprintf(runout,'select("base"); delete; \n');
%         %make copies
%         fprintf(runout,'select("C");');
%         fprintf(runout,'for (n=2:%i){ \n',GAoptions.Nstructrepeat);
%         fprintf(runout,'copy(%e,%e,%e); \n',0,0,GAoptions.Z_T_shrunk*1e-6);
%         fprintf(runout,'} \n');
%         
%     else   %import n,k
%         fprintf(runout,'addimport; \n');
%         fprintf(runout,'importnk("%s","%s",%e,%e,%e,%i); \n',[GAoptions.dir,GAoptions.currentNKfile],...
%             'microns',0,0,-GAoptions.Z_T_shrunk*1e-6*(GAoptions.Nstructrepeat-1)/2,0);
%         fprintf(runout,'for (n=2:%i){ \n',GAoptions.Nstructrepeat);
%             fprintf(runout,'copy(%e,%e,%e); \n',0,0,GAoptions.Z_T_shrunk*1e-6);
%         fprintf(runout,'} \n');
%     end
%     fprintf(runout,'save("%s"); \n',[GAoptions.dir,GAoptions.currentLumSave]);
% %fprintf(runout,'runparallel; \n');
%     fprintf(runout,'transmission_left = transmission("transmission"); \n');
%     fprintf(runout,'reflection_left = transmission("reflection"); \n');
% 
%     fprintf(runout,'switchtolayout; \n');
%     fprintf(runout,'select("source2"); \n');
%     fprintf(runout,'set("phase",-90); \n');  %Set to right-handed polarization
% %fprintf(runout,'runparallel; \n');
%     fprintf(runout,'transmission_right = transmission("transmission"); \n');
%     fprintf(runout,'reflection_right = transmission("reflection"); \n');
% 
%     fprintf(runout,'matlabsave("%s",transmission_right,reflection_right,transmission_left,reflection_left); \n',[GAoptions.dir,GAoptions.currentLumResultsFile]);
% 
%     %fprintf(runout,'exit; \n');
%     fclose(runout);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if ~addCubesDirectly
%     if GAoptions.hexagonalGrating
%         if GAoptions.useRealisticInversion
%             outputLumericalFile(finalE,[GAoptions.period,2*GAoptions.yperiod,GAoptions.Z_T_shrunk],GAoptions.n_Au, GAoptions.k_Au, 0, [GAoptions.dir,GAoptions.currentNKfile]); %0=don't invert
%         else
%             outputLumericalFile(structureSU8,[GAoptions.period,2*GAoptions.yperiod,GAoptions.Z_T_shrunk], GAoptions.n_Au, GAoptions.k_Au, 1, [GAoptions.dir,GAoptions.currentNKfile]); %1=invert
%         end
%     else %square grating
%         if GAoptions.useRealisticInversion
%             outputLumericalFile(finalE,[GAoptions.period,GAoptions.period,GAoptions.Z_T_shrunk],GAoptions.n_Au, GAoptions.k_Au, 0, [GAoptions.dir,GAoptions.currentNKfile]); %0=don't invert
%         else
%             outputLumericalFile(structureSU8,[GAoptions.period,GAoptions.period,GAoptions.Z_T_shrunk], GAoptions.n_Au, GAoptions.k_Au, 1, [GAoptions.dir,GAoptions.currentNKfile]); %1=invert
%         end
%     end
% end
% 
%     %Remove old Lumerical output file
%     if exist([GAoptions.dir,GAoptions.currentLumResultsFile], 'file')
%         delete([GAoptions.dir, GAoptions.currentLumResultsFile]);
%     end
%     %NEED TO GET LUMERICAL TO RUN THE SCRIPT
%     system(['fdtd-solutions -run ',GAoptions.dir,GAoptions.LumRunScript]);
%     %NEED TO WAIT FOR THE APPEARANCE OF THE LUMERICAL OUTPUT FILE
%     while (~exist([GAoptions.dir,GAoptions.currentLumResultsFile], 'file'))         %%UNCOMMENT!!!!   
%         pause(0.1);
%     end
%     %NEED TO ANALYZE THE OUTPUT FILE
%     LumResults = load([GAoptions.dir,GAoptions.currentLumResultsFile]);
%     transmission_right = LumResults.transmission_right;
%     transmission_left = LumResults.transmission_left;
%     trans_right = sum(transmission_right(1:end));
%     trans_left = sum(transmission_left(1:end));
%     trans_delta = abs(trans_right-trans_left)
%        
%     %Determining wavelengths:
%     frequencyStart = 7.494e13; %0.5um
%     frequencyStop = 5.996e14; %4um
%     frequencyNum = 300;
%     frequencies = frequencyStart:(frequencyStop-frequencyStart)/(frequencyNum-1):frequencyStop;
%     wavelengths = 2.99e8./frequencies *1e6; %microns
%     
%     figure
%     plot(wavelengths,transmission_left,wavelengths,transmission_right);
%     
%     
%     fitness = -1 * trans_delta * bicontinuousFactor;



end
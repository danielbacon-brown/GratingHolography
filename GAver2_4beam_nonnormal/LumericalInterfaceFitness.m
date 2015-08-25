classdef LumericalInterfaceFitness
    %This is intended for using Lumerical as a part of the fitness function
    %Need to limit the amount of modifications made in lumerical to improve
    %efficiency.
    
    properties
%         dimensions;
%         cells;
%         SU8matrix; %boolean, whether to have metal lattice in a matrix
%         repeatUnits;
%         shrinkFactor;
%         minMeasWL;
%         maxMeasWL;
%         addCubesDirectly;
%         meshAccuracy;
%         simulationTime; %fs
%         
%         n_exposed;
%         k_exposed;
%         n_inverted;
%         k_inverted;
%         
%         dir;
%         baseScriptFile;
%         baseLumSave;
%         modLumSave;
%         runScriptFile;
%         resultsFile;
%         NKfile;
        
    options;
    structureThickness;
    end
    
    methods
        
        function L = LumericalInterfaceFitness(inoptions)
%             L.dimensions = options.dimensions;
%             L.cells = options.cells;
%             L.SU8matrix = options.SU8matrix;
%             L.repeatUnits = options.repeatUnits;
%             L.shrinkFactor = options.shrinkFactor;
%             L.minMeasWL = options.minMeasWL;
%             L.maxMeasWL = options.maxMeasWL;
%             L.addCubesDirectly = options.addCubesDirectly;
%             L.meshAccuracy = options.meshAccuracy;
%             L.simulationTime = options.simulationTime;
%             
%             L.n_exposed = options.n_exposed;
%             L.k_exposed = options.k_exposed;
%             L.n_inverted = options.n_inverted;
%             L.k_inverted = options.k_inverted;
%             
%             L.dir = options.dir;
%             L.baseScriptFile = options.baseScriptFile;
%             L.baseLumSave = options.baseLumSave;
%             L.modLumSave = options.modLumSave;
%             L.runScriptFile = options.runScriptFile;
%             L.NKfile = options.NKfile;
%             L.resultsFile = options.resultsFile;
            L.options = inoptions;
            
            
            L.options.dimensions(2) = L.options.dimensions(2) 
            
            setupLumerical(L)
        end
        
        
        function setupLumerical(L)
            %Assume that if the structure is hexagonal, the structure includes 2
            %periods in Y.
            
            setupout = fopen([L.options.dir,L.options.baseScriptFile], 'w');
            
            
%             Nx = L.cells(1);
%             Ny = L.cells(2);
%             Nz = L.cells(3);

            
            
            %Create repeat structure
            %structureSU8 = repmat(structureSU8,[1,1,L.options.fdtd.repeatUnits]);
            
            
            
            
            %SETUP
            fprintf(setupout,'newproject; \n');
            
            L.structureThickness = L.options.shrinkFactor*L.options.repeatUnits*L.options.dimensions(3);
            
            %Mesh
            fprintf(setupout,'addfdtd; \n');
            fprintf(setupout,L.setXYdimensionStr() );
            fprintf(setupout,'set("z", 0); \n');
            fprintf(setupout,'set("simulation time",%i); \n', L.options.simulationTime*1e-15);
            fprintf(setupout,'set("z span", %e ); \n', L.structureThickness*3*1e-6);
            fprintf(setupout,'set("mesh type","auto non-uniform"); \n');
            fprintf(setupout,'set("mesh accuracy",%i); \n',L.options.meshAccuracy);
            fprintf(setupout,'set("use early shutoff",0); \n');
            fprintf(setupout,'set("x min bc","periodic"); \n');
            fprintf(setupout,'set("y min bc","periodic"); \n');
            
            %Plane Source
            fprintf(setupout,'addplane; \n');
            fprintf(setupout,'set("name","source1"); \n');
            fprintf(setupout,L.setXYdimensionStr() );
            fprintf(setupout,'set("z", %e); \n', (-L.structureThickness/2-0.3)*1e-6);
            fprintf(setupout,'set("wavelength start",%e); \n', L.options.minSourceWL);
            fprintf(setupout,'set("wavelength stop", %e); \n', L.options.maxSourceWL);
            
            %Plane Source for Circular Polarization  %Right-handed (I think)
            fprintf(setupout,'copy(0,0,0); \n');
            fprintf(setupout,'set("name","source2"); \n');
            fprintf(setupout,'set("polarization angle",90); \n');
            fprintf(setupout,'set("phase",90); \n');
            
            %Transmission detector
            fprintf(setupout,'addprofile; \n');
            fprintf(setupout,'set("name","transmission"); \n');
            fprintf(setupout,L.setXYdimensionStr() );
            fprintf(setupout,'set("z", %e); \n', (L.structureThickness/2+0.5)*1e-6); %transmission detector is 2 microns above structure
            fprintf(setupout,'set("override global monitor settings",1); \n');
            fprintf(setupout,'set("use source limits",0); \n');
            fprintf(setupout,'set("minimum wavelength",%e); \n',L.options.minMeasWL);
            fprintf(setupout,'set("maximum wavelength",%e); \n',L.options.maxMeasWL);
            fprintf(setupout,'set("frequency points",%i); \n',L.options.numMeasWL);
            fprintf(setupout,'set("use linear wavelength spacing",0); \n');
            
%             %Reflection detector
%             fprintf(setupout,'addprofile; \n');
%             fprintf(setupout,'set("name","reflection"); \n');
%             setXYdimensionStr(setupout,L.options);
%             fprintf(setupout,'set("z", %e); \n', (-L.structureThickness/2-0.5)*1e-6); %reflection detector is 2 microns below structure
%             fprintf(setupout,'set("override global monitor settings",1); \n');
%             fprintf(setupout,'set("use source limits",0); \n');
%             fprintf(setupout,'set("minimum wavelength",%e); \n',L.options.minMeasWL);
%             fprintf(setupout,'set("maximum wavelength",%e); \n',L.options.maxMeasWL);
%             fprintf(setupout,'set("frequency points",%i); \n', L.options.numMeasWL);
%             fprintf(setupout,'set("use linear wavelength spacing",0); \n');
            
            %Add background material (matrix)
            if L.options.SU8matrix && L.options.addCubesDirectly
                fprintf(setupout,'addrect; \n');
                fprintf(setupout,L.setXYdimensionStr() );
                fprintf(setupout,'set("z",%e); \n', 0);
                fprintf(setupout,'set("z span",%e); \n', L.structureThickness*1e-6);
                fprintf(setupout,'set("material","<Object defined dielectric>"); \n');
                fprintf(setupout,'set("index",1.58); \n');
                fprintf(setupout,'set("override mesh order from material database", 1); \n');
                fprintf(setupout,'set("mesh order",3); \n');
                
            end
            
            %Save the base file
            fprintf(setupout,'save("%s"); \n',[L.options.dir,L.options.baseLumSave]);
            fprintf(setupout,'exit; \n');
            
            fclose(setupout);
            
            %Run setup
            system(sprintf('fdtd-solutions -run %s',  [L.options.dir, L.options.baseScriptFile] ) );
            
            
            
            %Can create a common run file if an external structure file is
            %used
            if ~L.options.addCubesDirectly
                
                runout = fopen([L.options.dir,L.options.runScriptFile], 'w');
                
                %fprintf(runout,'open("%s"); \n',[L.options.dir,L.baseLumSave]);
               
                
                %Dimensions of the structure for import
                Dx = L.options.dimensions(1)*1e-6;
                if strcmp(L.options.lattice,'square')
                    Dy = L.options.dimensions(2)*1e-6;
                elseif strcmp(L.options.lattice,'hexagonal')
                    Dy = L.options.dimensions(2)*2*1e-6;
                end
                Dz = L.options.dimensions(3)*L.options.shrinkFactor*1e-6;
               
                %import n k material
                fprintf(runout,'addimport; \n');
                fprintf(runout,'importnk( "%s" , "microns", %f, %f, %f, 0); \n ', [L.options.dir,L.options.NKfile], Dx, Dy, Dz );
                
                %Need to change filename to avoid overwriting file
                fprintf(runout,'save("%s"); \n',[L.options.dir,L.options.modLumSave]);
                
                fprintf(runout,'runparallel; \n');

                
                fprintf(runout,'transmission_left = transmission("transmission"); \n');
                %fprintf(runout,'reflection_left = transmission("reflection"); \n');
                
                fprintf(runout,'switchtolayout; \n');
                fprintf(runout,'select("source2"); \n');
                fprintf(runout,'set("phase",-90); \n');  %Set to right-handed polarization
                %if L.options.runSim
                fprintf(runout,'runparallel; \n');
                %end
                fprintf(runout,'transmission_right = transmission("transmission"); \n');
                %fprintf(runout,'reflection_right = transmission("reflection"); \n');
                
                fprintf(runout,'matlabsavelegacy("%s",transmission_right,transmission_left); \n',[L.options.dir,L.options.resultsFile]);
                
                fprintf(runout,'exit; \n');
                fclose(runout);
                
            end
            
        end
        
        
        function writeLumericalRunFileCubes(L, exposedStruct )
            
            
            
            %outputLumericalFile(exposedStruct, L.dimensions, L.n_exposed, L.k_exposed, L.n_inverse, L.k_inverse, 0, [L.dir,L.currentNKfile] );
            
            
            Nx = size(exposedStruct,1);
            Ny = size(exposedStruct,2);
            Nz = size(exposedStruct,3);
            
            %Create repeat structure in z
            exposedStructRep = repmat(exposedStruct,[1,1,L.options.repeatUnits]);
            
            fprintf(setupout,'redrawoff; \n');
            
            
            dx = L.options.dimensions(1)/Nx*1e-6;
            if strcmp(L.options.lattice,'square')
                dy = L.options.dimensions(2)/Ny*1e-6;
            elseif strcmp(L.options.lattice,'hexagonal')
                dy = L.options.dimensions(2)*2/Ny*1e-6;
            end
            dz = L.options.dimensions(3)*L.options.shrinkFactor/Nz*1e-6;
            %initial rect that will be copied
            fprintf(setupout,'addrect; set("name","base"); set("x span",%e); set("y span",%e); set("z span",%e); \n',...
                dx, dy, dz);
            %fprintf(setupout,'set("material","Au (Gold) - CRC"); \n');
            fprintf(setupout,'set("material","Ag (Silver) - CRC"); \n');
            
            if strcmp(L.options.lattice,'square')
                fprintf(setupout,'set("x",%e); set("y",%e); set("z",%e); \n', -L.options.dimensions(1)/2*1e-6 - dx/2, -L.options.dimensions(2)/2*1e-6 - dy/2, -L.options.structureThickness/2*1e-6 - dz/2);
            elseif strcmp(L.options.lattice,'hexagonal')
                fprintf(setupout,'set("x",%e); set("y",%e); set("z",%e); \n', -L.options.dimensions(1)/2*1e-6 - dx/2, -L.options.dimensions(2)*1e-6 - dy/2, -L.options.structureThickness/2*1e-6 - dz/2);
                
            end
            %fprintf(setupout,'set("x",%e); set("y",%e); set("z",%e); \n', -L.options.dimensions(1)/2*1e-6 - dx/2, -2*L.options.dimensions(2)/2*1e-6 - dy/2, -param.Z_T_shrunk/2*1e-6 - dz/2 - param.Z_T_shrunk*1e-6*(L.options.Nstructrepeat-1)/2);
            for i = 1:Nx
                for j = 1:Ny
                    for k = 1:Nz*L.options.repeatUnits
                        if ~exposedStructRep(i,j,k)
                            fprintf(setupout,'select("base"); copy(%e,%e,%e); set("name","%s"); \n',i*dx,j*dy,k*dz, 'C');
                        end
                    end
                end
            end
            fprintf(setupout,'select("base"); delete; \n');
            %make copies
            %fprintf(setupout,'select("C");');
            %fprintf(setupout,'for (n=2:%i){ \n',L.options.Nstructrepeat);
            %fprintf(setupout,'copy(%e,%e,%e); \n',0,0,L.options.Z_T_shrunk*1e-6);
            %fprintf(setupout,'} \n');
                
%             else %Import n,k 
% 
%                 dimLum = L.dimensions;
%                 dimLum(3) = L.shrinkFactor*L.options.repeatUnits*L.dimensions(3);
%                 %n_SU8 = 1.6;
%                 %k_SU8 = 0;
%                 %n_Ag = 0.10858;  %1.33um  Babar and Weaver
%                 %k_g = 9.6590;    % "
%                 outputLumericalFile( exposedStructRep, dimLum, L.n_exposed,L.k_exposed,L.n_inverted,L.k_inverted, 0, [L.dir,L.LumNKfile] );
%                 
%                 sprintf('addimport; \n');
%                 sprintf('importnk("%s","%s",%e,%e,%e,%i); \n',[L.dir,L.currentNKfile],...
%                     'microns',0,0,0*1e-6,0);
%                 
%             end
%         
        
        
        end
        
        
               
        
        function fitness = calcCDfitness(L,exposedStruct)
            
            
            if L.options.addCubesDirectly
                L.writeLumericalRunFile(exposedStruct);
            else
                outputLumericalFile(exposedStruct, L.options.dimensions, L.options.n_exposed, L.options.k_exposed, L.options.n_inverted, L.options.k_inverted, 0, [L.options.dir,L.options.NKfile] );
            end
            
            
            %Load base file and run sim
            system(sprintf('fdtd-solutions %s -run %s', [L.options.dir, L.options.baseLumSave] ,  [L.options.dir, L.options.runScriptFile] ) );
            %Import results
            while(~exist([L.options.dir,L.options.resultsFile,'.mat'],'file'))
                pause(0.1)
            end
            LumResults = load([L.options.dir,L.options.resultsFile]);
            %Calculate CD
            transmissionRight = sum(LumResults.transmission_right)/sqrt(2);
            transmissionLeft = sum(LumResults.transmission_left)/sqrt(2);
            fitness = -1*abs( transmissionRight - transmissionLeft );    
            
        end
        
            
%         function writeLumericalRunFileSquare( L.options, structureSU8)
%             %This writes the Lumerical Script that imports the structure, runs the
%             %simulation, and exports the data.
%             
%             
%             setXYdimensionStr = setXYdimensionStr();
%             setWLlimStr = setWavelengthLimStr();
%             
%             %Create repeat structure
%             structureSU8 = repmat(structureSU8,[1,1,L.options.repeatUnits]);
%             
%             
%             %SETUP
%             sprintf('newproject; \n');
%             
%             structureThickness = L.options.fdtd.shrinkFactor*L.options.fdtd.repeatUnits*L.options.dimensions(3);
%             
%             %Mesh
%             setupStr = [ setupStr, '\n' ...
%                 , 'addfdtd; \n' ...
%                 , setXYstr ...
%                 , sprintf('set("z", 0); \n') ...
%                 , sprintf('set("simulation time",%i); \n', L.options.fdtd.simulationTime*1e-15) ...
%                 , sprintf('set("z span", %e ); \n', (structureThickness*2+2)*1e-6) ...
%                 , sprintf('set("mesh type","auto non-uniform"); \n') ...
%                 , sprintf('set("mesh accuracy",%i); \n',L.options.fdtd.meshAccuracy) ...
%                 , sprintf('set("use early shutoff",0); \n') ...
%                 , sprintf('set("x min bc","periodic"); \n') ...
%                 , sprintf('set("y min bc","periodic"); \n') ...
%                 ];
%             
%             %Plane Source
%             setupStr = [ setupStr, '\n' ...
%                 'addplane; \n' ...
%                 , sprintf('set("name","source1"); \n') ...
%                 , setXYstr ...
%                 , sprintf('set("z", %e); \n', (-structureThickness/2-0.3)*1e-6) ...
%                 , setWLlimStr ...
%                 ];
%             
%             %Plane Source for Circular Polarization
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('copy(0,0,0); \n') ...
%                 , sprintf('set("name","source2"); \n') ...
%                 , sprintf('set("polarization angle",90); \n') ...
%                 , sprintf('set("phase",90); \n') ...
%                 ];
%             
%             %Transmission detector
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('addprofile; \n') ...
%                 , sprintf('set("name","transmission"); \n') ...
%                 , setXYstr ...
%                 , sprintf('set("z", %e); \n', (structureThickness/2+0.5)*1e-6) ... %transmission detector is 2 microns above structure
%                 , sprintf('set("override global monitor settings",1); \n') ...
%                 , sprintf('set("use source limits",0); \n') ...
%                 , setWLlimStr ...
%                 , sprintf('set("frequency points",%i); \n',L.options.fdtd.numMeasWL) ...
%                 , sprintf('set("use linear wavelength spacing",0); \n') ...
%                 ];
%             
%             %Reflection detector
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('addprofile; \n') ...
%                 , sprintf('set("name","reflection"); \n') ...
%                 , setXYstr ...
%                 , sprintf('set("z", %e); \n', (-structureThickness/2-0.5)*1e-6) ... %reflection detector is 2 microns below structure
%                 , sprintf('set("override global monitor settings",1); \n') ...
%                 , sprintf('set("use source limits",0); \n') ...
%                 , setWLlimStr ...
%                 , sprintf('set("frequency points",%i); \n', L.options.fdtd.numMeasWL) ...
%                 , sprintf('set("use linear wavelength spacing",0); \n') ...
%                 ];
%             
%  
%             if L.addSU8matrix
%                 setupStr = [ setupStr, '\n' ...
%                     , sprintf('addrect; \n')
%                     , sprintf('set("name","base"); \n') ...
%                     , setXYstr ...
%                     , sprintf('set("z span",%e); \n', structureThickness*1e-6) ];
%                 
%             end
%             
%             
%             if L.options.fdtd.addCubesDirectly
%                 
%                 setupStr = [ setupStr, '\n' ...
%                     sprintf('redrawoff; \n') ... %prevents redrawing of elements during the creation of the cubes
%                     ];
%                 
%                 
%                 dx = L.options.dimensions(1)/L.options.cells(1)*1e-6;
%                 dy = L.options.dimensions(2)/L.options.cells(2)*1e-6;
%                 dz = L.options.dimensions(3)*L.options.fdtd.shrinkFactor/L.options.cells(3)*1e-6;
%                 %initial rect that will be copied
%                 setupStr = [setupStr, sprintf('addrect; set("name","base"); set("x span",%e); set("y span",%e); set("z span",%e); \n', ...
%                     dx, dy, dz) ];
%                 
%                 switch L.metal
%                     case 'gold'
%                         setupStr = [setupStr, sprintf('set("material","Au (Gold) - CRC"); \n') ];
%                     case 'silver'
%                         setupStr = [setupStr, sprintf('set("material","Ag (Silver) - CRC"); \n') ];
%                 end
%                 
%                 setupStr = [setupStr, sprintf('set("x",%e); set("y",%e); set("z",%e); \n', -L.options.dimensions(1)/2*1e-6 - dx/2, -L.options.dimensions(2)/2*1e-6 - dy/2, -structureThickness/2*1e-6 - dz/2)];
%                 for i = 1:L.options.cells(1)
%                     for j = 1:L.options.cells(2)
%                         for k = 1:L.options.cells(3)*L.options.fdtd.repeatUnits
%                             if ~structureSU8(i,j,k)
%                                 setupStr = [setupStr, sprintf('select("base"); copy(%e,%e,%e); set("name","%s"); \n',i*dx,j*dy,k*dz, 'C')];
%                             end
%                         end
%                     end
%                 end
%                 setupStr = [setupStr, sprintf('select("base"); delete; \n')];
%                 
%                 
%             else
%                 dimLum = L.dimensions;
%                 dimLum(3) = L.options.fdtd.shrinkFactor*L.options.fdtd.repeatUnits*L.options.dimensions(3);
%                 n_SU8 = 1.6;
%                 k_SU8 = 0;
%                 n_Ag = 0.10858;  %1.33um  Babar and Weaver
%                 k_g = 9.6590;    % "
%                 outputLumericalFile( structureSU8, dimLum, n_SU8,k_SU8,n_Ag,k_Ag, 1, [L.options.dir,L.options.LumNKfile] );
%                 
%                 sprintf('addimport; \n');
%                 sprintf('importnk("%s","%s",%e,%e,%e,%i); \n',[L.options.dir,L.options.currentNKfile],...
%                     'microns',0,0,0*1e-6,0);
%                 
%                 %sprintf('for (n=2:%i){ \n',L.options.Nstructrepeat);
%                 %    sprintf('copy(%e,%e,%e); \n',0,0,param.Z_T_shrunk*1e-6);
%                 %sprintf('} \n');
%             end
%             
%             %Save and run
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('save("%s"); \n',[L.options.dir,L.options.currentLumSave]) ...
%                 , sprintf('runparallel; \n') ...
%                 ];
%             %Record LCP polarization data
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('transmission_left = transmission("transmission"); \n') ...
%                 , sprintf('reflection_left = transmission("reflection"); \n') ...
%                 ];
%             %Change to RCP polarization
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('switchtolayout; \n') ...
%                 , sprintf('select("source2"); \n') ...
%                 , sprintf('set("phase",-90); \n') ...  %Set to right-handed polarization
%                 ];
%             %Save and run again
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('save("%s"); \n',[L.options.dir,L.options.currentLumSave]) ...
%                 , sprintf('runparallel; \n') ...
%                 ];
%             %Record RCP polarization data
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('transmission_right = transmission("transmission"); \n') ...
%                 , sprintf('reflection_right = transmission("reflection"); \n') ...
%                 ];
%             
%             setupStr = [ setupStr, '\n' ...
%                 sprintf('matlabsavelegacy("%s",transmission_right,reflection_right,transmission_left,reflection_left); \n',[L.options.dir,L.options.currentLumResultsFile]) ...
%                 ];
%             
%             if L.options.fdtd.end
%                 setupStr = [setupStr, sprintf('exit; \n')];
%             end
%             
%             %Open file
%             setupout = fopen([L.options.dir,L.options.LumRunScript], 'w');
%             fprintf(setupout, setupStr);
%             fclose(setupout);
%             
%             
%         end
        
        
        function setXYstr = setXYdimensionStr(L)
            setXYstr = [ ...
                sprintf('set("x",%e); \n', 0) ...
                , sprintf('set("y",%e); \n', 0) ...
                , sprintf('set("x span",%e); \n', L.options.dimensions(1)*1e-6) ...
                , sprintf('set("y span",%e); \n', L.options.dimensions(2)*1e-6) ...
                ];
        end
        
        function setWLstr = setWavelengthLimStr(L)
            setWLstr = [ ...
                sprintf('set("minimum wavelength",%e); \n', L.options.minMeasWL) ...
                , sprintf('set("maximum wavelength",%e); \n', L.options.maxMeasWL) ...
                ];
        end
        
        
        
    end
    
    
end

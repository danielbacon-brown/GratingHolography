classdef MPBinterface
   
    properties
    hostname;
    n_fill;
    n_void;
    scriptFilename;
    resultFilename;
    end
    
    methods
        function M = MPBInterface(options)
            M.hostname = options.hostname;
            M.n_fill = options.n_fill;
            M.n_void = options.n_void;
            M.scriptFilename = 'mpbRun.ctl';
            M.resultFilename = 'thisResult.out';
            M.dataFilename = 'thisResult.dat';
        end
        
        
        function fitness = doMPB(M,GAoptions,structureFill)
            %simple tetragonal
            latticeStr = sprintf('(set! geometry-lattice (make lattice (basis-size %f %f %f)  (basis1 1 0 0) (basis2 0 1 0) (basis3 0 0 1))) \r\n', ...
                GAoptions.dimensions(1),GAoptions.dimensions(2),GAoptions.dimensions(3) );
            latticeStr = [latticeStr, '(set! resolution 25) \r\n', ...
                '(set! mesh-size 7) \r\n'];
            
            kVectorStr = sprintf( '(set! k-points (list (vector3 0 0 0) (vector3 0.5 0 0) (vector3 0 0.5 0) (vector3 0 0 0.5))) \r\n' ...
                );
            kVectorStr = [kVectorStr, '(set! k-points (interpolate 4 k-points)) \r\n'];
            kVectorStr = [kVectorStr, '(set! num-bands 4) \r\n']; 
%             (set! k-points (list (vector3 0 0 0)     ; Gamma
%                      (vector3 0.5 0 0)   ; X
%                      (vector3 0.5 0.5 0) ; M
%                      (vector3 0 0 0)))   ; Gamma
            %set the dielectric constant
            setDielectricStr = sprintf('(define diel (make dielectric (epsilon %s))) \r\n', ...
                M.n_fill^2);
            
            
            %Go through every point in the structure - create a block there
            %if filled:
            dX = GAoptions.dimensions(1)/GAoptions.cells(1);
            dY = GAoptions.dimensions(2)/GAoptions.cells(2);
            dZ = GAoptions.dimensions(3)/GAoptions.cells(3);
            addStructStr = '(set! geometry (list \r\n';
            for i = 1:size(structureFill,1)
                for j = 1:size(structureFill,2)
                    for k = 1:size(structureFill,3)
                        if ~structureFill(i,j,k)
                            addStructStr = [addStructStr, sprintf( '(make block (center %f %f %f) (size %f %f %f)  (material diel) ) \r\n ', ...
                                dX*i,dY*j,dZ*k,   dX,dY,dZ )];
                        end
                    end
                end
            end
            addStructStr = [addStructStr,') )\r\n'];
            
            runTEstr = '(run-te) \r\n';
            runTMstr = '(run-tm) \r\n';
            
            %retreiveBGstr = '(retrieve-gap 1)) \r\n';

            %Get frequencies:
            getFreqStr = ['grep freqs ',M.resultFilename,' > ',M.dataFilename,' \r\n'];
             
            
            %Write to file
            runMPB = fopen([GAoptions.dir,M.scriptFilename,'w');
            fprintf(runMPB, [latticeStr, kVectorStr, setDielectricStr, addStructStr, runTEstr, retreiveBGstr, runMPB]);
            fclose(runMPB);
            
            %Run script
            system(['mpb ', GAoptions.dir,M.scriptFilename, ' >& ', M.resultFilename,' \r\n']);
            
            %Import the data.  freqVal: [kindex,k1,k2,k3,kmag/2pi,band1,band2,band3m, ...] 
            fullData = importdata(M.scriptFilename);
            freqBand1 = fullData(:,6);
            freqBand2 = fullData(:,7);
            
            if GAoptions.runSingle == 1  %Do plotting for single test
                figure 
                plot(fullData(:,1),freqBand1,fullData(:,1),freqBand2)
            end
            
        end
        
    end
    
%CODE
% ; Dielectric spheres in a diamond (fcc) lattice.  This file is used in
% ; the "Data Analysis Tutorial" section of the MPB manual.
% 
% (set! geometry-lattice (make lattice 
% 			 (basis-size (sqrt 0.5) (sqrt 0.5) (sqrt 0.5))
% 			 (basis1 0 1 1)
% 			 (basis2 1 0 1)
% 			 (basis3 1 1 0)))
% 
% ; Corners of the irreducible Brillouin zone for the fcc lattice,
% ; in a canonical order:
% (set! k-points (interpolate 4 (list
% 			       (vector3 0 0.5 0.5)            ; X
% 			       (vector3 0 0.625 0.375)        ; U
% 			       (vector3 0 0.5 0)              ; L
% 			       (vector3 0 0 0)                ; Gamma
% 			       (vector3 0 0.5 0.5)            ; X
% 			       (vector3 0.25 0.75 0.5)        ; W
% 			       (vector3 0.375 0.75 0.375))))  ; K
% 
% ; define a couple of parameters (which we can set from the command-line)
% (define-param eps 11.56) ; the dielectric constant of the spheres
% (define-param r 0.25)    ; the radius of the spheres
% 
% (define diel (make dielectric (epsilon eps)))
% 
% ; A diamond lattice has two "atoms" per unit cell:
% (set! geometry (list (make sphere (center 0.125 0.125 0.125) (radius r) 
% 			   (material diel))
% 		     (make sphere (center -0.125 -0.125 -0.125) (radius r) 
% 			   (material diel))))
% 
% ; (A simple fcc lattice would have only one sphere/object at the origin.)
% 
% (set-param! resolution 16) ; use a 16x16x16 grid
% (set-param! mesh-size 5)
% (set-param! num-bands 5)
% 
% ; run calculation, outputting electric-field energy density at the U point:
% (run (output-at-kpoint (vector3 0 0.625 0.375) output-dpwr))


    
end
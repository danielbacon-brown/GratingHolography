classdef GratingGridHexagonalKnownHelix
    
    %Defines a grating according to a mesh, wherein the block in the mesh
    %can freely change in width and height.
    
    
    %This makes small variations on the typical helix mask
    
    %ASSUMES BINARY GRATING
    
    properties
            period; %Periodicity of grating %um
            NblockX; % Num blocks in x-dir
            NblockY; % Num blocks in y-dir
            constantGrating; %Beginning of initialization of grating (common to all gratings)
            pmtIndices; %index for permittivities (1=air, 2=grating/SU8)
            %chromNthickness; %Num chromosome positions needed for grating thickness
            chromNcellArr; %Num chromosome positions for on-off of grid
            chromNspacingX; %Num chromosome positions for x-lengths of the blocks
            chromNspacingY; %Num chromosome positions for y-lengths of the blocks
            chromNCX;
            chromNCY;
            %chromNSU8thickness; %thickness of the SU8 interference layer
            %thicknessMax; %Maximum value of chromosome thickness  %um
            spacingMin; %Minimum width of each block %um
            %SU8thicknessMin; %Limits of the thickness of the SU8 layer %um
            %SU8thicknessMax;
    end
   
    methods
        function G = GratingGridHexagonalKnownHelix(options)
            G.NblockX = options.NblockX;
            G.NblockY = options.NblockY;
            
            G.chromNspacingX = options.chromNspacingX;
            G.chromNspacingY = options.chromNspacingY;
            G.chromNcellArr = G.NblockX * G.NblockY;
            G.chromNCX = 8;
            G.chromNCY = 8;
            
            G.period = options.periodicity;   
            
            G.spacingMin = options.spacingMin;
            
            
            G.constantGrating.pmt={options.n_void^2, options.n_filled^2};  %Permitivity index
            G.pmtIndices.void = 1;
            G.pmtIndices.filled = 2;
            
            G.constantGrating.pmt_sub_index=G.pmtIndices.filled;  %Index of the substrate
            G.constantGrating.pmt_sup_index=G.pmtIndices.void;  %Index of the superstrate
            G.constantGrating.d21 = G.period;  %First periodicity vector in x-direction
            G.constantGrating.d31 = 0;  %First periodicity vector in y-direction
                
            %G.constantGrating.d22 = 0;  %Second periodicity vector x-direction
            %G.constantGrating.d32 = G.period;  %Second periodicity vector y-direction
            G.constantGrating.d22 = G.period/2;  %Second periodicity vector x-direction
            G.constantGrating.d32 = G.period*sqrt(3)/2;  %Second periodicity vector y-direction
            
                
            G.constantGrating.stratum={}; %container that will hold all the strata
            G.constantGrating.stratum{1}.type=2;%biperiodic
            G.constantGrating.stratum{1}.h11=0;%stratum periodicity is same as grating periodicity
            G.constantGrating.stratum{1}.h12=1;
            G.constantGrating.stratum{1}.h21=1;
            G.constantGrating.stratum{1}.h22=0; 
            
                        
            
            
        end
        
        function grating = generateGrating(G, chromosomeSection)
            
%             [spacingXchrom, spacingYchrom, cellArrChrom ] = splitChromosome(chromosomeSection,  ...
%                 [G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX-1), G.chromNcellArr]);
%             [spacingXchrom, spacingYchrom ] = splitChromosome(chromosomeSection,  ...
%                 [G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX)]);
            [spacingXchrom, spacingYchrom, CXchrom, CYchrom ] = splitChromosome(chromosomeSection,  ...
                [G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX), 4*G.chromNCX, 4*G.chromNCY]);
            
            cellArrChrom = [1,0;1,0;1,0;1,0]; %for typical helix pattern
            
            %Generate spacings of the grid
            A = G.NblockX;
            B = G.NblockY;
            spacingX = zeros(1,A+1);  % the block width as a fraction of the periodicity
            %spacingY = zeros(1,A+1);
            spacingY = zeros(1,A);  % one spacingY value for each stripe
            
            %CHANGE TO LIMIT THE SPACING TO WITHIN MINSPACINGX, etc
            %xspacing:
%             switch G.NblockX %This distributes the chromosome into 
%                 case 2
%                     [spacingXfrac(2)] = convertChrom_gc( spacingXchrom, G.chromNspacingX); 
%                 case 3
%                     [spacingXfrac(2),spacingXfrac(3)] = convertChrom_gc( spacingXchrom, G.chromNspacingX, G.chromNspacingX);
%                 case 4
%                     [spacingXfrac(2),spacingXfrac(3),spacingXfrac(4)] = convertChrom_gc( spacingXchrom, G.chromNspacingX, G.chromNspacingX, G.chromNspacingX)
%             end

            
            spacingXfrac = convertChrom_gc( spacingXchrom, ones(1,A-1)*G.chromNspacingX) ;  %Converts the chromosome data into an array of fractions 
            %spacingYfrac = convertChrom_gc( spacingYchrom, ones(1,B-1)*G.chromNspacingY) ;
            spacingYfrac = convertChrom_gc( spacingYchrom, ones(1,A)*G.chromNspacingY) ;
            CXfrac = convertChrom_gc(CXchrom, ones(1,4)*G.chromNCX);  %Used for chamfers and bevels
            CYfrac = convertChrom_gc(CYchrom, ones(1,4)*G.chromNCY);
            
            
            %Base values (will give previously used helix pattern)
            spacingXbase = [0.22185;0.2851;0.5];
            spacingYbase = [0.6391;0.6391;0.2717;0.2717];
            
%             spacingXmodMin = 0.9;
%             spacingXmodMax = 1.1;
%             spacingYmodMin = 0.9;
%             spacingYmodMax = 1.1;
            spacingXmodMin = 0.999;
            spacingXmodMax = 1.001;
            spacingYmodMin = 0.999;
            spacingYmodMax = 1.001;
            
            for ix = 1:size(spacingXbase,1)
                spacingXfrac(ix) = spacingXfrac(ix)*(spacingXmodMax-spacingXmodMin) + spacingXmodMin
            end
            for iy = 1:size(spacingYbase,1)
                spacingYfrac(iy) = spacingYfrac(iy)*(spacingYmodMax-spacingYmodMin) + spacingYmodMin;
            end
            
            
            %spacingX describes the positions of the edge between blocks
            remaining = 1; %describes thickness between most recent spacingX and the end of line
            for a = 1:(A-1)
                spacingX(a+1) = spacingXfrac(a)*spacingXbase(a)*( remaining) + (1-remaining)   %relative to origin
                %spacingX(a+1) = spacingXfrac(a)*( (remaining - G.spacingMin*(A-a)) - (G.spacingMin) ) + G.spacingMin+(1-remaining);   %frac*( maxT-minT)+minT
                remaining = 1-spacingX(a+1)
            end
            spacingX(A+1) = 1
            
            %remaining = 1; 
            for b = 1:A
                spacingY(b) = spacingYfrac(b) * spacingYbase(b);
                %spacingY(b+1) = spacingYfrac(b)*( (remaining - G.spacingMin*(B-b)) - (G.spacingMin) ) + G.spacingMin+(1-remaining);   %frac*( maxT-minT)+minT
                %remaining = 1-spacingY(b+1);
            end
            %spacingYout = spacingY;
            %spacingY(B+1) = 1

            %SPACING Y WORKS DIFFERENTLY FROM SPACINGX
            %already a good value, does not need scaling by remaining space
            %Just needs scaling by base values
            

            
            %Converts chrom into map of blocks to raise
            %map = reshape(cellArrChrom,A,B);
            %map = cellArrChrom;
            
            periodX = G.period
            periodY = G.period*sqrt(3)/2
            
            grating.W1 = spacingX(2)*periodX;
            grating.W2 = (spacingX(3)-spacingX(2)) * periodX;
            grating.W3 = (spacingX(4)-spacingX(3)) * periodX;
            grating.W4 = (spacingX(5)-spacingX(4)) * periodX;
            
            grating.L1 = spacingY(1)*periodY;
            grating.L2 = spacingY(2)*periodY;
            grating.L3 = spacingY(3)*periodY;
            grating.L4 = spacingY(4)*periodY;
            
            grating.CX1 = grating.W1*CXfrac(1);
            grating.CX2 = grating.W2*CXfrac(2);
            grating.CX3 = grating.W3*CXfrac(3);
            grating.CX4 = grating.W4*CXfrac(4);
            
            grating.CY1 = (grating.L1-grating.L4)*CYfrac(1);
            grating.CY2 = (grating.L2-grating.L3)*CYfrac(2);
            grating.CY3 = (grating.L2-grating.L3-grating.CY2)*CYfrac(3);
            grating.CY4 = (grating.L1-grating.L4-grating.CY1)*CYfrac(4)
            
            
            
            %With more complex masks, can no longer rely on this grating
            %model
%             %Make grating
%             grating = G.constantGrating;
%             
%             grating.stratum{1}.thick = 1; %thickness set by layer object
%             
%             
%             for ix = 1:G.NblockX %stripes run along y-direction
%                 grating.stratum{1}.stripe{ix}.type = 1;  %inhomogeneous stripe
%                 %grating.stratum{1}.stripe{j}.c1 = spacingYfrac(j+1);  %width of stripe
%                 grating.stratum{1}.stripe{ix}.c1 = spacingX(ix+1)
%                 %blockNum = 1; %the order of the current block
%                 for iy = 1:G.NblockY
% 
%                     %Only create a block right before it switches, to enhance the calculation speed
%                     %if  (iy == G.NblockX) | (map(iy,ix) ~= map(iy+1,ix))
%                     if iy == 1
%                         grating.stratum{1}.stripe{ix}.block{iy}.c2 = spacingY(ix); %ending position of block
%                     elseif iy == 2
%                         grating.stratum{1}.stripe{ix}.block{iy}.c2 = 1;
%                     end
%                         if map(ix,iy) == 1
%                             grating.stratum{1}.stripe{ix}.block{iy}.pmt_index = G.pmtIndices.filled;
%                         else
%                             grating.stratum{1}.stripe{ix}.block{iy}.pmt_index = G.pmtIndices.void;
%                         end
%                         %blockNum = blockNum+1;
%                     %end
%                 end
%             end
            
            
%             for j = 1:G.NblockY %stripes run along x-direction
%                 grating.stratum{1}.stripe{j}.type = 1;  %inhomogeneous stripe
%                 %grating.stratum{1}.stripe{j}.c1 = spacingYfrac(j+1);  %width of stripe
%                 grating.stratum{1}.stripe{j}.c1 = spacingY(j+1); 
%                 blockNum = 1; %the order of the current block
%                 for i = 1:G.NblockX
% 
%                     %A block is filled at the nth level and below if there
%                     %are n '1's in that column of the chromosome.
% 
% 
%                     %Only create a block right before it switches, to enhance the calculation speed
%                     if  (i == G.NblockX) | (map(i,j) ~= map(i+1,j))
%                         grating.stratum{1}.stripe{j}.block{blockNum}.c2 = spacingX(i+1); %ending position of block
%                         if map(i,j) == 1
%                             grating.stratum{1}.stripe{j}.block{blockNum}.pmt_index = G.pmtIndices.filled;
%                         else
%                             grating.stratum{1}.stripe{j}.block{blockNum}.pmt_index = G.pmtIndices.void;
%                         end
%                         blockNum = blockNum+1;
%                     end
%                 end
%             end
            
            
        end
        
        
        
        function chromosomeSize = getChromosomeSize(G)
            %chromosomeSize =  G.chromNspacingX*(G.NblockX-1) + G.chromNspacingY*(G.NblockY-1) + G.chromNcellArr;
            %chromosomeSize =  G.chromNspacingX*(G.NblockX-1) + G.chromNspacingY*(G.NblockX);
            chromosomeSize =  G.chromNspacingX*(G.NblockX-1) + G.chromNspacingY*(G.NblockX) + G.chromNCX*4 + G.chromNCY*4;
        end
        
        
        function plotGrating(G,grating)
            clear pmt_display   %Set colors
            pmt_display(1).name='void';
            pmt_display(1).color=[255/255,255/255,240/255];
            pmt_display(1).alpha=1;
            pmt_display(2).name ='filled';
            pmt_display(2).color=[0/255,0/255,128/255];
            pmt_display(2).alpha=1;
            
            x_limit=[-0.5,0,0;1.0,G.period*2,G.period*2];
            h_plot=gdc_plot(grating,1,pmt_display,x_limit);
            view(0,90);
        end
       
    end
    
end
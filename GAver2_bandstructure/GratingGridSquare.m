classdef GratingGridSquare
    %Defines a grating according to a mesh, wherein the block in the mesh
    %can freely change in width and height.
    
    
    %ASSUMES SQUARE
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
            %chromNSU8thickness; %thickness of the SU8 interference layer
            %thicknessMax; %Maximum value of chromosome thickness  %um
            spacingMin; %Minimum width of each block %um
            %SU8thicknessMin; %Limits of the thickness of the SU8 layer %um
            %SU8thicknessMax;
    end
   
    methods
        function G = GratingGridSquare(options)
            G.NblockX = options.NblockX;
            G.NblockY = options.NblockY;
            
            %G.chromNthickness = options.chromNthickness;
            G.chromNspacingX = options.chromNspacingX;
            G.chromNspacingY = options.chromNspacingY;
            G.chromNcellArr = G.NblockX * G.NblockY;
            %G.chromNSU8thickness = options.chromNSU8thickness;
            
            G.period = options.periodicity;   
            
            %G.thicknessMax = options.thicknessMax;
            G.spacingMin = options.spacingMin;
            
            %G.SU8thicknessMin = options.SU8thicknessMin;
            %G.SU8thicknessMax = options.SU8thicknessMax;
            
            G.constantGrating.pmt={options.n_void^2, options.n_filled^2};  %Permitivity index
            G.pmtIndices.void = 1;
            G.pmtIndices.filled = 2;
            
            G.constantGrating.pmt_sub_index=G.pmtIndices.filled;  %Index of the substrate
            G.constantGrating.pmt_sup_index=G.pmtIndices.void;  %Index of the superstrate
            G.constantGrating.d21 = G.period;  %First periodicity vector in x-direction
            G.constantGrating.d31 = 0;  %First periodicity vector in y-direction
                
            G.constantGrating.d22 = 0;  %Second periodicity vector x-direction
            G.constantGrating.d32 = G.period;  %Second periodicity vector y-direction
            
                
            G.constantGrating.stratum={}; %container that will hold all the strata
            G.constantGrating.stratum{1}.type=2;%biperiodic
            G.constantGrating.stratum{1}.h11=0;%stratum periodicity is same as grating periodicity
            G.constantGrating.stratum{1}.h12=1;
            G.constantGrating.stratum{1}.h21=1;
            G.constantGrating.stratum{1}.h22=0; 
            
                        
            
            
        end
        
        function grating = generateGrating(G, chromosomeSection)
            
%             thicknessChrom = chromosomeSection(1:G.chromNthickness);
%             remainingChrom = chromosomeSection( (G.chromNthickness+1):end);
            
%             [thicknessChrom, spacingXchrom, spacingYchrom, cellArrChrom, SU8thicknessChrom ] = splitChromosome(chromosomeSection,  ...
%                 [G.chromNthickness, G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX-1), G.chromNcellArr, G.chromNSU8thickness]);
%            [thicknessChrom, spacingXchrom, spacingYchrom, cellArrChrom ] = splitChromosome(chromosomeSection,  ...
%                [G.chromNthickness, G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX-1), G.chromNcellArr]);
            [spacingXchrom, spacingYchrom, cellArrChrom ] = splitChromosome(chromosomeSection,  ...
                [G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX-1), G.chromNcellArr]);

            %Generate spacings of the grid
            A = G.NblockX;
            B = G.NblockY;
            %spacingXfrac = zeros(1,A-1); % the block as a fraction of the range of possible thicknesses
            %spacingYfrac = zeros(1,B-1);
            spacingX = zeros(1,A+1);  % the block width as a fraction of the periodicity
            spacingY = zeros(1,A+1);
            
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

            %TEST
            %ones(1,A-1)*G.chromNspacingX
            %[vartest] = {convertChrom_gc( spacingXchrom, ones(1,A-1)*G.chromNspacingX)}
            
            spacingXfrac = convertChrom_gc( spacingXchrom, ones(1,A-1)*G.chromNspacingX) ;  %Converts the chromosome data into an array of fractions 
            spacingYfrac = convertChrom_gc( spacingYchrom, ones(1,B-1)*G.chromNspacingY) ;
            
            %spacingX describes the positions of the edge between blocks
            remaining = 1; %describes thickness between most recent spacingX and the end of line
            for a = 1:(A-1)
                %spacingX(a+1) = spacingXfrac(a)*( (a+1)/A *remaining) + (1-remaining);  %relative to origin
                spacingX(a+1) = spacingXfrac(a)*( (remaining - G.spacingMin*(A-a)) - (G.spacingMin) ) + G.spacingMin+(1-remaining);   %frac*( maxT-minT)+minT
                remaining = 1-spacingX(a+1);
            end
            spacingX(A+1) = 1;
            
            remaining = 1; 
            for b = 1:(B-1)
                spacingY(b+1) = spacingYfrac(b)*( (remaining - G.spacingMin*(B-b)) - (G.spacingMin) ) + G.spacingMin+(1-remaining);   %frac*( maxT-minT)+minT
                remaining = 1-spacingY(b+1);
            end
            spacingY(B+1) = 1;
            
            %Converts chrom into map of blocks to raise
            map = reshape(cellArrChrom,A,B);
            
            
            %Make grating
            grating = G.constantGrating;
            
            %grating.stratum{1}.thick = convertChrom_gc(thicknessChrom,G.chromNthickness)*G.thicknessMax;  %Sets thickness of grating
            grating.stratum{1}.thick = 1; %thickness is set by the layer object
            
            for j = 1:G.NblockY %stripes run along x-direction
                grating.stratum{1}.stripe{j}.type = 1;  %inhomogeneous stripe
                grating.stratum{1}.stripe{j}.c1 = spacingY(j+1);  %width of stripe
                blockNum = 1; %the order of the current block
                for i = 1:G.NblockX

                    %A block is filled at the nth level and below if there
                    %are n '1's in that column of the chromosome.


                    %Only create a block right before it switches, to enhance the calculation speed
                    if  (i == G.NblockX) | (map(i,j) ~= map(i+1,j))
                        grating.stratum{1}.stripe{j}.block{blockNum}.c2 = spacingX(i+1); %ending position of block
                        if map(i,j) == 1
                            grating.stratum{1}.stripe{j}.block{blockNum}.pmt_index = G.pmtIndices.filled;
                        else
                            grating.stratum{1}.stripe{j}.block{blockNum}.pmt_index = G.pmtIndices.void;
                        end
                        blockNum = blockNum+1;
                    end
                end
            end
            
            %grating.SU8thickness = convertChrom_gc( SU8thicknessChrom, G.chromNSU8thickness)*(G.SU8thicknessMax-G.SU8thicknessMin) + G.SU8thicknessMin;
            %grating.ITOthickness = convertChrom_gc( ITOthicknessChrom, G.chromNITOthickness)*(G.ITOthicknessMax-G.ITOthicknessMin) + G.ITOthicknessMin;
            
        end
        
        
        
        function chromosomeSize = getChromosomeSize(G)
            %chromosomeSize = G.chromNthickness + G.chromNspacingX*(G.NblockX-1) + G.chromNspacingY*(G.NblockY-1) + G.chromNcellArr + G.chromNSU8thickness + G.chromNITOthickness;
            %chromosomeSize = G.chromNthickness + G.chromNspacingX*(G.NblockX-1) + G.chromNspacingY*(G.NblockY-1) + G.chromNcellArr;
            chromosomeSize = G.chromNspacingX*(G.NblockX-1) + G.chromNspacingY*(G.NblockY-1) + G.chromNcellArr;
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
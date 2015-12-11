classdef GratingGridRhombus
    %Defines a grating according to a mesh, wherein the block in the mesh
    %can freely change in width and height.
    
    
    %The lattice is rhombohedral, but we still use a rectangular unit cell.
    %When GDC plots this grating, each stripe of the unit cell will be
    %offset, but this is not true for S4.
    
    %ASSUMES BINARY GRATING
    
    properties
            rhombohedralAngle;
            period;
            periodX; %Periodicity of grating %um
            periodY; %um  %The y-component of second periodicty vector
            NblockX; % Num blocks in x-dir
            NblockY; % Num blocks in y-dir
            constantGrating; %Beginning of initialization of grating (common to all gratings)
            pmtIndices; %index for permittivities (1=air, 2=grating/SU8)
            chromNcellArr; %Num chromosome positions for on-off of grid
            chromNspacingX; %Num chromosome positions for x-lengths of the blocks
            chromNspacingY; %Num chromosome positions for y-lengths of the blocks
            spacingMin; %Minimum width of each block %um
    end
   
    methods
        function G = GratingGridRhombus(options)
            G.rhombohedralAngle = optinos.rhombohedralAngle;
            
            G.NblockX = options.NblockX;
            G.NblockY = options.NblockY;
            
            G.chromNspacingX = options.chromNspacingX;
            G.chromNspacingY = options.chromNspacingY;
            G.chromNcellArr = G.NblockX * G.NblockY;
            
            %G.periodX = options.periodX;
            %G.periodY = options.periodY;
            G.period = options.period;
            G.periodX = G.period;
            G.periodY = G.period*sin(G.rhombohedralAngle);
            
            G.spacingMin = options.spacingMin;
            
            
            G.constantGrating.pmt={options.n_void^2, options.n_filled^2};  %Permitivity index
            G.pmtIndices.void = 1;
            G.pmtIndices.filled = 2;
            
            G.constantGrating.pmt_sub_index=G.pmtIndices.filled;  %Index of the substrate
            G.constantGrating.pmt_sup_index=G.pmtIndices.void;  %Index of the superstrate
            G.constantGrating.d21 = G.periodX;  %First periodicity vector in x-direction
            G.constantGrating.d31 = 0;  %First periodicity vector in y-direction
                
            G.constantGrating.d22 = 0;  %Second periodicity vector x-direction
            G.constantGrating.d32 = G.periodY;  %Second periodicity vector y-direction
            
                
            G.constantGrating.stratum={}; %container that will hold all the strata
            G.constantGrating.stratum{1}.type=2;%biperiodic
            G.constantGrating.stratum{1}.h11=0;%stratum periodicity is same as grating periodicity
            G.constantGrating.stratum{1}.h12=1;
            G.constantGrating.stratum{1}.h21=1;
            G.constantGrating.stratum{1}.h22=0; 
            
                        
            
            
        end
        
        function grating = generateGrating(G, chromosomeSection)
            
            [spacingXchrom, spacingYchrom, cellArrChrom ] = splitChromosome(chromosomeSection,  ...
                [G.chromNspacingX*(G.NblockX-1), G.chromNspacingY*(G.NblockX-1), G.chromNcellArr]);

            %Generate spacings of the grid
            A = G.NblockX;
            B = G.NblockY;
            spacingX = zeros(1,A+1);  % the block width as a fraction of the periodicity
            spacingY = zeros(1,A+1);
            

            
            spacingXfrac = convertChrom_gc( spacingXchrom, ones(1,A-1)*G.chromNspacingX) ;  %Converts the chromosome data into an array of fractions 
            spacingYfrac = convertChrom_gc( spacingYchrom, ones(1,B-1)*G.chromNspacingY) ;
            
            %spacingX describes the positions of the edge between blocks
            remaining = 1; %describes thickness between most recent spacingX and the end of line
            for a = 1:(A-1)
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
            
            grating.stratum{1}.thick = 1; %thickness is set by the layer object
            
            for j = 1:G.NblockY %stripes run along x-direction
                grating.stratum{1}.stripe{j}.type = 1;  %inhomogeneous stripe
                grating.stratum{1}.stripe{j}.c1 = spacingY(j+1);  %width of stripe
                blockNum = 1; %the order of the current block
                for i = 1:G.NblockX
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
            
            
        end
        
        
        
        function chromosomeSize = getChromosomeSize(G)
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
            
            x_limit=[-0.5,0,0;1.0,G.periodX*2,G.periodY*2];
            h_plot=gdc_plot(grating,1,pmt_display,x_limit);
            view(0,90);
        end
       
    end
    
end
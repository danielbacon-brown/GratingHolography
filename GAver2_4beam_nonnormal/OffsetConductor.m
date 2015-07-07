classdef OffsetConductor
    %Discretely translates an 3D array of values, assuming the array is periodic at
    %every edge.
    %'hexagonal' assumes that x has ordinary periodicity
    
    properties
        offsetType;  %Type of periodicity ('hexagonal', 'square', 'parallelogram')
        chromNoffsetX;  %# of positions used in the chromosome
        chromNoffsetY;
        chromNoffsetZ;
        cells;  %Number of discrete units of matrix
    end
    
    methods
        
        function O = OffsetConductor(options)
            O.offsetType = options.offsetType;
            O.chromNoffsetX = options.chromNoffsetX;
            O.chromNoffsetY = options.chromNoffsetY;
            O.chromNoffsetZ = options.chromNoffsetZ;
            O.cells = options.cells;
        end
        
        
        function matrixOut = doOffset(O,matrixIn, offsetChromosome)
            
            [offset_x, offset_y, offset_z] = convertChrom_gc(offsetChromosome,[O.chromNoffsetX,O.chromNoffsetY,O.chromNoffsetZ]); %Returns offsets as fractions of unit cell
            %shift remaining cells to fill in the missing space)
            offset_x = floor( offset_x*O.cells(1) ); %discrete number of cells to offset
            rem_x = O.cells(1) - offset_x;  %remaining cells that need to be shifted down
            offset_y = floor( offset_y*O.cells(2) );
            rem_y = O.cells(2) - offset_y;
            offset_z = floor( offset_z*O.cells(3) )
            rem_z = O.cells(3) - offset_z
            
            
            switch O.offsetType
                case 'hexagonal'
                    
                    
                    %offsetX
                    tempVol = matrixIn( 1:offset_x, :,:);  %left section to be moved to right end
                    matrixIn(1:rem_x, :,:) = matrixIn((offset_x+1):end,:,:); %shift right section over
                    matrixIn( rem_x+1:end, :,:) = tempVol;  %fill in right end
                    
                    
                    %offsetY:
                    tempVol1 = matrixIn( 1:floor(end/2), 1:offset_y, :); %bottom left corner
                    tempVol2 = matrixIn( floor(end/2)+1:end, 1:offset_y, :); % bottom right corner
                    
                    matrixIn(:, 1:rem_y, :) = matrixIn(:, offset_y+1:end, :); %shift bottom section over
                    
                    matrixIn( 1:floor(end/2), rem_y+1:end, :) = tempVol2; % top left corner
                    matrixIn( floor(end/2)+1:end, rem_y+1:end, :) = tempVol1; %top right corner
                    
                    
                    %offsetZ:
                    tempVol = matrixIn(:,:, 1:offset_z);  %lower section to be moved to upper end
                    matrixIn(:,:, 1:rem_z) = matrixIn(:,:, (offset_z+1):end); %shift upper section down
                    matrixIn(:,:, rem_z+1:end) = tempVol; %fill in upper section
                    
                case {'square','parallogram'}
                    
                    %offsetX
                    tempVol = matrixIn( 1:offset_x, :,:);  %left section to be moved to right end
                    matrixIn(1:rem_x, :,:) = matrixIn((offset_x+1):end,:,:); %shift right section over
                    matrixIn( rem_x+1:end, :,:) = tempVol;  %fill in right end
                    
                    
                    %offsetY:
                    tempVol = matrixIn( :, 1:offset_y, :);  %left section to be moved to right end
                    matrixIn(:,1:rem_y, :) = matrixIn(:,(offset_y+1):end,:); %shift right section over
                    matrixIn(:, rem_y+1:end,:) = tempVol;  %fill in right end
                    
                    
                    %offsetZ:
                    sizeall = size(matrixIn)
                    tempVol = matrixIn(:,:, 1:offset_z);  %lower section to be moved to upper end
                    sizetemp = size(tempVol)
                    sizeL = size( matrixIn(:,:, 1:rem_z) )
                    sizeR = size( matrixIn(:,:, (offset_z+1):end) )
                    %testout = matrixIn(:,:, (offset_z+1):end)
                    matrixIn(:,:, 1:rem_z) = matrixIn(:,:, (offset_z+1):end); %shift upper section down
                    matrixIn(:,:, rem_z+1:end) = tempVol; %fill in upper section
                    
                    
                    
                    
            end
            
            matrixOut = matrixIn;
            
        end
        
        
        function chromosomeSize = getChromosomeSize(O)
            chromosomeSize = O.chromNoffsetX + O.chromNoffsetY + O.chromNoffsetZ;
        end
    end
    
end
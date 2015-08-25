%Outputs the structure in n,k format.

function [] = outputLumericalFile( structure,mapdim, n_exposed,n_exposed_cplx,n_inverse,n_inverse_cplx, inverse,filename )
%structure = 3d matrix of the SU8 structure (1 is SU8)
%threshold = the minimum intensity that will cause crosslinking
%mapdim = dimensions of map [xdim,ydim,zdim] (preferably microns)
%N = real refractive index of structure
%k = imaginary refractive index of structure
%inverse = boolean;      inverse=1 --> inverted structure
%filename = file to output to


mapsize = size(structure);  %number of cells in the map in each dimension

ef=fopen(filename,'w');  %open/create file
fprintf(ef,'%G %G %G \n',mapsize(1),-mapdim(1)/2,mapdim(1)/2); % #blocks in x direction , min x coordinate, max x coordinate
fprintf(ef,'%G %G %G \n',mapsize(2),-mapdim(2)/2,mapdim(2)/2); % #blocks in y direction , min y coordinate, max y coordinate
fprintf(ef,'%G %G %G \n',mapsize(3),-mapdim(3)/2,mapdim(3)/2); % #blocks in z direction , min z coordinate, max z coordinate

if inverse==1  %If 1 corresponds to void/inverted and 0 corresponds to SU8
    for k=1:mapsize(3)
        for j=1:mapsize(2)
            for i=1:mapsize(1)
                %if structure(i,j,k)<threshold
                if ~structure(i,j,k) %if void
                    fprintf(ef,'%G %G \n', n_exposed, n_exposed_cplx);
                else %if SU8
                    fprintf(ef,'%G %G \n', n_inverse, n_inverse_cplx);
                end
            end
        end
    end
    
else %If 0 corresponds to void/inverted and 1 corresponds to SU8
    
    for k=1:mapsize(3)
        for j=1:mapsize(2)
            for i=1:mapsize(1)
                %if structure(i,j,k)>threshold
                if structure(i,j,k)
                    fprintf(ef,'%G %G \n', n_exposed, n_exposed_cplx);
                else
                    fprintf(ef,'%G %G \n', n_inverse, n_inverse_cplx);
                end
            end
        end
    end
end


end


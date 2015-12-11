function structure = generateTwistedCrossStructure(options)
%Generates a twisted cross structure, as described in M. Decker, M. Ruther, C.E. Kriegler, J. Zhou, C.M. Soukoulis, S. Linden, and M. Wegener, Optics Letters 34, 2501 (2009).

%fields: cells [Nx,Ny,Nz]
%All units are um
%   radb = big radius / distance of center of wire to axis 
%   radl = little radius / distance of center of wire to surface of wire 
%   u = first periodicity vector [x,y,z] %um
%   v = second periodicity vector [x,y,z] %um
%   w = third periodicity vector [x,y,z] %um
% 1 = in helix

struct = zeros(options.cells); 

%Grid describing location of each point on
linU = linspace(0,1,options.cells(1));  
linV = linspace(0,1,options.cells(2));
linW = linspace(0,1,options.cells(3));
[gridU, gridV, gridW] = ndgrid(linU,linV,linW); %um

gridX = gridU*options.u(1) + gridV*options.v(1) + gridW*options.w(1);
gridY = gridU*options.u(2) + gridV*options.v(2) + gridW*options.w(2);
gridZ = gridU*options.u(3) + gridV*options.v(3) + gridW*options.w(3);

relativeZ = options.relativeZ; %Describes z-dimension of sphere as fraction of diameters

for i_offset_u = -1:1
    for i_offset_v = -1:1
        for i_offset_w = -1:1

            offset_x = i_offset_u*options.u(1) + i_offset_v*options.v(1) + i_offset_w*options.w(1);
            offset_y = i_offset_u*options.u(2) + i_offset_v*options.v(2) + i_offset_w*options.w(2);
            offset_z = i_offset_u*options.u(3) + i_offset_v*options.v(3) + i_offset_w*options.w(3);
            
            %Central pillar
            for theta = linspace(0,2*pi,100)  %iterate by z (copied from helical)
                %Center of circle
                circX =  (options.u(1)+options.v(1)+options.w(1))/2 + offset_x;  %um
                circY =  (options.u(2)+options.v(2)+options.w(2))/2 + offset_y; %um
                circZ = theta/(2*pi)*options.w(3) + offset_z; %um
                %Adds circles to struct
                struct = struct | ( (circX-gridX).^2 + (circY-gridY).^2 + ((circZ-gridZ)/relativeZ).^2  <= options.radl^2 );  %Sets all points within circle to 1
            end
            
            
            %First line - at zero
            
            
            
            
        end
    end
end




end
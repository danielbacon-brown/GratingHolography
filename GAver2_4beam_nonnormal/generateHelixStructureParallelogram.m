function struct = generateHelixStructureParallelogram(options)
%Generates a 3-D helix matrix according to a parallelogramic grid
%fields: cells [Nx,Ny,Nz]
%All units are um
%   radb = big radius / distance of center of wire to axis 
%   radl = little radius / distance of center of wire to surface of wire 
%   u = first periodicity vector [x,y,z] %um
%   v = second periodicity vector [x,y,z] %um
%   w = third periodicity vector [x,y,z] %um

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

angleOffset = 135;

%Makes 9 helices on each side (to make sure the targetstructure is periodic
for i_offset_u = -1:1
    for i_offset_v = -1:1
        for i_offset_w = -1:1

            offset_x = i_offset_u*options.u(1) + i_offset_v*options.v(1) + i_offset_w*options.w(1);
            offset_y = i_offset_u*options.u(2) + i_offset_v*options.v(2) + i_offset_w*options.w(2);
            offset_z = i_offset_u*options.u(3) + i_offset_v*options.v(3) + i_offset_w*options.w(3);
            
            for theta = linspace(0,2*pi,100)  %iterate by angle
                %Center of circle
                circX = options.radb * cos(theta-angleOffset) + (options.u(1)+options.v(1)+options.w(1))/2 + offset_x;  %um
                circY = options.radb * sin(theta-angleOffset) + (options.u(2)+options.v(2)+options.w(2))/2 + offset_y; %um
                circZ = theta/(2*pi)*options.w(3) + offset_z; %um
                %Adds circles to struct
                struct = struct | ( (circX-gridX).^2 + (circY-gridY).^2 + ((circZ-gridZ)/relativeZ).^2  <= options.radl^2 );  %Sets all points within circle to 1
            end

        end
    end
end


% %TEST
figure
patched = patch(isosurface(padarray(struct,[1,1,1]),0.5));
set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
view(3);
camlight
axis equal
lighting gouraud
xlim([1 options.cells(1)+2])
ylim([1 options.cells(2)+2])
zlim([1 options.cells(3)+2])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[])



end
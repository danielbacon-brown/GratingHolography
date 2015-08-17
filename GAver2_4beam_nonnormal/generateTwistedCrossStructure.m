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






end
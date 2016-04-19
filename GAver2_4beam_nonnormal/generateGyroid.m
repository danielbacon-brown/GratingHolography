function gyroidStruct = generateGyroid(cells,fill)

%Generate gyroid structure:
%options.cells = [25,25,25];
options.cells = cells;

%struct = zeros(options.cells); 

%Grid describing location of each point on
linU = linspace(0,2*pi,options.cells(1));  
linV = linspace(0,2*pi,options.cells(2));
linW = linspace(0,2*pi,options.cells(3));
[gridU, gridV, gridW] = ndgrid(linU,linV,linW); %um

%[meshx,meshy,meshz] = ndmeshgrid(x,y,z)

gyroidI = sin(gridU).*cos(gridV) + sin(gridV).*cos(gridW) + sin(gridW).*cos(gridU);

%targetfill = 0.50;
targetfill = fill;
threshold = fixfill(gyroidI, 256, targetfill);
gyroidStruct = gyroidI < threshold;

%TEST
figure
patched = patch(isosurface(padarray(gyroidStruct,[1,1,1]),0.5));
set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');
view(3);
camlight
axis equal
lighting gouraud
xlim([1 options.cells(1)+2])
ylim([1 options.cells(2)+2])
zlim([1 options.cells(3)+2])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[]);


end


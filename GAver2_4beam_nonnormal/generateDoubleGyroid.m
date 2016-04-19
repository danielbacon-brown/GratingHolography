function gyroidStruct = generateDoubleGyroid(options)

%Generate gyroid structure:
%options.cells = [25,25,25];
options.cells = [50,50,50];


%struct = zeros(options.cells); 

%Grid describing location of each point on
linU = linspace(0,3*pi,options.cells(1));  
linV = linspace(0,3*pi,options.cells(2));
linW = linspace(0,3*pi,options.cells(3));
[gridU, gridV, gridW] = ndgrid(linU,linV,linW); %um

%Reshape grid so that it can be rotated:
gridUr = reshape(gridU,1,[]);
gridVr = reshape(gridV,1,[]);
gridWr = reshape(gridW,1,[]);
gridUVWr = [gridUr;gridVr;gridWr];

%Rotation transformation for grid so that [1 1 1] direction is pointing vertical
theta1 = 45*pi/180;
R1 = [  1, 0, 0; ...
        0, cos(theta1), -sin(theta1); ...
        0, sin(theta1), cos(theta1) ];
theta2 = -1*atan(1/sqrt(2));
R2 = [  cos(theta2), 0, sin(theta2); ...
        0, 1, 0; ...
        -sin(theta2), 0, cos(theta2) ];

    %Display new periodicity parameters:
    Lx = [ 1;0;0];
    Ly = [ 0;1;0];
    Lz = [ 0;0;1];
    LxPrime = R2*R1*Lx
    LyPrime = R2*R1*Ly
    LzPrime = R2*R1*Lz
    figure
    hold on
    %origin
    quiver3(0,0,0,LxPrime(1),LxPrime(2),LxPrime(3))
    quiver3(0,0,0,LyPrime(1),LyPrime(2),LyPrime(3))
    quiver3(0,0,0,LzPrime(1),LzPrime(2),LzPrime(3))

    originUnits = [0,0,0;  1,0,0;  0,1,0;  0,0,1;  1,1,0; 1,0,1; 0,1,1; 1,1,1]';
    originVectors = R2*R1*originUnits;
            
    for ii = 1:size(originUnits,2)
        %quiver3(originPoints(ii,1)*LxPrime(1), originPoints(ii,2)*LxPrime, originPoints(ii,3), 
        quiver3(originVectors(1,ii),originVectors(2,ii),originVectors(3,ii),LxPrime(1),LxPrime(2),LxPrime(3))
        quiver3(originVectors(1,ii),originVectors(2,ii),originVectors(3,ii),LyPrime(1),LyPrime(2),LyPrime(3))
        quiver3(originVectors(1,ii),originVectors(2,ii),originVectors(3,ii),LzPrime(1),LzPrime(2),LzPrime(3))
        %quiver3(originVector(1),originVector(2),originVector(3),
    end
    
    %Generating U,V,W quivers
    originUnitsUVW = [0,0,0;  1,0,0;  0,1,0;  0,0,1;  1,1,0; 1,0,1; 0,1,1; 1,1,1]';
    Lu = [ 1.2247, 0.7071, 0];
    Lv = [ 0, 1.4142, 0];
    Lw = [ 0,0, 1.7332];
    Luvw = [Lu;Lv;Lw]';
    originVectorsUVW = Luvw * originUnitsUVW;
    for ii = 1:size(originUnits,2)
    	quiver3(originVectorsUVW(1,ii)+LzPrime(1),originVectorsUVW(2,ii)+LzPrime(2),originVectorsUVW(3,ii)+LzPrime(3),Lu(1),Lu(2),Lu(3))
        quiver3(originVectorsUVW(1,ii)+LzPrime(1),originVectorsUVW(2,ii)+LzPrime(2),originVectorsUVW(3,ii)+LzPrime(3),Lv(1),Lv(2),Lv(3))
        quiver3(originVectorsUVW(1,ii)+LzPrime(1),originVectorsUVW(2,ii)+LzPrime(2),originVectorsUVW(3,ii)+LzPrime(3),Lw(1),Lw(2),Lw(3))
    end
    xlim([-1 1.5])
    ylim([-1 1.5])
    zlim([0 3])

%MAKE UNIT CELL PATTERN
    %Create grid using UVW vectors (corresponds to unit cell
    Uunits = linspace(0,1.4142*2*pi,options.cells(1)+1);
    Uunits = Uunits(1:(end-1));  %otherwise last point is doublecounted (due to periodicity)
    Vunits = linspace(0,1.4142*2*pi,options.cells(2)+1);
    Vunits = Vunits(1:(end-1));
    Wunits = linspace(0,1.4142*2*pi,options.cells(3)+1);
    Wunits = Wunits(1:(end-1));
    [UunitGrid, VunitGrid, WunitGrid] = ndgrid(Uunits,Vunits,Wunits);
    
    %Transform grid into cartesian coordinates (where 1,1,1 is vertical)
    UVWunitGrid = [ reshape(UunitGrid,1,[]);  reshape(VunitGrid,1,[]); reshape(WunitGrid,1,[]) ];
    UVWgrid = Luvw * UVWunitGrid;
    %Calculate 'potential' of pattern
    UVWpotential = sin(UVWgrid(1,:)).*cos(UVWgrid(2,:)) + sin(UVWgrid(2,:)).*cos(UVWgrid(3,:)) + sin(UVWgrid(3,:)).*cos(UVWgrid(1,:));
    size(UVWpotential)
    %Reshape grid
    potentialGridR = reshape( UVWpotential(1,:), options.cells(1),options.cells(2),options.cells(3)); 
    
%TEST
figure
patched = patch(isosurface(padarray(potentialGridR,[1,1,1]),0.5));
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
    
    
    
    
%Apply rotations:
gridUVWrPrime = R1*R2*gridUVWr;

%Separate and reshape coordinate vectors:
gridUprime = reshape(gridUVWrPrime(1,:), options.cells(1),options.cells(2),options.cells(3));
gridVprime = reshape(gridUVWrPrime(2,:), options.cells(1),options.cells(2),options.cells(3));    
gridWprime = reshape(gridUVWrPrime(3,:), options.cells(1),options.cells(2),options.cells(3));   


%Calculate 'potential' for gyroid function
gyroidI = sin(gridUprime).*cos(gridVprime) + sin(gridVprime).*cos(gridWprime) + sin(gridWprime).*cos(gridUprime);
%gyroidI = sin(gridU).*cos(gridV) + sin(gridV).*cos(gridW) + sin(gridW).*cos(gridU);
% gyroidI2 = rot90(gyroidI,2);
% gyroidI3 = gyroidI + gyroidI2;
%gyroidI2 = sin(gridU2).*cos(gridV2) + sin(gridV2).*cos(gridW2) + sin(gridW2).*cos(gridU2);

targetfill = 0.5;
threshold = fixfill(gyroidI, 256, targetfill);
gyroidStruct1 = (gyroidI) > threshold;
%gyroidStruct2 = rot90(gyroidStruct1,2);
threshold2 = fixfill(gyroidI, 256, 1-targetfill);
gyroidStruct2 =  gyroidI < threshold2;
gyroidStruct3 = gyroidStruct1 + gyroidStruct2;

%TEST
figure
patched = patch(isosurface(padarray(gyroidStruct3,[1,1,1]),0.5));
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


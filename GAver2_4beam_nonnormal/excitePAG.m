function acidDist = excitePAG(intensityDist,dimensions,sensDens,absCrossSection)
%Returns the number of excited PAG in each cell
%intensityDist: based off of power of 1/area using S4's units
%sensDens: scalar density of the PAG in mol/um^2
%dimensions: [x,y,z] in microns
%absCrossSection: scalar in 1/um^2

N_a = 6.022e23; %Avogadro's Number


%%% First figure out a random distribution of PAG:
vol_unit = dimensions(1)*dimensions(2)*dimensions(3); %um^3
moleculesPerUnit = N_a*sensDens*vol_unit %molecules/unit cell

%lambda_poisson = average # molecules per cell < 10 %in order for Poisson distribution to be approximately right 
%So a^3 > #molecules/unit / #cells/unit
a = ceil( (moleculesPerUnit/(cells(1)*cells(2)*cells(3)))^(1/3) );

cellsA = cells*a; %increase number of cells to better approximate Poisson's distribution

%populate sensitizer molecules:
sensCount = poissrnd(poisson_lambda, cellsA,cellsA,cellsA);  %Create higher density cell array for more accurate distribution
sensSum = sum(sum(sum(sensCount))) %Should match moleculesPerUnit

if a > 1
    %Reduce number of cells in sensCount to match 
    sensCount_reduced = zeros(Nx,Ny,Nz);
    for i_i = 1:Nx
        for j_i = 1:Ny
            for k_i = 1:Nz
                sensCount_reduced(i_i,j_i,k_i) = sum(sum(sum( sensCount(a*i_i-2:a*i_i, a*j_i-2:3*j_i, a*k_i-2:a*k_i) ))); %gets total number of molecules in larger cell
            end
        end
    end
    sensCount = sensCount_reduced;
end







end
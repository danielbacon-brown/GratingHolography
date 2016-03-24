%function acidCount = excitePAG(intensityDist,dimensions,sensDens,absCrossSection, t_exposure)
function acidCount = excitePAG(intensityDist,dimensions,sensDens,QYtimesAbsCrossSection, t_exposure)

plotHeatMaps = 0;

%Returns the number of excited PAG in each cell
%intensityDist: W/m^2
%sensDens: scalar density of the PAG in molecules/um^3
%dimensions: [x,y,z] in microns
%%%%absCrossSection: scalar in 1/nm^2
%beamPowerDens: power of the laser (W/m^2)
%QYtimesAbsorptionCrossSection: m^2


INmax = max(max(max(intensityDist)))
INmin = min(min(min(intensityDist)))


N_a = 6.022e23; %Avogadro's Number
c = 2.99e8; %speed of light %m/s
eps_0 = 8.854e-12; %vacuum permittivity %F/m  %length units of c and eps_0 will cancel
n = 1.58; %refractive index of SU8
q = 1.6e-19; %Coulombs % fundamental charge
h = 6.626e-34; %Planck's constant %J*s

%quantum_efficiency = 0.07; %Fraction of absorbed photons leading to 
%t_exposure = 100;
E_photon = h*c/(0.532e-6); %J   (=2.33eV) %532nm

%%% First figure out a random distribution of PAG:
vol_unit = dimensions(1)*dimensions(2)*dimensions(3); %um^3
moleculesPerUnit = sensDens*vol_unit %molecules/unit cell   %sensdens in molecules/um^3

cells = size(intensityDist);

%lambda_poisson = average # molecules per cell < 10 %in order for Poisson distribution to be approximately right 
%So a^3 > #molecules/unit / #cells/unit
a = ceil( (moleculesPerUnit/(cells(1)*cells(2)*cells(3)))^(1/3) )

cellsA = cells*a; %increase number of cells to better approximate Poisson's distribution
poisson_lambda = moleculesPerUnit/(cellsA(1)*cellsA(2)*cellsA(3)) %lambda used in poisson distribution

%populate sensitizer molecules:
sensCount = poissrnd(poisson_lambda, cellsA(1),cellsA(2),cellsA(3));  %Create higher density cell array for more accurate distribution
sensSum = sum(sum(sum(sensCount))) %Should match moleculesPerUnit

if a > 1
    %Reduce number of cells in sensCount to match 
    sensCount_reduced = zeros(cells(1),cells(2),cells(3));
    for i_i = 1:cells(1)
        for j_i = 1:cells(2)
            for k_i = 1:cells(3)
                sensCount_reduced(i_i,j_i,k_i) = sum(sum(sum( sensCount((a*(i_i-1)+1):(a*i_i), (a*(j_i-1)+1):(a*j_i), (a*(k_i-1)+1):(a*k_i)) ))); %gets total number of molecules in larger cell
            end
        end
    end
    sensCount = sensCount_reduced;
end
maxSenscount = max(max(max(sensCount)))
minSensCount = min(min(min(sensCount)))
sensSum = sum(sum(sum(sensCount))) %Should match moleculesPerUnit

%PLOT HEATMAP:
if plotHeatMaps == 1
    figure
    colormap(hot)
    imagesc(sensCount(:,:,20), [0,max(max(sensCount(:,:,20)))])
    colorbar
end


%Calculate the real electric field
%  S4 gives electric field in units of V/sqrt(A), assuming incident beam power of 1
%  W/A, where A is the area of the unit cell (S4 later divides this power by
%  cos(phi) to account for cosine law)
% E_r = E_s*sqrt(I_r/I_s)
%I_dist = c*n*eps_0/2*(I_r/I_s)*abs(E_s)^2

%I_r_over_I_s = beamPowerDens * (dimensions(1)*1e-6*dimensions(2)*1e-6) %Scale by real beam power %W/A %(I_s is 1W/unitcellarea)
%I_r_over_I_s = beamPowerDens / (1 / (dimensions(1)*1e-6*dimensions(2)*1e-6) )
%intensityDist_r = c*n*eps_0/2*I_r_over_I_s*intensityDist;  %Calc real interference intensity %W/um^2

%%% E is in units of V/sqrt(A), so intensityDist is given in W/A
%intensityDist_r = I_r_over_I_s*intensityDist;  %Calc real interference intensity %W/m^2
%intensityDist_r = intensityDist/(dimensions(1)*dimensions(2))  * 1e12;  %W/A * (1A/um^2) * (um^2/m^2) => W/m^2
%Actually avoiding normalization by unit cell area
intensityDist_r = intensityDist;
INrmax = max(max(max(intensityDist_r)))
INrmin = min(min(min(intensityDist_r)))

%flux_photon = intensityDist_r*t_exposure/E_photon/1e18;   %Distribution of photons at each points %#photons/nm^2        %W/m^2 * s / (J/photon) * (nm/m)^2
%pfmax = max(max(max(flux_photon)))
%pfmin = min(min(min(flux_photon)))

%prob_abs = 1 - exp(-absCrossSection.*quantum_efficiency.*flux_photon); %Probability at least one photon is absorbed at each point   %nm^2/molecule * 1 * #photons/nm^2
%prob_abs = 1 - exp(-QYtimesAbsCrossSection.*flux_photon);    % ~8e-7nm^2  * #photons/nm^2 -> unitless
prob_abs = 1 - exp(-QYtimesAbsCrossSection.*intensityDist_r.*t_exposure./E_photon);    % m^2/molecule  * #W/m^2 * s / J ->  1/molecules
PAmax = max(max(max(prob_abs)))
PAmin = min(min(min(prob_abs)))


acidCount = zeros(size(intensityDist_r));
for i_x = 1:cells(1)
    for i_y = 1:cells(2)
        for i_z = 1:cells(3)
            acidCount(i_x,i_y,i_z) = binornd( sensCount(i_x,i_y,i_z), prob_abs(i_x,i_y,i_z) ); %Binomial distribution to give # of excited acid in each volume
        end
    end
    disp(i_x)
end

percent_excited = sum(sum(sum(acidCount))) / sum(sum(sum(sensCount)))


    %Plot histogram of intensities
    figure
    hist(reshape(intensityDist_r,1,[]),100)
    
    
    %PLOT HEATMAP:
    if plotHeatMaps == 1
        figure
        colormap(hot)
        imagesc(acidCount(:,:,20), [0,max(max(acidCount(:,:,20)))])
        colorbar
    end
    
    
    


end
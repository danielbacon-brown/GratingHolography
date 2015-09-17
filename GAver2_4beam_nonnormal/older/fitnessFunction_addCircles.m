function fitness = fitnessFunction_addCircles( chromosome, GAoptions )
%This is the fitness function for the GA

    tfitness = tic;
    
    %Read chromosomal data
    chromosomalData = readChromosomalData(chromosome,GAoptions);
    
    
    
    %MAKE GRATING
    grating = makeGrating_addCircles(chromosomalData, GAoptions);
    
    intensityDist = calcIntensities(grating,chromosomalData,GAoptions);
    
    
    %Normalizes the threshold
    maxIntensity = max(real(reshape(intensityDist,GAoptions.Ncell,1)));
    minIntensity = min(real(reshape(intensityDist,GAoptions.Ncell,1)));
    threshold = chromosomalData.thresholdfraction*(maxIntensity-minIntensity)+minIntensity;
    
    %Convert the intensity distribution into filling data of SU8
    structureSU8 = intensityDist > threshold;  
    
    %Determine whether the structure is bicontinuous (returns 1) and if not, returns approx. how bicontinuous it is    
    [bicontinuousFactor, finalF, finalE] = ContinuityTest(structureSU8, GAoptions.useRealisticInversion);  
      
    %Calculate the bicontinuous factors of structures with different thresholds
    thresholdLow = (2*threshold+minIntensity)/3  %between actual and minimum-possible thresholds, but weighted towards actual
    thresholdHigh = (2*threshold+maxIntensity)/3  
    bicontinuousFactorLowThresh = ContinuityTest(intensityDist > thresholdLow, GAoptions.useRealisticInversion);
    bicontinuousFactorHighThresh = ContinuityTest(intensityDist > thresholdHigh, GAoptions.useRealisticInversion);
        
    
    %Save the structure to a file
    if GAoptions.hexagonalGrating
        if GAoptions.useRealisticInversion
            outputLumericalFile(finalE,[GAoptions.period,2*GAoptions.yperiod,GAoptions.Z_T_shrunk],GAoptions.n_Ag, GAoptions.k_Ag, 0, [GAoptions.dir,GAoptions.currentNKfile]); %0=don't invert
        else
            outputLumericalFile(structureSU8,[GAoptions.period,2*GAoptions.yperiod,GAoptions.Z_T_shrunk], GAoptions.n_Ag, GAoptions.k_Ag, 1, [GAoptions.dir,GAoptions.currentNKfile]); %1=invert
        end
    else %square grating
        if GAoptions.useRealisticInversion
            outputLumericalFile(finalE,[GAoptions.period,GAoptions.period,GAoptions.Z_T_shrunk],GAoptions.n_Au, GAoptions.k_Au, 0, [GAoptions.dir,GAoptions.currentNKfile]); %0=don't invert
        else
            outputLumericalFile(structureSU8,[GAoptions.period,GAoptions.period,GAoptions.Z_T_shrunk], GAoptions.n_Au, GAoptions.k_Au, 1, [GAoptions.dir,GAoptions.currentNKfile]); %1=invert
        end
    end
    
    
    
    %If the 'pause' file is found, shift control to keyboard, pausing the
    %simulation.  Useful if someone else needs to use Lumerical.
    if exist('/home/danielbacon-brown/pause','dir') == 7
        keyboard;
    end

    
    
    %Remove old Lumerical output file (if it exists)
    if exist([GAoptions.dir,GAoptions.currentLumResultsFile], 'file')
        delete([GAoptions.dir, GAoptions.currentLumResultsFile]);
    end
    
    %Tell Lumerical to run the script
    system(['fdtd-solutions -nw -run ',GAoptions.dir,GAoptions.LumRunScript]);
    %Wait for the appearance of the Lumerical output file
    while (~exist([GAoptions.dir,GAoptions.currentLumResultsFile], 'file'))
        pause(0.1);
    end
    
    %Analyze the output data
    LumResults = load([GAoptions.dir,GAoptions.currentLumResultsFile]);
    transmission_right = LumResults.transmission_right;
    transmission_left = LumResults.transmission_left;
    trans_right = sum(transmission_right(1:end));
    trans_left = sum(transmission_left(1:end));
    trans_delta = abs(trans_right-trans_left)
    
    if bicontinuousFactor == 1
        bicontinuousFitness = 5^2;  %max possible value (because it being centered at 1 is all you need
    else
        bicontinuousFitness = (3*bicontinuousFactor + bicontinuousFactorLowThresh + bicontinuousFactorHighThresh)^2;
    end
    %bicontinuousFitness = (3*bicontinuousFactor + bicontinuousFactorLowThresh + bicontinuousFactorHighThresh)^2
    fitness = -1 * trans_delta * bicontinuousFitness %More negative fitness = better fit

    toc(tfitness)
end


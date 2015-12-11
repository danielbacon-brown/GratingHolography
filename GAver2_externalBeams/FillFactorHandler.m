classdef FillFactorHandler
   %This keeps track of a variable fill factor determined by chromosomal inputs 
   
    properties
        constFill;  % <0 for variable
        fillMin;
        fillMax;
        chromNfill;
        nbins;
    end
    
    
    methods
        function F = FillFactorHandler( options)
            if options.constFill < 0 %Variable fill
                F.constFill = options.constFill;
                F.fillMin = options.fillMin;
                F.fillMax = options.fillMax;
                F.chromNfill = options.chromNfill;
            else %constant fill
                F.constFill = options.constFill;
                F.chromNfill = 0;
            end
            F.nbins = 256;
            
        end
        
        function chromosomeSize = getChromosomeSize(F)
            chromosomeSize = F.chromNfill;
        end
        
        function fill = getFill(F,chromosome)
            if F.constFill >= 0
                fill = F.constFill;
            else
                fill = convertChrom_gc( chromosome, F.chromNfill)*(F.fillMax-F.fillMin) + F.fillMin;
            end
                
        end
        
        
        function [exposedStruct,fill,threshold] = applyFill(F,chromosome,intensityDist) 
            fill = getFill(F,chromosome);
            threshold = fixfill(intensityDist,F.nbins,fill); %Calculates the threshold value that will yield desired fill fraction
            exposedStruct = intensityDist>threshold;  %1 if intensity above threshold (SU8), 0 if below (void)
            
        end
    end
    
end
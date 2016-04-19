classdef Material
    
    properties
        materialName;
        constRI; % <0  =>  constant thickness
        RIMin;
        RIMax;
        chromNRI;
    end
    
    
    methods
        function M = Material( materialName, inconstRI, RIMin, RIMax, chromNRI)
            if real(inconstRI) < 0 %Variable RI
                M.RIMin = RIMin;
                M.RIMax = RIMax;
                M.chromNRI = chromNRI;
                M.materialName = materialName;
            else %constant thickness
                M.materialName = materialName;
                M.constRI = inconstRI;
                M.chromNRI = 0;
            end
            
        end
        
        function chromosomeSize = getChromosomeSize(M)
            chromosomeSize = M.chromNRI;
        end
        
        function [n,k] = getRI(M,chromosome)
            if M.constRI >= 0
                n = real(M.constRI);
                k = imag(M.constRI);
            else
                n = convertChrom_gc( chromosome, M.chromNRI)*(M.RIMax-M.RIMin) + M.RIMin; %um
                k = 0;
            end
                
        end
        
    end
    
end
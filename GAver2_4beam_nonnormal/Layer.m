classdef Layer
    
    properties
        layerName;
        constThick; % <0  =>  constant thickness
        thickMin;
        thickMax;
        chromNthick;
        materialName;
    end
    
    
    methods
        function L = Layer( layerName, materialName, constThick, thickMin, thickMax, chromNthick)
            if constThick < 0 %Variable thickness
                L.layerName = layerName;
                L.thickMin = thickMin;
                L.thickMax = thickMax;
                L.chromNthick = chromNthick;
                L.materialName = materialName;
            else %constant thickness
                L.layerName = layerName;
                L.materialName = materialName;
                L.constThick = constThick;
                L.chromNthick = 0;
            end
            
        end
        
        function chromosomeSize = getChromosomeSize(L)
            chromosomeSize = L.chromNthick;
        end
        
        function thickness = getThickness(L,chromosome)
            if L.constThick >= 0
                thickness = L.constThick;
            else
                thickness = convertChrom_gc( chromosome, L.chromNthick)*(L.thickMax-L.thickMin) + L.thickMin; %um
            end
                
        end
        
    end
    
end
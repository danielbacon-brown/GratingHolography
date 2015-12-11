function writeLumericalRunFileSkin( options, structSkin)
%This writes the Lumerical Script that imports the structure, runs the
%simulation, and exports the data.
%Assume that if the structure is hexagonal, the structure includes 2
%periods in Y.
%Uses a "skin" structure, where 0=SU8, 1=void, 2=skin

runout = fopen([options.dir,options.LumRunScript], 'w');


Nx = size(structSkin,1);
Ny = size(structSkin,2);
Nz = size(structSkin,3);

%Create repeat structure
structSkin = repmat(structSkin,[1,1,options.fdtd.repeatUnits]);




%SETUP
fprintf(runout,'newproject; \n');

structureThickness = options.fdtd.shrinkFactor*options.fdtd.repeatUnits*options.dimensions(3);

%Mesh
fprintf(runout,'addfdtd; \n');
setXY(runout,options);
fprintf(runout,'set("z", 0); \n');
fprintf(runout,'set("simulation time",%i); \n', options.fdtd.simulationTime*1e-15);
fprintf(runout,'set("z span", %e ); \n', structureThickness*3*1e-6);
fprintf(runout,'set("mesh type","auto non-uniform"); \n');
fprintf(runout,'set("mesh accuracy",%i); \n',options.fdtd.meshAccuracy);
fprintf(runout,'set("use early shutoff",0); \n');
fprintf(runout,'set("x min bc","periodic"); \n');
fprintf(runout,'set("y min bc","periodic"); \n');

%Plane Source
fprintf(runout,'addplane; \n');
fprintf(runout,'set("name","source1"); \n');
setXY(runout,options);
fprintf(runout,'set("z", %e); \n', (-structureThickness/2-0.3)*1e-6);
fprintf(runout,'set("wavelength start",%e); \n', options.fdtd.minSourceWL);
fprintf(runout,'set("wavelength stop", %e); \n', options.fdtd.maxSourceWL);

%Plane Source for Circular Polarization  %Right-handed (I think)
fprintf(runout,'copy(0,0,0); \n');
fprintf(runout,'set("name","source2"); \n');
fprintf(runout,'set("polarization angle",90); \n');
fprintf(runout,'set("phase",90); \n');

%Transmission detector
fprintf(runout,'addprofile; \n');
fprintf(runout,'set("name","transmission"); \n');
setXY(runout,options);
fprintf(runout,'set("z", %e); \n', (structureThickness/2+0.5)*1e-6); %transmission detector is 2 microns above structure
fprintf(runout,'set("override global monitor settings",1); \n');
fprintf(runout,'set("use source limits",0); \n');
fprintf(runout,'set("minimum wavelength",%e); \n',options.fdtd.minMeasWL);
fprintf(runout,'set("maximum wavelength",%e); \n',options.fdtd.maxMeasWL);
fprintf(runout,'set("frequency points",%i); \n',options.fdtd.numMeasWL);
fprintf(runout,'set("use linear wavelength spacing",0); \n');

%Reflection detector
fprintf(runout,'addprofile; \n');  
fprintf(runout,'set("name","reflection"); \n');
setXY(runout,options);
fprintf(runout,'set("z", %e); \n', (-structureThickness/2-0.5)*1e-6); %reflection detector is 2 microns below structure
fprintf(runout,'set("override global monitor settings",1); \n');
fprintf(runout,'set("use source limits",0); \n');
fprintf(runout,'set("minimum wavelength",%e); \n',options.fdtd.minMeasWL);
fprintf(runout,'set("maximum wavelength",%e); \n',options.fdtd.maxMeasWL);
fprintf(runout,'set("frequency points",%i); \n', options.fdtd.numMeasWL);
fprintf(runout,'set("use linear wavelength spacing",0); \n');


if options.fdtd.SU8matrix
    fprintf(runout,'addrect; \n');
    setXY(runout,options);
    fprintf(runout,'set("z",%e); \n', 0);
    fprintf(runout,'set("z span",%e); \n', structureThickness*1e-6);
    fprintf(runout,'set("material","<Object defined dielectric>"); \n');
    fprintf(runout,'set("index",1.58); \n');
    fprintf(runout,'set("override mesh order from material database", 1); \n');
    fprintf(runout,'set("mesh order",3); \n');
end


if options.fdtd.addCubesDirectly
    fprintf(runout,'redrawoff; \n');
    
    
    dx = options.dimensions(1)/Nx*1e-6;
    if strcmp(options.lattice,'square')
        dy = options.dimensions(2)/Ny*1e-6;
    elseif strcmp(options.lattice,'hexagonal')
        dy = options.dimensions(2)*2/Ny*1e-6;
    end
    dz = options.dimensions(3)*options.fdtd.shrinkFactor/Nz*1e-6;
    %initial rect that will be copied
    fprintf(runout,'addrect; set("name","base"); set("x span",%e); set("y span",%e); set("z span",%e); \n',...
        dx, dy, dz);
    %fprintf(runout,'set("material","Au (Gold) - CRC"); \n');
    fprintf(runout,'set("material","Ag (Silver) - CRC"); \n');

    if strcmp(options.lattice,'square')
        fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -options.dimensions(1)/2*1e-6 - dx/2, -options.dimensions(2)/2*1e-6 - dy/2, -structureThickness/2*1e-6 - dz/2);
    elseif strcmp(options.lattice,'hexagonal')
        fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -options.dimensions(1)/2*1e-6 - dx/2, -options.dimensions(2)*1e-6 - dy/2, -structureThickness/2*1e-6 - dz/2);

    end
    %fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -options.dimensions(1)/2*1e-6 - dx/2, -2*options.dimensions(2)/2*1e-6 - dy/2, -param.Z_T_shrunk/2*1e-6 - dz/2 - param.Z_T_shrunk*1e-6*(options.Nstructrepeat-1)/2);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz*options.fdtd.repeatUnits
                if structSkin(i,j,k) == 2
                    fprintf(runout,'select("base"); copy(%e,%e,%e); set("name","%s"); \n',i*dx,j*dy,k*dz, 'C');
                end
            end
        end
    end
    fprintf(runout,'select("base"); delete; \n');
    
    
%%Do same rectangle generation with an etch material    
    %initial rect that will be copied
    fprintf(runout,'addrect; set("name","baseEtch"); set("x span",%e); set("y span",%e); set("z span",%e); \n',...
        dx, dy, dz);
    fprintf(runout,'set("material","etch"); \n');

    if strcmp(options.lattice,'square')
        fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -options.dimensions(1)/2*1e-6 - dx/2, -options.dimensions(2)/2*1e-6 - dy/2, -structureThickness/2*1e-6 - dz/2);
    elseif strcmp(options.lattice,'hexagonal')
        fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -options.dimensions(1)/2*1e-6 - dx/2, -options.dimensions(2)*1e-6 - dy/2, -structureThickness/2*1e-6 - dz/2);

    end
    %fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -options.dimensions(1)/2*1e-6 - dx/2, -2*options.dimensions(2)/2*1e-6 - dy/2, -param.Z_T_shrunk/2*1e-6 - dz/2 - param.Z_T_shrunk*1e-6*(options.Nstructrepeat-1)/2);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz*options.fdtd.repeatUnits
                if structSkin(i,j,k) == 1
                    fprintf(runout,'select("baseEtch"); copy(%e,%e,%e); set("name","%s"); \n',i*dx,j*dy,k*dz, 'E');
                end
            end
        end
    end
    fprintf(runout,'select("baseEtch"); delete; \n');
    
    %make copies
    %fprintf(runout,'select("C");');
    %fprintf(runout,'for (n=2:%i){ \n',options.Nstructrepeat);
    %fprintf(runout,'copy(%e,%e,%e); \n',0,0,options.Z_T_shrunk*1e-6);
    %fprintf(runout,'} \n');
else
    %NOT CURRENTLY SUPPORTED
    fprintf(runout,'addimport; \n');
    fprintf(runout,'importnk("%s","%s",%e,%e,%e,%i); \n',[options.dir,options.currentNKfile],...
        'microns',0,0,0*1e-6,0);
    %fprintf(runout,'for (n=2:%i){ \n',options.Nstructrepeat);
    %    fprintf(runout,'copy(%e,%e,%e); \n',0,0,param.Z_T_shrunk*1e-6);
    %fprintf(runout,'} \n');
end

fprintf(runout,'save("%s"); \n',[options.dir,options.currentLumSave]);
%if options.runSim
    fprintf(runout,'runparallel; \n');
%end
fprintf(runout,'transmission_left = transmission("transmission"); \n');
fprintf(runout,'reflection_left = transmission("reflection"); \n');

fprintf(runout,'switchtolayout; \n');
fprintf(runout,'select("source2"); \n');
fprintf(runout,'set("phase",-90); \n');  %Set to right-handed polarization
%if options.runSim
    fprintf(runout,'runparallel; \n');
%end
fprintf(runout,'transmission_right = transmission("transmission"); \n');
fprintf(runout,'reflection_right = transmission("reflection"); \n');

fprintf(runout,'matlabsavelegacy("%s",transmission_right,reflection_right,transmission_left,reflection_left); \n',[options.dir,options.currentLumResultsFile]);

fprintf(runout,'exit; \n');
fclose(runout);


end

%Sets the x-y dimensions to match the periodicity
function setXY(runout,options)
    fprintf(runout,'set("x", 0); \n');
    fprintf(runout,'set("y", 0); \n');
    if strcmp(options.lattice,'square')
        fprintf(runout,'set("x span", %e ); \n', options.dimensions(1)*1e-6);
        fprintf(runout,'set("y span", %e ); \n', options.dimensions(2)*1e-6);
    elseif strcmp(options.lattice,'hexagonal')
        fprintf(runout,'set("x span", %e ); \n', options.dimensions(1)*1e-6);
        fprintf(runout,'set("y span", %e ); \n', options.dimensions(2)*2*1e-6);        
    end

end
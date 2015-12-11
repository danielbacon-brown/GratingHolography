function writeLumericalRunFileSquare( GAoptions, structureSU8,  param)
%This writes the Lumerical Script that imports the structure, runs the
%simulation, and exports the data.

runout = fopen([GAoptions.dir,GAoptions.LumRunScript], 'w');

%SETUP
fprintf(runout,'newproject; \n');

%Mesh
fprintf(runout,'addfdtd; \n');
fprintf(runout,'set("x", 0); \n');
fprintf(runout,'set("y", 0); \n');
fprintf(runout,'set("z", 0); \n');
fprintf(runout,'set("x span", %e ); \n', GAoptions.dimensions(1)*1e-6);
%if GAoptions.hexagonalGrating
%    fprintf(runout,'set("y span", %e ); \n', 2*GAoptions.dimensions(2)*1e-6);
%else    %square grating
    fprintf(runout,'set("y span", %e ); \n', GAoptions.dimensions(1)*1e-6);
%end
fprintf(runout,'set("simulation time",%i); \n', GAoptions.simulationTime*1e-15);
%fprintf(runout,'set("z span", %e ); \n', (GAoptions.Nstructrepeat*param.Z_T_shrunk+5)*1e-6); %want 2.5 microns on both top and bottom of structure
fprintf(runout,'set("z span", %e ); \n', (param.su8thicknessShrunk+5)*1e-6);
fprintf(runout,'set("mesh type","auto non-uniform"); \n');
fprintf(runout,'set("mesh accuracy",%i); \n',GAoptions.meshAccuracy);
fprintf(runout,'set("use early shutoff",0); \n');
fprintf(runout,'set("x min bc","periodic"); \n');
fprintf(runout,'set("y min bc","periodic"); \n');

%Plane Source
fprintf(runout,'addplane; \n');
fprintf(runout,'set("name","source1"); \n');
fprintf(runout,'set("x", 0); \n');
fprintf(runout,'set("y", 0); \n');
fprintf(runout,'set("x span", %e ); \n', GAoptions.dimensions(1)*1e-6);
if GAoptions.hexagonalGrating
    fprintf(runout,'set("y span", %e ); \n', 2*GAoptions.dimensions(2)*1e-6);
else %square grating
    fprintf(runout,'set("y span", %e ); \n', GAoptions.dimensions(1)*1e-6);
end
%fprintf(runout,'set("z", %e); \n', (-GAoptions.Nstructrepeat*param.Z_T_shrunk/2-1)*1e-6); %plane source is 1 micron below bottom of structure
fprintf(runout,'set("z", %e); \n', (-param.su8thicknessShrunk/2-1)*1e-6);
fprintf(runout,'set("wavelength start",%e); \n', GAoptions.minSourceWL);
fprintf(runout,'set("wavelength stop", %e); \n', GAoptions.maxSourceWL);

%Plane Source for Circular Polarization  %Right-handed (I think)
fprintf(runout,'copy(0,0,0); \n');
fprintf(runout,'set("name","source2"); \n');
fprintf(runout,'set("polarization angle",90); \n');
fprintf(runout,'set("phase",90); \n');

%Transmission detector
fprintf(runout,'addprofile; \n');
fprintf(runout,'set("name","transmission"); \n');
fprintf(runout,'set("x",%e); \n', 0);
fprintf(runout,'set("y",%e); \n', 0);
fprintf(runout,'set("x span",%e); \n', GAoptions.dimensions(1)*1e-6);
if GAoptions.hexagonalGrating
    fprintf(runout,'set("y span", %e ); \n', 2*GAoptions.dimensions(2)*1e-6);
else
    fprintf(runout,'set("y span",%e); \n', GAoptions.dimensions(1)*1e-6);
end
fprintf(runout,'set("z", %e); \n', (param.su8thicknessShrunk/2+2)*1e-6); %transmission detector is 2 microns above structure
fprintf(runout,'set("override global monitor settings",1); \n');
fprintf(runout,'set("use source limits",0); \n');
fprintf(runout,'set("minimum wavelength",%e); \n',GAoptions.minMeasWL);
fprintf(runout,'set("maximum wavelength",%e); \n',GAoptions.maxMeasWL);
fprintf(runout,'set("frequency points",%i); \n',GAoptions.numMeasWL);
fprintf(runout,'set("use linear wavelength spacing",0); \n');

%Reflection detector
fprintf(runout,'addprofile; \n');  
fprintf(runout,'set("name","reflection"); \n');
fprintf(runout,'set("x",%e); \n', 0);
fprintf(runout,'set("y",%e); \n', 0);
fprintf(runout,'set("x span",%e); \n', GAoptions.dimensions(1)*1e-6);
if GAoptions.hexagonalGrating
    fprintf(runout,'set("y span", %e ); \n', 2*GAoptions.dimensions(2)*1e-6);
else
    fprintf(runout,'set("y span",%e); \n', GAoptions.dimensions(1)*1e-6);
end
fprintf(runout,'set("z", %e); \n', (-param.su8thicknessShrunk/2-2)*1e-6); %reflection detector is 2 microns below structure
fprintf(runout,'set("override global monitor settings",1); \n');
fprintf(runout,'set("use source limits",0); \n');
fprintf(runout,'set("minimum wavelength",%e); \n',GAoptions.minMeasWL);
fprintf(runout,'set("maximum wavelength",%e); \n',GAoptions.maxMeasWL);
fprintf(runout,'set("frequency points",%i); \n', GAoptions.numMeasWL);
fprintf(runout,'set("use linear wavelength spacing",1); \n');



if GAoptions.addCubesDirectly
    fprintf(runout,'redrawoff; \n');
    %ASSUMES HEXAGONAL
    dx = GAoptions.dimensions(1)/param.Ncellx*1e-6;
    dy = 2*GAoptions.dimensions(2)/param.Ncelly*1e-6;
    dz = param.su8thicknessShrunk/param.Ncellz*1e-6;
    %initial rect that will be copied
    fprintf(runout,'addrect; set("name","base"); set("x span",%e); set("y span",%e); set("z span",%e); \n',...
        dx, dy, dz);
    %fprintf(runout,'setmaterial("Au(Gold)-CRC");');
    %fprintf(runout,'set("material","Au (Gold) - CRC"); \n');
    fprintf(runout,'set("material","Ag (Silver) - CRC"); \n');
%%MODIFY FOLLOWING LINE?????
    fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -GAoptions.dimensions(1)/2*1e-6 - dx/2, -2*GAoptions.dimensions(2)/2*1e-6 - dy/2, -param.su8thicknessShrunk/2*1e-6 - dz/2);
    %fprintf(runout,'set("x",%e); set("y",%e); set("z",%e); \n', -GAoptions.dimensions(1)/2*1e-6 - dx/2, -2*GAoptions.dimensions(2)/2*1e-6 - dy/2, -param.Z_T_shrunk/2*1e-6 - dz/2 - param.Z_T_shrunk*1e-6*(GAoptions.Nstructrepeat-1)/2);
    for i = 1:param.Ncellx
        for j = 1:param.Ncelly
            for k = 1:param.Ncellz
                if ~structureSU8(i,j,k)
                    fprintf(runout,'select("base"); copy(%e,%e,%e); set("name","%s"); \n',i*dx,j*dy,k*dz, 'C');
                end
            end
        end
    end
    fprintf(runout,'select("base"); delete; \n');
    %make copies
    %fprintf(runout,'select("C");');
    %fprintf(runout,'for (n=2:%i){ \n',GAoptions.Nstructrepeat);
    %fprintf(runout,'copy(%e,%e,%e); \n',0,0,GAoptions.Z_T_shrunk*1e-6);
    %fprintf(runout,'} \n');
else
    fprintf(runout,'addimport; \n');
    fprintf(runout,'importnk("%s","%s",%e,%e,%e,%i); \n',[GAoptions.dir,GAoptions.currentNKfile],...
        'microns',0,0,0*1e-6,0);
    %fprintf(runout,'for (n=2:%i){ \n',GAoptions.Nstructrepeat);
    %    fprintf(runout,'copy(%e,%e,%e); \n',0,0,param.Z_T_shrunk*1e-6);
    %fprintf(runout,'} \n');
end

fprintf(runout,'save("%s"); \n',[GAoptions.dir,GAoptions.currentLumSave]);
if GAoptions.runSim
    fprintf(runout,'runparallel; \n');
end
fprintf(runout,'transmission_left = transmission("transmission"); \n');
fprintf(runout,'reflection_left = transmission("reflection"); \n');

fprintf(runout,'switchtolayout; \n');
fprintf(runout,'select("source2"); \n');
fprintf(runout,'set("phase",-90); \n');  %Set to right-handed polarization
if GAoptions.runSim
    fprintf(runout,'runparallel; \n');
end
fprintf(runout,'transmission_right = transmission("transmission"); \n');
fprintf(runout,'reflection_right = transmission("reflection"); \n');

fprintf(runout,'matlabsave("%s",transmission_right,reflection_right,transmission_left,reflection_left); \n',[GAoptions.dir,GAoptions.currentLumResultsFile]);

fprintf(runout,'exit; \n');
fclose(runout);


end


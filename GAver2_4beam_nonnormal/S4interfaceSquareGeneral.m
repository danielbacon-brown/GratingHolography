classdef S4interfaceSquareGeneral
    %Writes a script for creating a square grating and illumination conditions in S4
    %The grating is a reflection grating - no need for metallic coating due
    %to TIR at SU8 air interface
    
    %ASSUMES SQUARE GRATING
    %ASSUMES only first order diffraction (7 orders for hexagonal)
    
    properties
        cells;
        dimensions;
        %n_metal;
        %k_metal;
        gratingCoatingThickness;
        %isAirGap;
        setBasicScript;
        %setMaterialsScript;
        %setFrontLayersScript;
        %setBackLayersScript;
        collectDataScript;
        getAmplitudesScript;
        c;
        n_interference;
        eps_0;
    end
    
    methods
        
        function S = S4interfaceSquareGeneral(options) %.dimensions,cells,metal_name,thickness_metal, isAirGap,n_interference)
            
            S.c = 2.99e8; %m/s
            S.n_interference = 1.58;
            S.eps_0 = 8.854e-12; %vacuum permittivity %F/m 
            
            %S.isAirGap = options.isAirGap;
            
            S.cells = options.cells;
            S.dimensions = options.dimensions;
%             S.gratingCoatingThickness = options.gratingCoatingThickness;
%             switch options.gratingCoatingMetal
%                 case 'silver'
%                     S.n_metal = 0.054007; %at 532nm
%                     S.k_metal = 3.429;
%                 case 'gold'
%                     S.n_metal = 0.54386; %at 532nm
%                     S.k_metal = 2.2309;
%             end
            
            %Write string for first constant part of file %Constant
            %dimensions, etc

            periodX = options.dimensions(1);
            periodY = options.dimensions(2);
            
            
            %cs = [cs, 'time1 = os.clock() \r\n'];
            
            if strcmp(options.lattice,'square')
                latticeStr = sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',periodX,0,0,periodY)
            elseif strcmp(options.lattice,'hexagonal')
                %latticeStr = sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',-periodX/2,-periodY,periodX/2,-periodY)
                %latticeStr = sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',periodX,0,-periodX/2,periodY)
                latticeStr = sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',periodX,0,periodX/2,periodY)
            end
            
            S.setBasicScript = [ ...
                'S = S4.NewSimulation() \r\n' ...
                ...%Define periodicity
                ..., sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',periodX,0,0,periodY) ...
                , latticeStr ...
                , 'S:SetNumG(9) \r\n' ...
..., 'S:SetNumG(100) \r\n' ...
                , 'S:SetResolution(16) \r\n' ...
                , 'S:UseDiscretizedEpsilon(1) \r\n' ...
                ...%fprintf(ef, 'S:UseLanczosSmoothing(1) \n');
                , 'S:UsePolarizationDecomposition(1) \r\n' ...
                , 'S:UseJonesVectorBasis(1) \r\n' ...
                , 'S:SetLatticeTruncation(''Parallelogramic'') \r\n' ...
                ];
            
            
            
            %cs2 = [cs2, 'S:SetFrequency(1/0.532) \r\n']; %532nm
            %cs2 = [cs2, 'S:SetFrequency(n_glass/0.532) \r\n'];
            %fprintf(ef, 'S:UseNormalVectorBasis(1) \n');
            
            %Define materials
            %setMaterialsScript = '';
            %n_TCO = 1.94; %ITO  %ITO: 1.94 at 532nm   --AZO: 1.89 at 532nm   --FTO: ~2
            %k_TCO = 0.046;
            %n_glass = options.n_glass; %MAKE VARIABLE?
            
%             S.setMaterialsScript = [ ...
%                 'S:AddMaterial("Vacuum", {1,0}) \r\n' ...
%                 , 'S:AddMaterial("SU8", {1.58^2,0}) \r\n' ...
%                 , sprintf('n_TCO = %1.3f      \r\n',n_TCO) ...
%                 , sprintf('k_TCO = %1.3f  --at 532nm \r\n',k_TCO) ...
%                 , 'S:AddMaterial("TCO", {n_TCO^2-k_TCO^2,2*n_TCO*k_TCO}) \r\n' ...
%                 ..., sprintf('n_metal = %1.3f  \r\n',S.n_metal) ...
%                 ..., sprintf('k_metal = %1.3f  \r\n',S.k_metal) ...
%                 ..., 'S:AddMaterial("Metal", {n_metal^2-k_metal^2,2*n_metal*k_metal}) \r\n' ...
%                 , sprintf('n_glass = %1.3f \r\n',n_glass) ...
%                 , 'S:AddMaterial("Glass", {n_glass^2,0}) \r\n' ...
%                 , sprintf('n_prism = %1.3f \r\n',options.n_prism) ...
%                 , 'S:AddMaterial("Prism", {n_prism^2,0}) \r\n' ...
%                 ];
            
%             %Set Materials:
%             S.setMaterialsScript = '';
%             for i_m = 1:length(options.materials)
%                 [n,k] = options.material(i_m).getRI(materialChromosome);
%                 k = imag(options.materialRI(i_m));
%                 S.setMaterialsScript = [ S.setMaterialsScript, ...
%                     sprintf('S:AddMaterial(''%s'', {%2.5f,%2.5f}) \r\n', ...
%                     options.materialNames{i_m}, ...
%                     n^2-k^2, ...
%                     2*n*k ) ...
%                     ];
%             end
            
            
%             %For reflection grating(2), the layer order is: prism, glass
%             %substrat (should be the same), ITO, SU8, SU8-air grating, air 
%             S.setFrontLayersScript = [ ...
%                 'S:AddLayer(''Front'', 0, ''Prism'')  \r\n' ... 
%                 , 'S:AddLayer(''TCOLayer'',0.1, ''TCO'') \r\n' ...
%                 ... , 'S:AddLayer(''PrInterference'', 5, ''SU8'')  -- thick SU8 layer \r\n' ...
%                 ];
%             S.setBackLayersScript = [ ...
%                 'S:AddLayer(''Back'', 0, ''Vacuum'')  \r\n' ... 
%                 ];
            
            %Second constant part of script
            %cs2 = '';
            
            %Define SU8 layer thickness
            %cs2 = [cs2, 'S:AddLayer(''PrInterference'', 5, ''SU8'')  -- thick SU8 layer \r\n'];
            %cs2 = [cs2, 'S:AddLayer(''ITOlayer'', 0.01, ''ITO'') \r\n']; %low reflection thickness
            %           cs2 = [cs2, 'S:AddLayer(''SU8afterGrating'',3,''SU8'')  \r\n'];
            %           cs2 = [cs2, 'S:AddLayer(''Back'',0,''Glass'') -- substrate \r\n']; %glass substrate
            %cs2 = [cs2, 'S:AddLayer(''Back'',0,''SU8'') \n');  %SU8 substrate (in case you want to see the no-reflection case)
            
            
            
            
            
            
            %cs2 = [cs2, 'os.remove(''',GAoptions.dir,'HelixFieldVolOutput.E'') \r\n'];
            %cs2 = [cs2, 'os.remove(''',GAoptions.dir,'HelixFieldVolOutput.H,'') \r\n'];
            
            
            %cs2 = [cs2, 'Cx = math.floor(0.5334/0.02) \r\n'];
            %cs2 = [cs2, 'Cy = math.floor((0.5334*math.sqrt(3)/2)/0.02) \r\n'];
            %            cs2 = [cs2, 'Cx = ',cells(1),' \r\n'];
            %            cs2 = [cs2, 'Cy = ',cells(2),' \r\n'];
            %            cs2 = [cs2, 'Cz = ',cells(3),' \r\n'];
            %            cs2 = [cs2, 'print(Cx) \r\n'];
            %            cs2 = [cs2, 'print(Cy) \r\n'];
            %            cs2 = [cs2, 'print(Cz) \r\n'];
            
            z_offset = 2; %um %distance from air-side of grating to start measuring field
            
            %cs3 = '';
            %cs2 = [cs2, 'z_table = matrix.linspace(', num2str(z_offset),',',num2str(z_offset+dimensions(3) ),',',num2str(dimensions(3)),') \r\n'];
            

            S.collectDataScript = [
                sprintf( 'for z_i = 1, %i do \r\n',S.cells(3)) ...
                , sprintf('	S:GetFieldPlane(z_i * %1.7f + %1.7f, {%i,%i}, ''FileAppend'', dataFilename) \r\n',S.dimensions(3)/S.cells(3), z_offset, S.cells(1),S.cells(2) ) ...
                ...%cs2 = [cs2, 'for z=', num2str(z_offset) ,',' num2str(z_offset+dimensions(3) ),',',num2str(dimensions(3)/cells(3)) ,' do \r\n'];
                ...%cs2 = [cs2, '	S:GetFieldPlane(z, {',num2str(cells(1)),',',num2str(cells(2)),'}, ''FileAppend'', dataFilename) \r\n'];  %Uses the filename to be set in varScript
                ...%cs2 = [cs2, '	print( z_i*',num2str(dimensions(3)/cells(3)),'+',num2str(z_offset),') \r\n'];
                , 'end \r\n' ...
                , 'nG = S:GetNumG() \r\n' ...
                , 'print(nG) \r\n' ...
                ];
            
            
            
            %IF YOU WANT TO COLLECT MANY UNIT CELLS IN Z:
            % S.collectDataScript = [
            %                 sprintf( 'for z_i = 1, %i do \r\n',S.cells(3)) ...
            %                 , sprintf('	S:GetFieldPlane(z_i * %1.7f + %1.7f, {%i,%i}, ''FileAppend'', dataFilename) \r\n',S.dimensions(3)*4/S.cells(3), z_offset, S.cells(1),S.cells(2) ) ...
            %                 , 'end \r\n'];
            
            
            
            
            % %Get field in cartesian coordinate array
            % fprintf(ef, 'for z=(0.064 +1e-6),(0.064+0.532*2 -1e-6),Cz do \n'); %z
            %
            %     fprintf(ef, 'for y=(0 +1e-6),(periodY -1e-6),(periodY/Cy) do \n'); %y
            %         fprintf(ef, 'for x=(0 +1e-6),(periodX -1e-6),(periodX/Cx) do \n'); %y
            %
            %
            %         fprintf(ef, 'end \n');
            %     fprintf(ef, 'end \n');
            % fprintf(ef, 'end \n');
            
            
            %fprintf(ef, 'i = S:GetDiffractionOrder(1, 1) \n');
            %fprintf(ef, 'print(i) \n');
            
            %cs2 = [cs2, 'time2 = os.clock() \n' );
            S.getAmplitudesScript = [ ...
                'forw,back = S:GetAmplitudes(''PrInterference'',2.5) \n' ...
                , 'print(''forward waves:'') \n'...
                , 'for key,value in pairs(forw) do print(key, forw[key][1], forw[key][2]) end \n'...
                , 'print(''backward waves:'') \n'...
                , 'for key,value in pairs(back) do print(key, back[key][1], back[key][2]) end \n'...
                , 'print(''Num G: '') \n' ...
                , 'ng = S:GetNumG() \n' ...
                , 'print(ng) \n'...
                , 'print(''G list: '') \n' ...
                , 'glist = S:GetGList() \n'...
                , 'for key,value in pairs(glist) do print(key, glist[key][1], glist[key][2]) end \n'...
                , 'S:OutputStructurePOVRay(''HelixPOVrayScript.pov'') \n ' ...
                , 'Gu,Gv = S:GetReciprocalLattice() \n' ...
                , 'print(''Gu '', Gu[1][1]*2*math.pi,'' '',Gu[1][2]*2*math.pi, ''     Gv '',Gu[2][1]*2*math.pi,'' '',Gu[2][2]*2*math.pi ) \n' ... [1],'' '',Gv[2]) \n' ...
                ];
            %cs2 = [cs2, 'time3 = os.clock() \n' );
            
            
            
            
            
            %cs2 = [cs2, 'print(time2-time1) \n' );
            %cs2 = [cs2, 'print(time3-time2) \n' );
            %'automatedS4script.lua'
            %fprintf(ef, 'print(''each forward mode:'') \n');
            %fprintf(ef, 'for key,value in pairs(forw) do print(forw[key][1], forw[key][2]) end \n');
            %fprintf(ef, 'for key,value in pairs(forw) do print( forw[key][1]^2+ forw[key][2]^2 ) end \n');
            %fprintf(ef, 'print(''each backward mode:'') \n');
            %fprintf(ef, 'for key,value in pairs(back) do print(back[key][1], back[key][2]) end \n');
            
            %fclose(ef);
            
            
            
            
            
            
            
        end
        
        function intensityDist = doRCWA(S,GAoptions,grating,incidentFieldParams,layerChromosomes,materialChromosomes)
            
            %Use worker ID to different files intended for different
            %workers
            %t = getCurrentTask();
            %disp(['Task# = ', t.ID])
            
            %Since getCurrentTask fails, generate a random integer that
            %will be the file ID -> very low chance of attempting
            %simultaneous file-writes.
            t = num2str(randi(intmax()));
            
            if ~isempty(t)
                dataFilename = sprintf( 'fieldData_%s', t )
                scriptFilename = sprintf( 'automatedS4script_%s.lua', t )
            else
                dataFilename = 'fieldData'
                scriptFilename = 'automatedS4script.lua'
            end
            makeRunScript(S,GAoptions,grating,incidentFieldParams, dataFilename,scriptFilename,layerChromosomes,materialChromosomes); %Make script
            
            if strcmp(GAoptions.hostname,'Daniel-netbook')
                system(['C:/Users/daniel/S4-1.1.1-win32/S4 ', GAoptions.dir,scriptFilename]); %Run script
            elseif strcmp(GAoptions.hostname,'berzerk')
                message = ['running: ', '~/S4/build/S4 ', GAoptions.dir,scriptFilename]
                system(['~/S4/build/S4 ', GAoptions.dir,scriptFilename]);
            end
            disp(['importing data: ', GAoptions.dir,dataFilename,'.E'])
            if ~exist([GAoptions.dir,dataFilename,'.E'],'file') %If you can't find the file, ignore it and move on
                intensityDist = [];
                return;
            end
            A = importdata([GAoptions.dir,dataFilename,'.E']); %Load data from script
            delete([GAoptions.dir,dataFilename,'.E']); %Clear data file for reuse
            delete([GAoptions.dir,dataFilename,'.H']);
%delete([scriptFilename]);
            
            
            Ex = A(:,4) + 1i*A(:,5);
            Ey = A(:,6) + 1i*A(:,7);
            Ez = A(:,8) + 1i*A(:,9);
            I_linear = (Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez)) * S.c*S.n_interference*S.eps_0/2;
            I = reshape(I_linear, S.cells(2),S.cells(1),[]);
            intensityDist = permute(I,[2,1,3]);
            
            
        end
        
        %This adds an arbitrary field to the interference pattern
        function intensityDist = doRCWA_ExtraField(S,GAoptions,grating,incidentFieldParams,layerChromosomes,materialChromosomes,extraFieldParam)
            
            %Use worker ID to different files intended for different
            %workers
            %t = getCurrentTask();
            %disp(['Task# = ', t.ID])
            
            %Since getCurrentTask fails, generate a random integer that
            %will be the file ID -> very low chance of attempting
            %simultaneous file-writes.
            t = num2str(randi(intmax()));
            
            if ~isempty(t)
                dataFilename = sprintf( 'fieldData_%s', t )
                scriptFilename = sprintf( 'automatedS4script_%s.lua', t )
            else
                dataFilename = 'fieldData'
                scriptFilename = 'automatedS4script.lua'
            end
            makeRunScript(S,GAoptions,grating,incidentFieldParams, dataFilename,scriptFilename,layerChromosomes,materialChromosomes); %Make script
            
            if strcmp(GAoptions.hostname,'Daniel-netbook')
                system(['C:/Users/daniel/S4-1.1.1-win32/S4 ', GAoptions.dir,scriptFilename]); %Run script
            elseif strcmp(GAoptions.hostname,'berzerk')
                message = ['running: ', '~/S4/build/S4 ', GAoptions.dir,scriptFilename]
                system(['~/S4/build/S4 ', GAoptions.dir,scriptFilename]);
            end
            disp(['importing data: ', GAoptions.dir,dataFilename,'.E'])
            if ~exist([GAoptions.dir,dataFilename,'.E'],'file') %If you can't find the file, ignore it and move on
                intensityDist = [];
                return;
            end
            A = importdata([GAoptions.dir,dataFilename,'.E']); %Load data from script
            delete([GAoptions.dir,dataFilename,'.E']); %Clear data file for reuse
            delete([GAoptions.dir,dataFilename,'.H']);
delete([scriptFilename]);
            
            

%FINISH
            if strcmp(GAoptions.lattice,'square')
                %Create a 3d position matrix
                ticksX = linspace(0,S.dimensions(1)*(S.cells(1)-1)/S.cells(1), S.cells(1)); %um
                ticksY = linspace(0,S.dimensions(2)*(S.cells(2)-1)/S.cells(1), S.cells(1));
                ticksZ = linspace(0,S.dimensions(3)*(S.cells(3)-1)/S.cells(1), S.cells(1));
                [coorX,coorY,coorZ] = ndgrid(ticksX,ticksY,ticksZ)
                
                eikr_a = zeros(cells(1),cells(2),cells(3),8,8)
                
                for i_i = 1:7  %For each combination of vectors
                    for i_j = 1:7
                        eikr_x = exp(1i.*( k(1,i_i) - k(1,i_j) ).*coorX );  %exp(i*(k_ix-k_jx)*r_x)
                        eikr_y = exp(1i.*( k(2,i_i) - k(2,i_j) ).*coorY );
                        eikr_z = exp(1i.*( k(3,i_i) - k(3,i_j) ).*coorZ );
                        eikr_a(:,:,:,i_i,i_j) = eikr_x.*eikr_y.*eikr_z;   %exp(i*(k_i-k_z)*r)

                    end
                end
                
            end
            
            Ex = A(:,4) + 1i*A(:,5);
            Ey = A(:,6) + 1i*A(:,7);
            Ez = A(:,8) + 1i*A(:,9);
            I_linear = (Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez)) * S.c*S.n_interference*S.eps_0/2;
            I = reshape(I_linear, S.cells(2),S.cells(1),[]);
            intensityDist = permute(I,[2,1,3]);
            
            
        end
        
        
        function makeRunScript(S,GAoptions,grating,incidentFieldParams,dataFilename,scriptFilename,layerChromosomes,materialChromosomes)
            strat = grating.stratum{1,1};
            
            periodX = grating.d21;
            periodY = grating.d32;
            

  
            
            
            %Set Layers:
            setLayerScript = '';
            for i_l = 1:length(GAoptions.S4interfaceOptions.layers)
                setLayerScript = [ setLayerScript, ...
                    sprintf('S:AddLayer(''%s'', %2.3f, ''%s'') \r\n', ...
                    GAoptions.S4interfaceOptions.layers(i_l).layerName, ...
                    GAoptions.S4interfaceOptions.layers(i_l).getThickness( layerChromosomes{i_l} ), ...
                    GAoptions.S4interfaceOptions.layers(i_l).materialName) ...
                    ];
            end

            %Set Materials:
            setMaterialsScript = '';
            for i_m = 1:length(GAoptions.S4interfaceOptions.materials)
                [n,k] = GAoptions.S4interfaceOptions.materials(i_m).getRI( materialChromosomes{i_m} );
                setMaterialsScript = [ setMaterialsScript, ...
                    sprintf('S:AddMaterial(''%s'', {%2.5f,%2.5f}) \r\n', ...
                        GAoptions.S4interfaceOptions.materials(i_m).materialName, ...
                        n^2-k^2, ...
                        2*n*k ) ...
                    ];
            end
            
            
            figure
            hold on
            axis([-periodX, periodX*2, -periodY, periodY*2])
            
            if strcmp(GAoptions.lattice,'square')
                %Define block dimensions and positions:
                %           vertOffset = -periodX;  %Need to shift the block positions up and to the right as you move up in stripe#
                %           horiOffset = -periodY;
                lastY = 0; %marks the position of the end of the last stripe
                for s = 1:length(strat.stripe) %iterate by stripe
                    stripe = strat.stripe{1,s};
                    stripeEnd = stripe.c1; %top edge of stripe
                    lastX = 0; %marks the position of the end of the last block
                    for b = 1:length( stripe.block ) %iterate by block
                        blockEnd = stripe.block{1,b}.c2; %right edge of current block
                        if stripe.block{1,b}.pmt_index == 2 %if it's marked as SU8
                            %                       centerX = ((lastX + blockEnd)/2 + lastY/2) * periodX + horiOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
                            centerX = (lastX + blockEnd)/2 * periodX
                            widthX = (blockEnd - lastX) * periodX
                            %                       centerY = (stripeEnd + lastY)/2 * periodY + vertOffset;
                            centerY = (stripeEnd + lastY)/2 * periodY
                            widthY = (stripeEnd - lastY) * periodY
                            %varScript = [varScript, 'S:SetLayerPatternRectangle(''Grating'', ''SU8'', {', num2str(centerX),',',num2str(centerY),'}, 0, {',num2str(widthX/2),',',num2str(widthY/2),'}) \r\n'];
                            setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];
                            %varScript = [varScript, 'S:SetLayerPatternRectangle(''SU8AirGrating'', ''SU8'', {', num2str(centerX),',',num2str(centerY),'}, 0, {',num2str(widthX/2),',',num2str(widthY/2),'}) \r\n'];
                            %layerScript = [layerScript, 'S:SetLayerPatternRectangle(''SU8MetalGrating'', ''SU8'', {', num2str(centerX),',',num2str(centerY),'}, 0, {',num2str(widthX/2),',',num2str(widthY/2),'}) \r\n'];

                        rectangle('Position',[centerX-widthX/2, centerY-widthY/2, widthX, widthY] )
                        
                        
                        
                        end
                        lastX = blockEnd;
                    end
                    lastY = stripeEnd;
                end
            
            elseif strcmp(GAoptions.lattice,'hexagonal')
                %Define block dimensions and positions:
                vertOffset = -periodX/2;  %Need to shift the block positions up and to the right as you move up in stripe#
                horiOffset = -periodY/2;
                lastY = 0; %marks the position of the end of the last stripe
                for s = 1:length(strat.stripe) %iterate by stripe
                    stripe = strat.stripe{1,s};
                    stripeEnd = stripe.c1; %top edge of stripe
                    lastX = 0; %marks the position of the end of the last block
                    for b = 1:length( stripe.block ) %iterate by block
                        blockEnd = stripe.block{1,b}.c2; %right edge of current block
                        if stripe.block{1,b}.pmt_index == 2 %if it's marked as SU8
                            centerX = ((lastX + blockEnd)/2 + lastY/2) * periodX + horiOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
                            %centerX = (lastX + blockEnd)/2 * periodX;
                            widthX = (blockEnd - lastX) * periodX;
                            centerY = (stripeEnd + lastY)/2 * periodY + vertOffset;
                            %centerY = (stripeEnd + lastY)/2 * periodY;
                            widthY = (stripeEnd - lastY) * periodY;
                            %varScript = [varScript, 'S:SetLayerPatternRectangle(''Grating'', ''SU8'', {', num2str(centerX),',',num2str(centerY),'}, 0, {',num2str(widthX/2),',',num2str(widthY/2),'}) \r\n'];
                            setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];
                            %varScript = [varScript, 'S:SetLayerPatternRectangle(''SU8AirGrating'', ''SU8'', {', num2str(centerX),',',num2str(centerY),'}, 0, {',num2str(widthX/2),',',num2str(widthY/2),'}) \r\n'];
                            %layerScript = [layerScript, 'S:SetLayerPatternRectangle(''SU8MetalGrating'', ''SU8'', {', num2str(centerX),',',num2str(centerY),'}, 0, {',num2str(widthX/2),',',num2str(widthY/2),'}) \r\n'];

                           rectangle('Position',[centerX-widthX/2, centerY-widthY/2, widthX, widthY] )
                        
                        
                        
                        end
                        lastX = blockEnd;
                    end
                    lastY = stripeEnd;
                    %horiOffset = lastY/2
                end                
                
            end
            
            
            %Set incident light parameters
            Es = incidentFieldParams.Esp(2,1);
            Ep = incidentFieldParams.Esp(1,1);
            
            
            %Set the output data filename variable to be used in cs2
            setDataFilenameScript = sprintf('dataFilename = "%s" \r\n',GAoptions.dir,dataFilename);
            
            %excitationScript = '';
            excitationScript = [ ...
                'S:SetExcitationPlanewave( \r\n' ...
                , sprintf('	{%3.8f,%3.8f},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360]) \r\n',incidentFieldParams.phi,incidentFieldParams.theta) ...
                , sprintf('	{%8.8f,%3.8f},  -- s-polarization amplitude and phase (in degrees) \r\n', abs(Es),-1*angle(Es)*180/pi ) ...
                , sprintf('	{%8.8f,%3.8f})  -- p-polarization amplitude and phase \r\n',abs(Ep),-1*angle(Ep)*180/pi ) ...
                , ' \r\n' ...
                , 'S:SetFrequency(1/0.532) \r\n' ...
                ];
            %For visualizing the grating
            %varScript = [varScript, 'S:OutputStructurePOVRay(''HelixPOVrayScript.pov'') \n'];
            
            
            %Write strings to file and run:
            ef=fopen([GAoptions.dir,scriptFilename],'w');  %open/create file
            %fprintf(ef,[S.setBasicScript,S.setMaterialsScript,setLayerScript,excitationScript, setDataFilenameScript, S.collectDataScript]);
            fprintf(ef,[S.setBasicScript,setMaterialsScript,setLayerScript,excitationScript, setDataFilenameScript, S.collectDataScript, S.getAmplitudesScript]);
            %fprintf(ef,[S.constScript1,varScript,S.constScript2,excitationScript,S.constScript3]);
            fclose(ef);
            
            
        end
        
        
        
        
        
    end
    
end
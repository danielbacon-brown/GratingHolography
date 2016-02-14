classdef S4interfaceSquareGeneral
    
    
    properties
        cells;
        dimensions;
        gratingCoatingThickness;
        setBasicScript;
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
            
            
            S.cells = options.cells;
            S.dimensions = options.dimensions;

            


            periodX = options.dimensions(1);
            periodY = options.dimensions(2);
            
            
            %Define periodicity
            if strcmp(options.lattice,'square')
                latticeStr = sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',periodX,0,0,periodY);
            elseif strcmp(options.lattice,'hexagonal')
                latticeStr = sprintf( 'S:SetLattice({%1.5f,%1.5f},{%1.5f,%1.5f}) \r\n',periodX,0,periodX/2,periodY);
            end
            
            S.setBasicScript = [ ...
                'S = S4.NewSimulation() \r\n' ...
                , latticeStr ...
                , 'S:SetNumG(9) \r\n' ...
..., 'S:SetNumG(100) \r\n' ...
                , 'S:UsePolarizationDecomposition(1) \r\n' ...
                , 'S:SetResolution(64) \r\n' ...
                ..., 'S:UseDiscretizedEpsilon(1) \r\n' ...
                ..., 'S:UseJonesVectorBasis(1) \r\n' ...
                ..., 'S:SetLatticeTruncation(''Parallelogramic'') \r\n' ...
                ];
            
            
            
            
            
           
            
            z_offset = 2; %um %distance from air-side of grating to start measuring field
            

            S.collectDataScript = [
                sprintf( 'for z_i = 1, %i do \r\n',S.cells(3)) ...
                , sprintf('	S:GetFieldPlane(z_i * %1.7f + %1.7f, {%i,%i}, ''FileAppend'', dataFilename) \r\n',S.dimensions(3)/S.cells(3), z_offset, S.cells(1),S.cells(2) ) ...
                ...%cs2 = [cs2, 'for z=', num2str(z_offset) ,',' num2str(z_offset+dimensions(3) ),',',num2str(dimensions(3)/cells(3)) ,' do \r\n'];
                ...%cs2 = [cs2, '	S:GetFieldPlane(z, {',num2str(cells(1)),',',num2str(cells(2)),'}, ''FileAppend'', dataFilename) \r\n'];  %Uses the filename to be set in varScript
                ...%cs2 = [cs2, '	print( z_i*',num2str(dimensions(3)/cells(3)),'+',num2str(z_offset),') \r\n'];
                , 'end \r\n' ...
                ..., 'nG = S:GetNumG() \r\n' ...
                ..., 'print(nG) \r\n' ...
                ];
            
            
            
            %IF YOU WANT TO COLLECT MANY UNIT CELLS IN Z:
            if options.repeatingUnits ~= 1
                S.collectDataScript = [
                            sprintf( 'for z_i = 1, %i do \r\n',S.cells(3)) ...
                            , sprintf('	S:GetFieldPlane(z_i * %1.7f + %1.7f, {%i,%i}, ''FileAppend'', dataFilename) \r\n',S.dimensions(3)*options.repeatingUnits/S.cells(3), z_offset, S.cells(1),S.cells(2) ) ...
                            , 'end \r\n'];
            end
            
            
            
            S.getAmplitudesScript = '';  %Don't want to bother getting the amplitudes unless doing a check
            
% KEEP FOR DEBUGGING DIFFRACTION ORDERS            
%             S.getAmplitudesScript = [ ...
%                 'forw,back = S:GetAmplitudes(''PrInterference'',2.5) \n' ...
%                 , 'print(''forward waves:'') \n'...
%                 , 'for key,value in pairs(forw) do print(key, forw[key][1], forw[key][2]) end \n'...
%                 , 'print(''backward waves:'') \n'...
%                 , 'for key,value in pairs(back) do print(key, back[key][1], back[key][2]) end \n'...
%                 , 'print(''Num G: '') \n' ...
%                 , 'ng = S:GetNumG() \n' ...
%                 , 'print(ng) \n'...
%                 , 'print(''G list: '') \n' ...
%                 , 'glist = S:GetGList() \n'...
%                 , 'for key,value in pairs(glist) do print(key, glist[key][1], glist[key][2]) end \n'...
%                 , 'S:OutputStructurePOVRay(''HelixPOVrayScript.pov'') \n ' ...
%                 , 'Gu,Gv = S:GetReciprocalLattice() \n' ...
%                 , 'print(''Gu '', Gu[1][1]*2*math.pi,'' '',Gu[1][2]*2*math.pi, ''     Gv '',Gu[2][1]*2*math.pi,'' '',Gu[2][2]*2*math.pi ) \n' ... [1],'' '',Gv[2]) \n' ...
%                 ];



            
            
            
          
            
            
            
            
            
            
            
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
                dataFilename = sprintf( 'fieldData_%s', t );
                scriptFilename = sprintf( 'automatedS4script_%s.lua', t );
            else
                dataFilename = 'fieldData';
                scriptFilename = 'automatedS4script.lua';
            end
            makeRunScript(S,GAoptions,grating,incidentFieldParams, dataFilename,scriptFilename,layerChromosomes,materialChromosomes); %Make script
            
            if strcmp(GAoptions.hostname,'Daniel-netbook')
                system(['C:/Users/daniel/S4-1.1.1-win32/S4 ', GAoptions.dir,scriptFilename]); %Run script
            elseif strcmp(GAoptions.hostname,'berzerk')
                %message = ['running: ', '~/S4/build/S4 ', GAoptions.dir,scriptFilename];
                system(['~/S4/build/S4 ', GAoptions.dir,scriptFilename]);
            elseif strcmp(GAoptions.hostname,'lotus-bud')
                system(['~/S4mod/build/S4 ', GAoptions.dir,scriptFilename]);
            end
            %disp(['importing data: ', GAoptions.dir,dataFilename,'.E'])
            if ~exist([GAoptions.dir,dataFilename,'.E'],'file') %If you can't find the file, ignore it and move on
                intensityDist = [];
                return;
            end
            %A = importdata([GAoptions.dir,dataFilename,'.E']); %Load data from script
            A = dlmread([GAoptions.dir,dataFilename,'.E']); %Load data from script
            delete([GAoptions.dir,dataFilename,'.E']); %Clear data file for reuse
            delete([GAoptions.dir,dataFilename,'.H']);
            if GAoptions.runSingle == 0
                delete([scriptFilename]);
            else
                disp(['Saving: ', scriptFilename])
            end
            
            if size(A,2)<9 %If data is full of NaN, skip it and move on
                intensityDist = [];
                return;
            end
            
            Ex = A(:,4) + 1i*A(:,5);
            Ey = A(:,6) + 1i*A(:,7);
            Ez = A(:,8) + 1i*A(:,9);
            I_linear = (Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez)) * S.c*S.n_interference*S.eps_0/2;
%             if strcmp(GAoptions.lattice,'square')
%                 unitCellArea = GAoptions.period*GAoptions.period*1e-12;  %m^2
%             elseif strcmp(GAoptions.lattice,'hexagonal')
%                 unitCellArea = GAoptions.period*GAoptions.period*sqrt(3)/2*1e-12; %m^2
%             end
%             I_linear = (Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez)) * S.c*S.n_interference*S.eps_0/2 * unitCellArea;  %Need to multiply by unit cell area in m^2 to get I in W/m^2
%Avoiding scaling by unit cell area, E is in V/m, so I is in W/m^2
            I = reshape(I_linear, S.cells(2),S.cells(1),[]);
            intensityDist = permute(I,[2,1,3]);
            
            
        end
%         
%         %This adds an arbitrary field to the interference pattern
%         function intensityDist = doRCWA_ExtraField(S,GAoptions,grating,incidentFieldParams,layerChromosomes,materialChromosomes,extraFieldParam)
%             
%             %Use worker ID to different files intended for different
%             %workers
%             %t = getCurrentTask();
%             %disp(['Task# = ', t.ID])
%             
%             %Since getCurrentTask fails, generate a random integer that
%             %will be the file ID -> very low chance of attempting
%             %simultaneous file-writes.
%             t = num2str(randi(intmax()));
%             
%             if ~isempty(t)
%                 dataFilename = sprintf( 'fieldData_%s', t );
%                 scriptFilename = sprintf( 'automatedS4script_%s.lua', t );
%             else
%                 dataFilename = 'fieldData';
%                 scriptFilename = 'automatedS4script.lua';
%             end
%             makeRunScript(S,GAoptions,grating,incidentFieldParams, dataFilename,scriptFilename,layerChromosomes,materialChromosomes); %Make script
%             
%             if strcmp(GAoptions.hostname,'Daniel-netbook')
%                 system(['C:/Users/daniel/S4-1.1.1-win32/S4 ', GAoptions.dir,scriptFilename]); %Run script
%             elseif strcmp(GAoptions.hostname,'berzerk')
%                 message = ['running: ', '~/S4/build/S4 ', GAoptions.dir,scriptFilename];
%                 system(['~/S4/build/S4 ', GAoptions.dir,scriptFilename]);
%             end
%             disp(['importing data: ', GAoptions.dir,dataFilename,'.E'])
%             if ~exist([GAoptions.dir,dataFilename,'.E'],'file') %If you can't find the file, ignore it and move on
%                 intensityDist = [];
%                 return;
%             end
%             A = importdata([GAoptions.dir,dataFilename,'.E']); %Load data from script
%             delete([GAoptions.dir,dataFilename,'.E']); %Clear data file for reuse
%             delete([GAoptions.dir,dataFilename,'.H']);
% delete([scriptFilename]);
%             
%             if size(A,2)<9 %If data is full of NaN, skip it and move on
%                 intensityDist = [];
%                 return;
%             end
% 
% %FINISH
%             if strcmp(GAoptions.lattice,'square')
%                 %Create a 3d position matrix
%                 ticksX = linspace(0,S.dimensions(1)*(S.cells(1)-1)/S.cells(1), S.cells(1)); %um
%                 ticksY = linspace(0,S.dimensions(2)*(S.cells(2)-1)/S.cells(1), S.cells(1));
%                 ticksZ = linspace(0,S.dimensions(3)*(S.cells(3)-1)/S.cells(1), S.cells(1));
%                 [coorX,coorY,coorZ] = ndgrid(ticksX,ticksY,ticksZ);
%                 
%                 eikr_a = zeros(S.cells(1),S.cells(2),S.cells(3),8,8);
%                 
%                 for i_i = 1:7  %For each combination of vectors
%                     for i_j = 1:7
%                         eikr_x = exp(1i.*( k(1,i_i) - k(1,i_j) ).*coorX );  %exp(i*(k_ix-k_jx)*r_x)
%                         eikr_y = exp(1i.*( k(2,i_i) - k(2,i_j) ).*coorY );
%                         eikr_z = exp(1i.*( k(3,i_i) - k(3,i_j) ).*coorZ );
%                         eikr_a(:,:,:,i_i,i_j) = eikr_x.*eikr_y.*eikr_z;   %exp(i*(k_i-k_z)*r)
% 
%                     end
%                 end
%                 
%             end
%             
%             Ex = A(:,4) + 1i*A(:,5);
%             Ey = A(:,6) + 1i*A(:,7);
%             Ez = A(:,8) + 1i*A(:,9);
%             I_linear = (Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez)) * S.c*S.n_interference*S.eps_0/2;
%             I = reshape(I_linear, S.cells(2),S.cells(1),[]);
%             intensityDist = permute(I,[2,1,3]);
%             
%             
%         end
        
        
        function makeRunScript(S,GAoptions,grating,incidentFieldParams,dataFilename,scriptFilename,layerChromosomes,materialChromosomes)
            %strat = grating.stratum{1,1};
            
            save('lastgrating','grating');
            
            
            %periodX = grating.d21;
            %periodY = grating.d32;
            periodX = 0.5334;
            periodY = 0.4619;

            
            
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
                            centerX = (lastX + blockEnd)/2 * periodX;
                            widthX = (blockEnd - lastX) * periodX;
                            %                       centerY = (stripeEnd + lastY)/2 * periodY + vertOffset;
                            centerY = (stripeEnd + lastY)/2 * periodY;
                            widthY = (stripeEnd - lastY) * periodY;
                            setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];
                        
                        end
                        lastX = blockEnd;
                    end
                    lastY = stripeEnd;
                end
                
            elseif strcmp(GAoptions.lattice,'hexagonal')
                %                 %Define block dimensions and positions:
                %                 vertOffset = -periodX/2;  %Need to shift the block positions up and to the right as you move up in stripe#
                %                 horiOffset = -periodY/2;
                %                 lastY = 0; %marks the position of the end of the last stripe
                %                 for s = 1:length(strat.stripe) %iterate by stripe
                %                     stripe = strat.stripe{1,s};
                %                     stripeEnd = stripe.c1; %top edge of stripe
                %                     lastX = 0; %marks the position of the end of the last block
                %                     for b = 1:length( stripe.block ) %iterate by block
                %                         blockEnd = stripe.block{1,b}.c2; %right edge of current block
                %                         if stripe.block{1,b}.pmt_index == 2 %if it's marked as SU8
                %                             centerX = ((lastX + blockEnd)/2 + lastY/2) * periodX + horiOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
                %                             %centerX = (lastX + blockEnd)/2 * periodX;
                %                             widthX = (blockEnd - lastX) * periodX;
                %                             centerY = (stripeEnd + lastY)/2 * periodY + vertOffset;
                %                             %centerY = (stripeEnd + lastY)/2 * periodY;
                %                             widthY = (stripeEnd - lastY) * periodY;
                %                             setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];
                %                         end
                %                         lastX = blockEnd;
                %                     end
                %                     lastY = stripeEnd;
                %                 end
                
%                 %Define block dimensions and positions:
%                 %horiOffset = 0;
%                 vertOffset = 0;
%                 horiOffset = -periodX/2;  %Need to shift the block positions up and to the right as you move up in stripe#
%                 %vertOffset = -periodY/2;
%                 lastX = 0; %marks the position of the end of the last stripe
%                 for s = 1:length(strat.stripe) %iterate by stripe
%                     stripe = strat.stripe{1,s};
%                     stripeEnd = stripe.c1; %top edge of stripe
%                     lastY = 0; %marks the position of the end of the last block
%                     for b = 1:length( stripe.block ) %iterate by block
%                         blockEnd = stripe.block{1,b}.c2; %right edge of current block
%                         if stripe.block{1,b}.pmt_index == 2 %if it's marked as SU8
%                             %centerY = ((lastY + blockEnd)/2 + lastY/2) * periodX + horiOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
%                             centerY = (lastY + blockEnd)/2 * periodY + vertOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
%                             %centerX = (lastX + blockEnd)/2 * periodX;
%                             widthY = (blockEnd - lastY) * periodY;
%                             centerX = (stripeEnd + lastX)/2 * periodX + horiOffset;
%                             %centerY = (stripeEnd + lastY)/2 * periodY;
%                             widthX = (stripeEnd - lastX) * periodX;
%                             setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];
%                             
%                             %rectangle('Position',[centerX-widthX/2, centerY-widthY/2, widthX,widthY]);
%                             
%                         end
%                         lastX = blockEnd;
%                     end
%                     lastX = stripeEnd;
%                 end
                
                
            end


%TEST
grating.CX1 = 0            
grating.CX2 = 0
grating.CX3 = 0
grating.CX4 = 0
grating.CY1 = 0
grating.CY2 = 0
grating.CY3 = 0
grating.CY4 = 0


%Need to define the grating as a polygon
%(Layer,Material,center,angle,vertices)
%vertices= [x1,y1,x2,y2,x3,y3 ...

setLayerScript = [setLayerScript, 'S:SetLayerPatternPolygon(''Grating'', ''SU8'',{-0.2667,-.2309},0, {',...
            sprintf( ' %1.8f, %1.8f, ', 0, 0), ...
            sprintf( ' %1.8f, %1.8f, ', 0.5334, 0), ...
            sprintf( ' %1.8f, %1.8f, ', 0.5334, grating.L4 + grating.CY4), ...
            sprintf( ' %1.8f, %1.8f, ', 0.5334-grating.CX4, grating.L4), ...
            sprintf( ' %1.8f, %1.8f, ', grating.W1+grating.W2 + grating.CX3, grating.L3), ...
            sprintf( ' %1.8f, %1.8f, ', grating.W1+grating.W2, grating.L3+grating.CY3), ...
            sprintf( ' %1.8f, %1.8f, ', grating.W1+grating.W2, grating.L2-grating.CY2), ...
            sprintf( ' %1.8f, %1.8f, ', grating.W1+grating.W2-grating.CX2, grating.L2), ...
            sprintf( ' %1.8f, %1.8f, ', grating.CX1, grating.L1), ...
            sprintf( ' %1.8f, %1.8f, ', 0, grating.L1-grating.CY1), ...
            '} ) \r\n']
            
            
            %Set incident light parameters
            Es = incidentFieldParams.Esp(2,1);
            Ep = incidentFieldParams.Esp(1,1);
            
%%%%%%TEST%%%%%%%            
Es = 3150.2
Ep = 293.7

            
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
            displayScript = ['S:OutputStructurePOVRay(''HelixPOVrayScript.pov'') \r\n'];
            displayScript = [displayScript,'S:OutputLayerPatternDescription(''Grating'',''patternDescription.ps'') \r\n'];
            displayScript = [displayScript,'S:OutputLayerPatternRealization(''Grating'',25,25,''patternRealization'') \r\n'];
            
            
            
            %Write strings to file and run:
            ef=fopen([GAoptions.dir,scriptFilename],'w');  %open/create file
            %fprintf(ef,[S.setBasicScript,setMaterialsScript,setLayerScript,excitationScript, setDataFilenameScript, S.collectDataScript, S.getAmplitudesScript]);
            fprintf(ef,[S.setBasicScript,setMaterialsScript,setLayerScript,excitationScript, setDataFilenameScript, S.collectDataScript, S.getAmplitudesScript,displayScript]);
            fclose(ef);
            
            
        end
        
        
        
        
        
    end
    
end
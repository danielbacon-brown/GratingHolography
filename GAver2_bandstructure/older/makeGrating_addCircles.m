function grating = makeGrating_addCircles( chromosomalData, GAoptions )
%This creates the grating in GDC format.


    clear stratum stripe stratum1 stratum2 stripe1 stripe2 block1 block2
    
    d = GAoptions.period; %grating period (unitless)
    Nx = GAoptions.Ngratingx;
    Ny = GAoptions.Ngratingy;
    Nlev = GAoptions.Nlev;
    
    ucell = zeros(Nx,Ny);
    
    %onesbase = ones(Nx,Ny)
    
    for iter = 1:GAoptions.Ncircles
        [circ_r,circ_c] = meshgrid(linspace(0,GAoptions.period,Nx),linspace(0,GAoptions.period,Ny)); %Creates coordinate list of all points
        C = sqrt((circ_r-chromosomalData.circ_x(iter)).^2+(circ_c-chromosomalData.circ_y(iter)).^2)<=GAoptions.circleradius; %Selects coordinates within a circular region
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)+d).^2+(circ_c-chromosomalData.circ_y(iter)+0).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)-d).^2+(circ_c-chromosomalData.circ_y(iter)+0).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)+0).^2+(circ_c-chromosomalData.circ_y(iter)+d).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)+0).^2+(circ_c-chromosomalData.circ_y(iter)-d).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)+d).^2+(circ_c-chromosomalData.circ_y(iter)+d).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)+d).^2+(circ_c-chromosomalData.circ_y(iter)-d).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)-d).^2+(circ_c-chromosomalData.circ_y(iter)+d).^2)<=GAoptions.circleradius;
        C = C | sqrt((circ_r-chromosomalData.circ_x(iter)-d).^2+(circ_c-chromosomalData.circ_y(iter)-d).^2)<=GAoptions.circleradius;
        ucell = ucell | C;
    end
    
    grating = GAoptions.constantGrating;  %gets grating with parameters common to all gratings

    if GAoptions.useDielectric == 1   %%%If there is a dielectric%%%
%         qc = 1;%current structure-layer
%         %SETTING UP STRATA
%         stratum{2*Nlev+1}.type = 2;  %preallocates array
%         for q=1:(2*Nlev)+1 %need two levels within each level because of dielectric
%             stratum{q}.type=2;%biperiodic **
%             stratum{q}.h11=1;%stratum periodicity is same as grating periodicity
%             stratum{q}.h12=0;
%             stratum{q}.h21=0;
%             stratum{q}.h22=1;    
%             if (rem(q,2))  %if q is odd
%                 stratum{q}.thick = chromosomalData.thick(Nlev+1);  %thickness of lay is equal to thickness of dielectric
%             else  %if q is even
%                 stratum{q}.thick = chromosomalData.thick(qc) - chromosomalData.thick(Nlev+1);  %thickness =  layer thickness - dielectric thickness
%                 qc = qc + 1;  %number of structure-layers completed
%             end      
%         end
% 
%         %SETTING RELIEF PROFILE
%         %1 = air; 2 = SU8; 3 = TiO2
%         %lowest layer
%         q = 1;
%         ut = (chromosomalData.ucell >= q);  %ut is equal to 1 where the corresponding ucell is atleast 1
%         for i = 1:GAoptions.Ngratingx
%             stratum{q}.stripe{i}.type = 1;
%             stratum{q}.stripe{i}.c1 = i/Nx;  %thickness of stripe
%             for j = 1:Ny
%                 stratum{q}.stripe{i}.block{j}.c2 = j/Ny;  %length of block
%                 stratum{q}.stripe{i}.block{j}.pmt_index =  3 - 1*(ut(i,j));  %2 if filled, 3 if empty
%             end
%         end 
% 
% 
%         for q = 1:Nlev;
%             %upper half of layer q and lower half of q+1
%             ut = (chromosomalData.ucell >= q);  %equal to 1 if the height is at least height of current layer
%             utm = (chromosomalData.ucell >= (q-1));  %equal to 1 if the height is at least height of (current layer - 1)  %not used
%             utp = (chromosomalData.ucell >= (q+1));  %equal to 1 if the height is at least height of (current layer + 1)
%             %Upper half of layer q
%             for i = 1:Nx
%                 stratum{(2*q)}.stripe{i}.type = 1;  %inhomogeneous stripe
%                 stratum{(2*q)}.stripe{i}.c1 = i/Nx;  %width of stripe
%                 for j = 1:Ny
%                     stratum{(2*q)}.stripe{i}.block{j}.c2 = j/Ny;  %length of block
%                     stratum{(2*q)}.stripe{i}.block{j}.pmt_index = 1 +(ut(i,j));  %2 if filled and 1 if empty
%                 end
%             end  
%             %Lower half of (q+1)
%             for i = 1:Nx
%                 stratum{(2*q) +1}.stripe{i}.type = 1;
%                 stratum{(2*q) +1}.stripe{i}.c1 = i/Nx;
%                 for j = 1:Ny
%                     stratum{(2*q) +1}.stripe{i}.block{j}.c2 = j/Ny;
%                     stratum{(2*q) +1}.stripe{i}.block{j}.pmt_index =  0*(utm(i,j)) + 2*(ut(i,j)) -1*(utp(i,j)) + 1 ;  %3 if filled and above is empty;  2 if filled and above is filled; 1 if empty
%                 end
%             end  
%         end
% 
%         for q = 1:(2*GAoptions.Nlev)+1  %copy the strata data to grating
%             grating.stratum{q} = stratum{q};
%         end
        
    else %%%%if there is no dielectric%%%
        
        stratum{Nlev}.type=2; %preallocates space for strata data
        for q=1:(Nlev) %need two levels withing each level because of dielectric
            stratum{q}.type=2;%biperiodic **
            stratum{q}.h11=0;%stratum periodicity is same as grating periodicity
            stratum{q}.h12=1;
            stratum{q}.h21=1;
            stratum{q}.h22=0;    
            stratum{q}.thick = chromosomalData.thick(q);  %thickness =  layer thickness
        end   
        for q = 1:Nlev;
            for i = 1:Nx
                stratum{q}.stripe{i}.type = 1;  %inhomogeneous stripe
                stratum{q}.stripe{i}.c1 = i/Nx;  %width of stripe
                for j = 1:Ny
                    stratum{q}.stripe{i}.block{j}.c2 = j/Ny;  %length of block
                    %A block is filled at the nth level and below if there
                    %are n '1's in that column of the chromosome.
                    if ucell(i,j)>=q
                        stratum{q}.stripe{i}.block{j}.pmt_index = GAoptions.pmtIndices.SU8;
                    else
                        stratum{q}.stripe{i}.block{j}.pmt_index = GAoptions.pmtIndices.air;
                    end
                end
            end
            grating.stratum{q} = stratum{q};
        end
           
    end    


end


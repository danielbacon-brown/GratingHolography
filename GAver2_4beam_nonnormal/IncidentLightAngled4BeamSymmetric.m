classdef IncidentLightAngled4BeamSymmetric %< IncidentLightGenerator
    %Abstract class for generating a grating consistent with GDCalc
    

   properties (SetAccess = protected, GetAccess = protected)
       incidentField;  %structure for GDcalc
       chromNpsi; %Number of chromosome positions for orientation angle
       chromNchi; %Number of chromosome positions for ellipticity angle
       n_interference; %Refractive index of material light is interfering in (photoresist)
       n_incidence; %Refractive index of material light is coming from
       wavelength;
       phi; %azimuthal angle of incidence
       theta; %polar angle of incidence
       Pdens; %incident power %W/m^2 
       eps_0; %vacuum permittivity %F/m
       c; %Speed of light
   end
   
   properties % Public access
   end 
   
   methods 
       
       function L = IncidentLightAngled4BeamSymmetric(options)
            L.wavelength = options.wavelength;
            L.n_interference = options.n_interference;
            L.n_incidence = options.n_incidence;
            %k_x = pi/options.period;
            k_0 = 2*pi/L.wavelength   %um^-1
            k_x = k_0/(sqrt(2 + 1/(options.C_over_A^2)))
            %k_y = pi/options.period;
            k_y=k_x
            k_z = (1/options.C_over_A) * k_x
            %L.phi = atan( (1/k_x^2+1/k_y^2)/(1/(k_z*L.n_interference/L.n_incidence)^2) ); %azimuthal angle
            %phi_interference = atan( sqrt(2) * options.C_over_A); %Gives azimuthal angle inside photoresist (relative to normal)
            phi_interference = atan( sqrt(k_x^2+k_y^2)/k_z)
            phi_glass = asin(L.n_interference/L.n_incidence*sin(phi_interference) ) %Gives azimuthal angle in glass substrate
            L.phi = phi_glass*180/pi %S$ wants angles in degrees
            L.theta = 45; %polar angle %relative to x-axis of grating
           
            %F.incidentField.wavelength = options.wavelength; %um %GDC will divide by refractive index of superstrate
            %F.incidentField.f2=sin(phi0)*cos(theta0)/options.wavelength; %These are projections of the incident spatial-frequency vector on x-y plane (in this case, zero)
            %F.incidentField.f3=sin(phi0)*sin(theta0)/options.wavelength;
            
            L.Pdens = options.beamPowerDens * options.period*options.period*1e-12; %W/m^2 * (um*um) * (m^2/um^2) %Gives power going through a unit cell (if normal) 
            L.eps_0 = 8.854e-12; %F/m
            L.c = 2.99e8; %m/s
            
            L.chromNpsi = options.chromNpsi;
            L.chromNchi = options.chromNchi;
       end
       
       %chi = ellipticity
       %psi = orientation
       function fieldParameters = generateField(L,chromosomeSection) %Generates a specific incident field and polarization base on chromosome
            
            fieldParameters.phi =  L.phi;
            fieldParameters.theta = L.theta;
            
            [chi,psi] = convertChrom_gc(chromosomeSection,[L.chromNchi,L.chromNpsi]); %Convert chromosome to fraction of min-max values
            fieldParameters.chi = chi*pi/2 - pi/4; %converts to within +-45degrees
            fieldParameters.psi = psi*pi - pi/2; %Converts to within +-90degrees
            
            
            if isfield(L,'Pdens')
                fieldParameters.Esp = [cos(psi)*cos(chi)-1i*sin(psi)*sin(chi); ...  %Calculate incident Field E-vector
                    sin(psi)*cos(chi)+1i*cos(psi)*sin(chi)] .* sqrt(2*L.Pdens/(L.n_incidence*L.eps_0*L.c) ) %/ cos(fieldParameters.phi) ); %Should be units of V/sqrt(unit-area)
            else
                fieldParameters.Esp = [cos(psi)*cos(chi)-1i*sin(psi)*sin(chi); ...  %Calculate incident Field E-vector
                    sin(psi)*cos(chi)+1i*cos(psi)*sin(chi)]  %/ cos(fieldParameters.phi) ); %Should be units of V/sqrt(unit-area)
            end
       end

       function chromosomeSize = getChromosomeSize(L) %Returns size of required chromosome
           chromosomeSize = L.chromNchi + L.chromNpsi;
       end
       
       function plotPolarization(L,chi,psi)  %Plots incident polarization
           
           %Jones vector
           Es = cos(psi)*cos(chi)-1i*sin(psi)*sin(chi);
           Ep = sin(psi)*cos(chi)+1i*cos(psi)*sin(chi);
                    
           t = linspace(0,2*pi,1000); %timesteps
           omega = 1; %relative wavenumber %value is unimportant
           
           ses = Es * exp(1i*-omega*t); %S-component of E-field
           pes = Ep * exp(1i*-omega*t); %P-component of E-field
           
           %plot
           figure
           plot(real(ses),real(pes))
           xlim([-1,1])
           ylim([-1,1])
           axis square
           set(gca,'XTick',[-1,1])
           set(gca,'YTick',[-1,1])
       end
       function plotPolarizationEsp(L,Esp)  %Plots incident polarization
           
           %Jones vector
           Es = Esp(1) 
           Ep = Esp(2) 
                    
           t = linspace(0,2*pi,1000); %timesteps
           omega = 1; %relative wavenumber %value is unimportant
           
           ses = Es * exp(1i*-omega*t); %S-component of E-field
           pes = Ep * exp(1i*-omega*t); %P-component of E-field
           
           %plot
           figure
           plot(real(ses),real(pes))
           xlim([-1,1])
           ylim([-1,1])
           axis square
           set(gca,'XTick',[-1,1])
           set(gca,'YTick',[-1,1])
       end
       
       
   end
   
   
   
   
end
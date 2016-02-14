function Icart = hex2cart(Ihex,cellsCart)
%Converts volumetric hexagonal grid to cartesian grid Icart with trilinear
%interpolation.
%Assumes u=(a,0), v=(sqrt(3)/2*a,a/2,0), w=(0,0,c)


%plotVolume(Ihex,0.5);

NxCart = cellsCart(1);
NyCart = cellsCart(2);
NzCart = cellsCart(3);
NxUV = size(Ihex,1);
NyUV = size(Ihex,2);
NzUV = size(Ihex,3);
%Icart = zeros(NxCart,NyCart,NzCart);
%Xq = zeros(NzCart*NyCart*NxCart,1); %Query points = cartesian grid points
%Yq = zeros(NzCart*NyCart*NxCart,1); 
%Zq = zeros(NzCart*NyCart*NxCart,1);
Xq = zeros(NxCart,NyCart,NzCart); %Query points = cartesian grid points
Yq = zeros(NxCart,NyCart,NzCart); 
Zq = zeros(NxCart,NyCart,NzCart);
cs = linspace(0,1,NzCart+1);  %Take each dimension of unit cell to be 1 (Nx+1 ticks for Nx lengths)
bs = linspace(0,1,NyCart+1);  
as = linspace(0,1,NxCart+1);
%Figure out coordinates of each query point in the UV coordinate system:
for c_i = 1:NzCart  
    for b_i = 1:NyCart
        for a_i = 1:NxCart
            c = cs(c_i);
            b = bs(b_i);
            a = as(a_i);
            
            %a = x; %Position in Cartesian space
            %b = y;
            %e = z;
            %d = b/(NyCart/NyUV); %Get position in U-V space
            %c = (b-d*NxCart/(2*NxUV))/(NxCart/NxUV);
            %f = e/(NzCart/NzUV);
            %d = a;
            
            d = a-b/2; %shift over x-values as y increases  %Should increase x by 1/2 for b=1
            e = b;
            f = c;
            
            if (d < 0)
                d = d + 1;
            end
            
            
%             Xq((c_i-1)*NyCart*NxCart+(b_i-1)*NxCart+a_i) = d;
%             Yq((c_i-1)*NyCart*NxCart+(b_i-1)*NxCart+a_i) = e;
%             Zq((c_i-1)*NyCart*NxCart+(b_i-1)*NxCart+a_i) = f;
            Xq(a_i,b_i,c_i) = d;
            Yq(a_i,b_i,c_i) = e;
            Zq(a_i,b_i,c_i) = f;


        end
    end
end



Xqr = reshape(Xq,[],1);
Yqr = reshape(Yq,[],1);
Zqr = reshape(Zq,[],1);


Xvec = linspace(0,1,NxUV);
Yvec = linspace(0,1,NyUV);
Zvec = linspace(0,1,NzUV);



%REPLICATE INPUT DATA TO MAKE BETTER INTERPOLATION AT PERIODICITY
IhexRep = repmat(Ihex,[3,3,3]);
XvecRep = linspace(-1,2,NxUV*3);
YvecRep = linspace(-1,2,NyUV*3);
ZvecRep = linspace(-1,2,NzUV*3);
Icartlin = interp3(XvecRep,YvecRep,ZvecRep,double(IhexRep),double(Xq),double(Yq),double(Zq)); %interpolates points
Icart = reshape(Icartlin,NxCart,NyCart,NzCart);
Icart = permute(Icart,[2,1,3]);



% %Need to shift Xgrid values up so that is it half shifted by 
% Icartlin = interp3(Xvec,Yvec,Zvec,double(Ihex),double(Xq),double(Yq),double(Zq)); %interpolates points
% Icart = reshape(Icartlin,NxCart,NyCart,NzCart);
% Icart = permute(Icart,[2,1,3]);
% %sizeIcart = size(Icart)


%plotVolume(Icart,0.5);

end
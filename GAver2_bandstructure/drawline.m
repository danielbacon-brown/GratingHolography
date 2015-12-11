function vol_border = drawline(vol, p0, p1, r)

r = round(r);

inc = .002;
[Nx,Ny,Nz] = size(vol);
%p1 = p0 + p1
%IF (~keyword_set(r)) THEN r=Nx/5

%el = sphere(2*r+1) %el is like a ball point pen
el = zeros(2*r+1,2*r+1,2*r+1);
for a = 1:2*r+1
    for b = 1:2*r+1
        for c = 1:2*r+1
            if (a-(r+1))^2 + (b-(r+1))^2 + (c-(r+1))^2 < r^2
                el(a,b,c) = 1;
            end
        end
    end
end
% elout = el;

% r = int32(r);
dx = ceil(Nx/4); %border
dy = ceil(Ny/4);
dz = ceil(Nz/4);
vol_border = zeros(Nx+2*(dx+r)+2,Ny+2*(dy+r)+2,Nz+2*(dz+r)+2);
%vol_border(dx+r+1:dx+r+Nx,dy+r+1:dy+r+Ny,dz+r+1:dz+r+Nz) = vol;
%vol_border = zeros(Nx+2*(dx+r),Ny+2*(dy+r),Nz+2*(dz+r));
vol_border(dx+r+1:dx+r+Nx,dy+r+1:dy+r+Ny,dz+r+1:dz+r+Nz) = vol;
A = p1(1)-p0(1);
B = p1(2)-p0(2);
C = p1(3)-p0(3);

for t=0:inc:1
    x = round(p0(1)+A*t) + dx+r+1;
    y = round(p0(2)+B*t) + dy+r+1;
    z = round(p0(3)+C*t) + dz+r+1;
    %[x-r+1, x+r+1, y-r+1, y+r+1, z-r+1, z+r+1]
%     xout = x;
%     yout = y;
%     zout = z;
%     rout = r;
%     size0 = size(vol_border);
%     size1 = size(vol_border(x-r:x+r,y-r:y+r,z-r:z+r));
%     size2 = size(el);
    vol_border((x-r):(x+r),(y-r):(y+r),(z-r):(z+r)) = vol_border((x-r):(x+r),(y-r):(y+r),(z-r):(z+r)) | el;
end



vol_border = vol_border((dx+r+1):(Nx+dx+r),(dy+r+1):(Ny+dy+r),(dz+r+1):(Nz+dz+r)) > 0;

message = 'line drawn'

end
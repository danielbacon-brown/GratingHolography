function skin = calcSkin(struct, numIterations)
%struct: 0=SU8, 1=void
%skin: 0=SU8, 1=void, 2=skin

%Assumes the skin starts to develop where the value is 1, and propagates
%inward into the 1-regions.
%Used for simple model of electroless deposition
%Sets to 1 at the surface of the structure, 0 elsewhere

skin = zeros(size(struct) + [2,2,2]);  %Creates a skin matrix ( 1 => is part of skin)
newskin = zeros(size(struct) + [2,2,2]);

struct = padarray( struct, [1,1,1] )

%Copy edges of the structure:
struct(1,:,:) = struct(end-1,:,:);
struct(end,:,:) = struct(2,:,:);
struct(:,1,:) = struct(:,end-1,:);
struct(:,end,:) = struct(:,2,:);
struct(:,:,1) = struct(:,:,end-1);
struct(:,:,end) = struct(:,:,2);

for i = 1:numIterations
    
    %     %First copy the edges to the other side for periodicity
    %     lx0 = struct(1,:,:);
    %     lx1 = struct(end,:,:);
    %     struct = cat(1,lx0,struct,lx1);
    %
    %     ly0 = struct(:,1,:);
    %     ly1 = struct(:,end,:);
    %     struct = cat(2,ly0,struct,ly1);
    %
    %     lz0 = struct(:,:,1);
    %     lz1 = struct(:,:,end);
    %     struct = cat(3,lz0,struct,lz1);
    
    %Copy edges of the skin in each iteration
    skin(1,:,:) = skin(end-1,:,:);
    skin(end,:,:) = skin(2,:,:);
    skin(:,1,:) = skin(:,end-1,:);
    skin(:,end,:) = skin(:,2,:);
    skin(:,:,1) = skin(:,:,end-1);
    skin(:,:,end) = skin(:,:,2);
    
    %Set skin to 1 where this point == 1 and  (either an adjacent point == 0 or
    %an adjacent point is a skin point
    %left,right,up,down,forward, or backward.
    
    newskin(2:(end-1),2:(end-1),2:(end-1)) = newskin(2:(end-1),2:(end-1),2:(end-1)) | ( struct(2:(end-1),2:(end-1),2:(end-1)) == 1 & (struct(1:(end-2),2:(end-1),2:(end-1)) == 0  | skin(1:(end-2),2:(end-1),2:(end-1)) == 1 ) );
    newskin(2:(end-1),2:(end-1),2:(end-1)) = newskin(2:(end-1),2:(end-1),2:(end-1)) | ( struct(2:(end-1),2:(end-1),2:(end-1)) == 1 & (struct(3:(end),2:(end-1),2:(end-1))  == 0  | skin(3:(end),2:(end-1),2:(end-1)) == 1 ) );
    newskin(2:(end-1),2:(end-1),2:(end-1)) = newskin(2:(end-1),2:(end-1),2:(end-1)) | ( struct(2:(end-1),2:(end-1),2:(end-1)) == 1 & (struct(2:(end-1),1:(end-2),2:(end-1)) == 0  | skin(2:(end-1),1:(end-2),2:(end-1)) == 1 ) );
    newskin(2:(end-1),2:(end-1),2:(end-1)) = newskin(2:(end-1),2:(end-1),2:(end-1)) | ( struct(2:(end-1),2:(end-1),2:(end-1)) == 1 & (struct(2:(end-1),3:(end),2:(end-1))  == 0  | skin(2:(end-1),3:(end),2:(end-1)) == 1 ) );
    newskin(2:(end-1),2:(end-1),2:(end-1)) = newskin(2:(end-1),2:(end-1),2:(end-1)) | ( struct(2:(end-1),2:(end-1),2:(end-1)) == 1 & (struct(2:(end-1),2:(end-1),1:(end-2)) == 0  | skin(2:(end-1),2:(end-1),1:(end-2)) == 1 ) );
    newskin(2:(end-1),2:(end-1),2:(end-1)) = newskin(2:(end-1),2:(end-1),2:(end-1)) | ( struct(2:(end-1),2:(end-1),2:(end-1)) == 1 & (struct(2:(end-1),2:(end-1),3:(end)) == 0  | skin(2:(end-1),2:(end-1),3:(end)) == 1 ) );

    %Need to wait until all adjacency calcs are done to set these new values
    %(Check whether each cell is adjacent to skin as of last turn)
    skin = newskin;
    
end


%Cut skin to remove excess edges used for periodicity
%skin = skin(2:(end-1), 2:(end-1), 2:(end-1) );

%To get final structure description (0=SU8, 1=void, 2=skin) add the skin
%(1=skin, 0=everything else) and struct (1=void, 0=SU8) matrices:
skin = skin(2:(end-1), 2:(end-1), 2:(end-1) ) + struct(2:(end-1), 2:(end-1), 2:(end-1) );

end
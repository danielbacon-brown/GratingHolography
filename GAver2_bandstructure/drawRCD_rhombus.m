%Draws an rod-connected diamond structure

function model = drawRCD_rhombus(mcells,r)

%unit cell = 25pixels - same as S4 max
%p = 25
p = mcells(1)


%model = zeros(3*p,3*p,3*p) %creates empty cube
model = zeros(p,p,p);

r = p/4  %radius based on sids diagram


%List of points to make type1 points from
%type1 = y-and-above, x-and-below
type1points = [ 1/2, 1/2, 0; ...
                1/2, 1/2, 1; ...
                0, 0, 1/2; ...
                1, 0, 1/2; ...
                0, 1, 1/2; ...
                1, 1, 1/2]*p;
%type2 = x-and above, y-and-below            
type2points = [ 1/2, 0, 1/4; ...
                1/2, 1, 1/4; ...
                0, 1/2, 3/4; ...
                1, 1/2, 3/4]*p;

type1directions = [ 0, 1/4, 1/8; ...
                    0,-1/4, 1/8; ...
                    1/4, 0,-1/8; ...
                   -1/4, 0,-1/8]*p; 
type2directions = [ 0, 1/4,-1/8; ...
                    0,-1/4,-1/8; ...
                    1/4, 0, 1/8; ...
                   -1/4, 0, 1/8]*p; 
               
            
%For each type1 point, draw 4lines in appropriate directions:
for i1p = 1:size(type1points,1)
    for i1d = 1:size(type1directions,1)
        model = drawline(model,type1points(i1p,:), type1points(i1p,:)+type1directions(i1d,:), r);
    end
end
for i2p = 1:size(type2points,1)
    for i2d = 1:size(type2directions,1)
        model = drawline(model,type2points(i2p,:), type2points(i2p,:)+type2directions(i2d,:), r);
    end
end


% %main 1
% model = drawline(model,[0,  0,  0],  [p/4,p/4,p/4],r);
% model = drawline(model,[p/2,p/2,0],  [p/4,p/4,p/4],r);
% model = drawline(model,[0,  p/2,p/2],[p/4,p/4,p/4],r);
% model = drawline(model,[p/2,0,  p/2],[p/4,p/4,p/4],r);
% 
% %main x,y
% model = drawline(model,[1/2*p,1/2*p,0],    [3/4*p,3/4*p,1/4*p],r);
% model = drawline(model,[p,    p,    0],    [3/4*p,3/4*p,1/4*p],r);
% model = drawline(model,[p,    1/2*p,1/2*p],[3/4*p,3/4*p,1/4*p],r);
% model = drawline(model,[1/2*p,p,    1/2*p],[3/4*p,3/4*p,1/4*p],r);
%  
% %main x,z
% model = drawline(model,[1/2*p,0,    1/2*p],[3/4*p,1/4*p,3/4*p],r);
% model = drawline(model,[p,    1/2*p,1/2*p],[3/4*p,1/4*p,3/4*p],r);
% model = drawline(model,[p,    0,    p],    [3/4*p,1/4*p,3/4*p],r);
% model = drawline(model,[1/2*p,1/2*p,p],    [3/4*p,1/4*p,3/4*p],r);
% 
% %main y,z
% model = drawline(model,[0,    1/2*p,1/2*p],[1/4*p,3/4*p,3/4*p],r);
% model = drawline(model,[1/2*p,p,    1/2*p],[1/4*p,3/4*p,3/4*p],r);
% model = drawline(model,[1/2*p,1/2*p,p],    [1/4*p,3/4*p,3/4*p],r);
% model = drawline(model,[0,    p,    p],    [1/4*p,3/4*p,3/4*p],r);
% 
% 
% %corner x
% model = drawline(model,[5/4*p,-1/4*p,-1/4*p],[p,1,1],r);
% model = drawline(model,[3/4*p, 1/4*p,-1/4*p],[p,1,1],r);
% model = drawline(model,[5/4*p, 1/4*p, 1/4*p],[p,1,1],r);
% model = drawline(model,[3/4*p,-1/4*p, 1/4*p],[p,1,1],r);
% 
% %corner y
% model = drawline(model,[-1/4*p,5/4*p,-1/4*p],[1,p,1],r);
% model = drawline(model,[ 1/4*p,3/4*p,-1/4*p],[1,p,1],r);
% model = drawline(model,[-1/4*p,3/4*p, 1/4*p],[1,p,1],r);
% model = drawline(model,[ 1/4*p,5/4*p, 1/4*p],[1,p,1],r);
% 
% %corner z
% model = drawline(model,[ 1/4*p,-1/4*p,3/4*p],[1,1,p],r);
% model = drawline(model,[-1/4*p, 1/4*p,3/4*p],[1,1,p],r);
% model = drawline(model,[ 1/4*p, 1/4*p,5/4*p],[1,1,p],r);
% model = drawline(model,[-1/4*p,-1/4*p,5/4*p],[1,1,p],r);
% 
% %corner x,y,z
% model = drawline(model,[ 3/4*p, 5/4*p,3/4*p],[p,p,p],r);
% model = drawline(model,[ 5/4*p, 3/4*p,3/4*p],[p,p,p],r);
% model = drawline(model,[ 5/4*p, 5/4*p,5/4*p],[p,p,p],r);
% model = drawline(model,[ 3/4*p, 3/4*p,5/4*p],[p,p,p],r);




%model = smooth3(model,'gaussian',[5,5,5])


    figure; 
    patched = patch(isosurface(shave(model),0.5));  %creates surface at threshold
    set(patched,'FaceColor', [255 127 80]/256, 'EdgeColor', 'none');  %display stuff
    view(3);
    axis tight;
    axis equal;
    camlight;
    lighting gouraud;
    xlim([1 p]);
    ylim([1 p]);
    zlim([1 p]);
    box on;
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    
    
    %STL Export:
%    M = isosurface(padarray(model,[1,1,1]),0.5);
%    tr=TriRep(M.faces,M.vertices);
%    stlwrite('perfect_RCD_single.stl',tr.Triangulation,tr.X);


    message = 'RCD made'

end

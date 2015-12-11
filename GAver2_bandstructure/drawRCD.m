%Draws an rod-connected diamond structure

function model = drawRCD()

%unit cell = 100pixels
p = 400


%model = zeros(3*p,3*p,3*p) %creates empty cube
model = zeros(p,p,p);

r = p/10  %radius based on sids diagram

%main 1
model = drawline(model,[0,  0,  0],  [p/4,p/4,p/4],r);
model = drawline(model,[p/2,p/2,0],  [p/4,p/4,p/4],r);
model = drawline(model,[0,  p/2,p/2],[p/4,p/4,p/4],r);
model = drawline(model,[p/2,0,  p/2],[p/4,p/4,p/4],r);

%main x,y
model = drawline(model,[1/2*p,1/2*p,0],    [3/4*p,3/4*p,1/4*p],r);
model = drawline(model,[p,    p,    0],    [3/4*p,3/4*p,1/4*p],r);
model = drawline(model,[p,    1/2*p,1/2*p],[3/4*p,3/4*p,1/4*p],r);
model = drawline(model,[1/2*p,p,    1/2*p],[3/4*p,3/4*p,1/4*p],r);
 
%main x,z
model = drawline(model,[1/2*p,0,    1/2*p],[3/4*p,1/4*p,3/4*p],r);
model = drawline(model,[p,    1/2*p,1/2*p],[3/4*p,1/4*p,3/4*p],r);
model = drawline(model,[p,    0,    p],    [3/4*p,1/4*p,3/4*p],r);
model = drawline(model,[1/2*p,1/2*p,p],    [3/4*p,1/4*p,3/4*p],r);

%main y,z
model = drawline(model,[0,    1/2*p,1/2*p],[1/4*p,3/4*p,3/4*p],r);
model = drawline(model,[1/2*p,p,    1/2*p],[1/4*p,3/4*p,3/4*p],r);
model = drawline(model,[1/2*p,1/2*p,p],    [1/4*p,3/4*p,3/4*p],r);
model = drawline(model,[0,    p,    p],    [1/4*p,3/4*p,3/4*p],r);


%corner x
model = drawline(model,[5/4*p,-1/4*p,-1/4*p],[p,1,1],r);
model = drawline(model,[3/4*p, 1/4*p,-1/4*p],[p,1,1],r);
model = drawline(model,[5/4*p, 1/4*p, 1/4*p],[p,1,1],r);
model = drawline(model,[3/4*p,-1/4*p, 1/4*p],[p,1,1],r);

%corner y
model = drawline(model,[-1/4*p,5/4*p,-1/4*p],[1,p,1],r);
model = drawline(model,[ 1/4*p,3/4*p,-1/4*p],[1,p,1],r);
model = drawline(model,[-1/4*p,3/4*p, 1/4*p],[1,p,1],r);
model = drawline(model,[ 1/4*p,5/4*p, 1/4*p],[1,p,1],r);

%corner z
model = drawline(model,[ 1/4*p,-1/4*p,3/4*p],[1,1,p],r);
model = drawline(model,[-1/4*p, 1/4*p,3/4*p],[1,1,p],r);
model = drawline(model,[ 1/4*p, 1/4*p,5/4*p],[1,1,p],r);
model = drawline(model,[-1/4*p,-1/4*p,5/4*p],[1,1,p],r);

%corner x,y,z
model = drawline(model,[ 3/4*p, 5/4*p,3/4*p],[p,p,p],r);
model = drawline(model,[ 5/4*p, 3/4*p,3/4*p],[p,p,p],r);
model = drawline(model,[ 5/4*p, 5/4*p,5/4*p],[p,p,p],r);
model = drawline(model,[ 3/4*p, 3/4*p,5/4*p],[p,p,p],r);


%NON-ROTATED:
%corner 0
% model = drawline(model, [ 1/4*p, 1/4*p, 


model = smooth3(model,'gaussian',[5,5,5])


%plot

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
    M = isosurface(padarray(model,[1,1,1]),0.5);
    tr=TriRep(M.faces,M.vertices);
    stlwrite('perfect_RCD_single.stl',tr.Triangulation,tr.X);


    message = 'RCD made'

end

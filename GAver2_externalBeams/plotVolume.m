function plotVolume(I,fill)

if exist('fill')
    threshold = fixfill(I,256,fill);
else
    threshold = fixfill(I,256,fill);
end

    figure
    patched = patch(isosurface(padarray(I,[1,1,1],100000),threshold));
    %patched = patch(isosurface(I,threshold));
    set(patched,'FaceColor', [30 255 30]/256, 'EdgeColor', 'none');
    view(3);
    camlight
    axis equal
    lighting gouraud
    xlim([2 size(I,2)+1])
    ylim([2 size(I,1)+1])
    zlim([2 size(I,3)+1])

end
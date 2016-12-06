% scratch. Used to test temporary commands.

fv = isosurface(unit_cell,4);
%fv = isosurface(ts111,8.899e4);
p=patch(fv);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); % axis tight
axis([1 53 1 53 1 73])
% axis equal
camlight 
lighting gouraud
box on


% another support with one layer with zeros on side.

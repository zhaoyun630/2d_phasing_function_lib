function showSupport(support)
support = support>0;
[x,y,z] = ind2sub(size(support), find(support));
plot3(x, y, z, 'k.');
axis equal
% axis off
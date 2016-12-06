% create lattice

diff_pattern = importdata('diff_pattern_512_demo.mat');
[s1,s2] = size(diff_pattern);
r1 = 16;
r2 = 16;
period_pixel_1 = s1/r1;
period_pixel_2 = s2/r2;

a = zeros(period_pixel_1,period_pixel_2);
a(1,1) = 1;
a(period_pixel_1,1)=1;
a(1,period_pixel_2)=1;
a(period_pixel_1,period_pixel_2)=1;

c = repmat(a,[r1,r2]);

h = fspecial('gaussian',[5 5],0.5);
d = imfilter(c,h);
d = d/max(d(:));
for i = 1:10
    d = imfilter(d,h)/max(d(:));
end

e = d>0.1;
save('lattice_demo_oversample_2.mat','e');

imagesc(e)
axis equal
axis off
% Create a sinc function

% not working good.

% yun zhao at Mar 14 2016.


%% It works
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



%% the following code doesn't work

% diff_pattern = importdata('Difpat.mat');
% [s1,s2] = size(diff_pattern);
% n_replica = 100;
% n_resolution_x = 257;
% n_resolution_y = 257;
% 
% lattice = zeros(n_resolution_x,n_resolution_y);
% a1 = linspace(-2,2,n_resolution_x);
% a2 = linspace(-2,2,n_resolution_y);
% [x,y]=meshgrid(a1,a2);
% 
% % c1 = (sin(n_replica*2*pi*x).*sin(n_replica*2*pi*y)).^2+0.0000000000000000000000001;
% % c2 = (sin(2*pi*x).*sin(2*pi*y)).^2+0.0000000000000000000000001;
% 
% c1 = sin(n_replica*2*pi*x).^2;
% c2 = sin(n_replica*2*pi*y).^2;
% 
% d1 = (sin(2*pi*x).*sin(2*pi*y)).^2;
% 
% for i = 1:n_resolution_x
%     for j = 1:n_resolution_y
%         if c1(i,j) ~=0 && c2(i,j) ~=0
%             lattice(i,j) = c1(i,j)*c2(i,j)/d1(i,j);
%         elseif c1(i,j) ==0 && c2(i,j) ==0
%                 lattice(i,j) = n_replica^4;
%         else
%             lattice(i,j) = n_replica^2;
%         end
%     end
% end
% 
% % lattice = c1./d1;
% lattice = lattice/max(lattice(:));
% imshow(lattice);
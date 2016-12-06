% load hkl file and rearrage it into 3D volumn array

 h_max = 30;
 k_max = 30;
 l_max = 30;
 M=dlmread('3ZEKhkl.txt'); % for the hkl file, you may need to delete the 
                           %  first and last line and save it as txt.
 hkl = int8(M(:,1:3)); % save hkl in a new matrix and convert it into integers.
                         % convert float hkl index to integer.
 F = M(:,4);  % extract intensity
 F1 = 100*F/norm(F);
 phase = M(:,4);
 F_comp = F.*exp(-1j.*phase*pi/180);
 % B =100* abs(log(F))/(norm(abs(log(F)))); %log scale
 [L_M, W_M]=size(M);
 
 Rec_lattice = zeros(h_max,k_max,l_max);
 
 for i=1:L_M
     Rec_lattice(hkl(i,1)+1,hkl(i,2)+1,hkl(i,3)+1) = F_comp(i);
 end
 
 obj = abs(ifftn(Rec_lattice));
 
 %[x,y,z] = meshgrid(1:h_max,1:k_max,1:l_max);
% scatter3(hkl(:,1),hkl(:,2),hkl(:,3), F1(:),'b','filled');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % this may work. you need to check details.
 % [xx,yy,zz] = meshgrid(1:30,1:30,1:30);
 % isosurface(xx,yy,zz,Rec_lattice,10);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  2nd way to do 3d plot
 %  scatter3(hkl(:,1),hkl(:,2),hkl(:,3),F2*100,'o','b','filled')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 % one way to do 3d plot
%  
%  v = Rec_lattice;
%  %# visualize the volume
% p = patch( isosurface(v,0) );                 %# create isosurface patch
% isonormals(v, p)                              %# compute and set normals
% set(p, 'FaceColor','r', 'EdgeColor','none')   %# set surface props
% daspect([1 1 1])                              %# axes aspect ratio
% view(3), axis vis3d tight, box on, grid on    %# set axes props
% camproj perspective                           %# use perspective projection
% camlight, lighting phong, alpha(0.8)           %# enable light, set transparency
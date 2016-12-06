% Used to compare diffraction intensity from protein and solvent

%% Input parameter

% Experimental structure
obj = importdata('inputID1_obj.mat');

objThreshold = 0.40; % solvent fraction


objSize = size(obj);
totalPixel = objSize(1)*objSize(2)*objSize(3);

sortDensity = sort(obj(:));
thresholdDensity = sortDensity(round(objThreshold*totalPixel));

% structure factor amplitudes and phases of experimental unit cell.
f_u = fftn(obj);
f_u_amp = abs(f_u);
f_u_ang = angle(f_u);


% support1 mask the area of protein
support1 = obj > thresholdDensity;
protein_region = obj.*support1;

% structure factor amplitudes and phases of protein in vacuum
f_p = fftn(protein_region);
f_p_amp = abs(f_p);
f_p_ang = angle(f_p);

% support2 mask the area of solvent 
support2 = 1 - support1;
solvent_region = obj.*support2;


% structure factor amplitudes and phases of solvent in vacuum
f_s = fftn(solvent_region);
f_s_amp = abs(f_s);
f_s_ang = angle(fftn(f_s));


% difference in structure factor amplitudes between protein and solvent
diff_f_amp = (f_p_amp-f_s_amp)./(f_p_amp+f_s_amp);

figure(1)
hist(diff_f_amp(:),200);


% difference in structure factor phases between protein and solvent
diff_f_ang = (f_p_ang - f_s_ang)*180/pi;
diff_f_cos = cos((f_p_ang - f_s_ang));

figure(2)
% hist(diff_f_ang(:),200);
hist(diff_f_cos(:),200);
% f2 = hist(diff_f_cos(:),200);
% f2 = f2/sum(f2);
% bar(f2);


% % histogram of structure factor of unit cell
% figure(3)
% hist(f_u_ang(:),200);

figure(4)
hist(f_u_amp(:),200);

% % histogram of structure factor of protein region
% figure(5)
% hist(f_u_ang(:),200);
% 
% figure(6)
% hist(f_u_amp(:),200);
% 
% % histogram of structure factor of solvent region
% figure(7)
% hist(f_u_ang(:),200);
% 
% figure(8)
% hist(f_u_amp(:),200);
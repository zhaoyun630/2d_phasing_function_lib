function obj = hkl2obj(hklFile, densityShift)

%% 

% This function is used to get obj from HKL file
% hklFile is the filename. For example: '3rdx_full_3A.hkl'
% densityShift: shift the overall density. Here I suggest to use a value in
% [0,1]
% Note that F(000) is always 0 from hkl file as it couldn't be directly
% measured.
% When densityShift = 0, it means no shift at all. The total charge is 0;
% When densityShift = 1, it means to shift the overall charges so that the
% minimum charge is 0.
% densityShift can be anyvalue in this program, but it is suggested to be
% within [0,1]


%% Part 1 get charge density in one unit cell from hkl file

%  for the hkl file, you may need to delete the
%  first and last line and save it as txt.
M=dlmread(hklFile); 

% save hkl in a new matrix and convert it into integers.
hkl = int32(M(:,1:3)); 

% Save phase
phase = M(:,5);

% convert float hkl index to integer.
F = M(:,4);  % extract intensity

%  F_000 = 50*max(F); % Here I use the maximum F to approaximate F_000.
F_000=0;

% Normalize Structure factors F.
% F1 = 100*F/norm(F);



F_comp = F.*exp(1j.*phase*pi/180);
% B =100* abs(log(F))/(norm(abs(log(F)))); %log scale


L_M=length(M); 

h_max = max(hkl(:,1));
k_max = max(hkl(:,2));
l_max = max(hkl(:,3));

% Creat a 3D array to save structure factors in reciprocal space.
Rec_lattice = zeros(2*h_max+1,2*k_max+1,2*l_max+1);

% Create a mask for futhur R factor calculation
Rec_lattice_mask = zeros(2*h_max+1,2*k_max+1,2*l_max+1);

for i=1:L_M
    h = int32(hkl(i,1)+h_max+1);
    k = int32(hkl(i,2)+k_max+1);
    l = int32(hkl(i,3)+l_max+1);
    Rec_lattice(h,k,l) = F_comp(i);
    Rec_lattice_mask(h,k,l) = 1;
end

% Keep an eye on the value of F000, it may be zero or some other constant.
Rec_lattice(h_max+1,k_max+1,l_max+1)=F_000;

% save('rec_lattice_3A.mat','Rec_lattice');
% save('rec_lattice_mask_3A.mat','Rec_lattice_mask');


%% Part 2 Inverse Fourier Transform to get density map
% Do inverse Fourier transform to get obj
 unit_cell_complex = ifftn(ifftshift(Rec_lattice));%for any FFT, the first number [111] in matrix is always (000) in real space.
 unit_cell_real = fftshift(real(unit_cell_complex));
 
%  unit_cell_imag = fftshift(imag(unit_cell_complex));
%  unit_cell_abs = fftshift(abs(unit_cell_complex));

%% Part 3 Shift the overall charge density
rhoMin = min(unit_cell_real(:));
% rhoMax = max(unit_cell_real(:));

obj = unit_cell_real - rhoMin*densityShift;

% save('unitCellRealShift1Res3A.mat','obj');

%% Part 4 generate triple cell
% Refer to function: getTriCell


%% Part 5 generate support
% Refer to function: getSupport; getTriSupport

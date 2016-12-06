
function [obj1,obj2] = hkl2obj_v3(hklFile1, hklFile2, densityShift)

%%
% This function is used to calculate density map from two hkl files with
% different resolution range. The size of density map will be scaled based
% on the higher resolution hkl file.
% The difference between this version and v2:
% V3 takes two hkl files as an input.
% It will automatically scale the two density map to the same size.

% I think I should use a different name for this function. Using V3 is a
% little bit confusing.    -> delete this line when you make the change.

% Created by Yun Zhao on 2016.4.12
% Last editted by Yun Zhao on 2016.4.12


%% Map hkl list into a 3-D array.
% New in V2:
% Compare with V1, this version take hkl from P1 symmetry. Namely, hkl
% lists half of structure factor amplitudes. 
% In V1, the hkl should include all +/- hkl values.

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
[hkl1,hklMask1] = hklReader(hklFile1);
[hkl2,hklMask2] = hklReader(hklFile2);

[x1,x2,x3] = size(hkl1);
[y1,y2,y3] = size(hkl2);

a1 = ceil(abs((x1-y1)/2))+1;
a2 = floor((x1+y1)/2);

b1 = ceil(abs((x2-y2)/2))+1;
b2 = floor((x2+y2)/2);

c1 = ceil(abs((x3-y3)/2))+1;
c2 = floor((x3+y3)/2);

if (x1+x2+x3)>(y1+y2+y3)
    hkl0 = zeros(x1,x2,x3);   
    hkl0(a1:a2,b1:b2,c1:c2) = hkl2;
    hkl2 = hkl0;
    
    hklMask0 = zeros(x1,x2,x3);
    hklMask0(a1:a2,b1:b2,c1:c2) = hklMask2;
    hklMask2 = hklMask0;
elseif (x1+x2+x3)<(y1+y2+y3)
    hkl0 = zeros(y1,y2,y3);
    hkl0(a1:a2,b1:b2,c1:c2) = hkl1;
    hkl1 = hkl0;
    
    hklMask0 = zeros(x1,x2,x3);
    hklMask0(a1:a2,b1:b2,c1:c2) = hklMask1;
    hklMask1 = hklMask0;
else
end

save('hklMask1.mat','hklMask1');
save('hklMask2.mat','hklMask2');


%% Part 2 Inverse Fourier Transform to get density map
% Do inverse Fourier transform to get obj
unit_cell_complex_1 = ifftn(ifftshift(hkl1));%for any FFT, the first number [111] in matrix is always (000) in real space.
unit_cell_real_1 = real(unit_cell_complex_1);

unit_cell_complex_2 = ifftn(ifftshift(hkl2));%for any FFT, the first number [111] in matrix is always (000) in real space.
unit_cell_real_2 = real(unit_cell_complex_2);
 
 
%  unit_cell_imag = fftshift(imag(unit_cell_complex));
%  unit_cell_abs = fftshift(abs(unit_cell_complex));

%% Part 3 Shift the overall charge density
rhoMin1 = min(unit_cell_real_1(:));
rhoMin2 = min(unit_cell_real_2(:));
% rhoMax = max(unit_cell_real(:));

obj1 = unit_cell_real_1 - rhoMin1*densityShift;
obj2 = unit_cell_real_2 - rhoMin2*densityShift;
% save('unitCellRealShift1Res3A.mat','obj');

%% Part 4 generate triple cell
% Refer to function: getTriCell


%% Part 5 generate support
% Refer to function: getSupport; getTriSupport

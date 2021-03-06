hklFile = '3rdu_full.hkl';

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
F_comp_c = conj(F_comp);
% B =100* abs(log(F))/(norm(abs(log(F)))); %log scale


L_M=length(M); 

h_max = max(hkl(:,1));
k_max = max(hkl(:,2));
l_max = max(hkl(:,3));

% Creat a 3D array to save structure factors in reciprocal space.
hkl = zeros(2*h_max+1,2*k_max+1,2*l_max+1);

% Create a mask for futhur R factor calculation
hklMask = zeros(2*h_max+1,2*k_max+1,2*l_max+1);

for i=1:L_M
    h1 = int32(hkl(i,1)+h_max+1);
    k1 = int32(hkl(i,2)+k_max+1);
    l1 = int32(hkl(i,3)+l_max+1);
    h2 = int32(-hkl(i,1)+h_max+1);
    k2 = int32(-hkl(i,2)+k_max+1);
    l2 = int32(-hkl(i,3)+l_max+1);
    hkl(h1,k1,l1) = F_comp(i);
    hkl(h2,k2,l2) = F_comp_c(i);
    hklMask(h1,k1,l1) = 1;
    hklMask(h2,k2,l2) = 1;
end
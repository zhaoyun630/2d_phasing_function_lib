% This piece of code used to generate 3D electron density contour
% from pdb file.

% First created by Yun Zhao on 04/14/14.
% Last editted by Yun Zhao on 04/16/14.


% In this part, the code extract the coordinate information from pdb.
% Then they are deposited to different element coordinate arrays.
Carbon_coord = [0 0 0];
Oxygen_coord = [0 0 0];
Sulfur_coord = [0 0 0];
Nitrogen_coord = [0 0 0];

% B-facor
Carbon_B_factor = 0;
Oxygen_B_factor = 0;
Sulfur_B_factor = 0;
Nitrogen_B_factor = 0;

pdbid= fopen('5LYZ_coord.pdb');
tline = fgetl(pdbid);

while ischar(tline)
        %disp(tline);
        line_column = strsplit(tline);
        coord = str2double(line_column(7:9));
        B_factor = str2double(line_column(11));
        atom_type = char(line_column(12));
        switch atom_type
            case 'C'
                Carbon_coord = [Carbon_coord; coord]; % append coordinates
                Carbon_B_factor = [Carbon_B_factor; B_factor]; % append B factors
            case 'O'
                Oxygen_coord = [Oxygen_coord; coord];
                Oxygen_B_factor = [Oxygen_B_factor; B_factor];
            case 'N'
                Nitrogen_coord = [Nitrogen_coord; coord];
                Nitrogen_B_factor = [Nitrogen_B_factor; B_factor];
            case 'S'
                Sulfur_coord = [Sulfur_coord; coord];
                Sulfur_B_factor = [Sulfur_B_factor; B_factor];
            otherwise
                sprintf('This line is not included: \n%s\n',tline);
        end
    tline = fgetl(pdbid);
    %pause;
end

fclose(pdbid);

% delete the first zero row in each matrix
% And add a column specifying number of electrons
% Carbon_coord = [Carbon_coord(2:end,:) 6*ones(Carbon_size(1)-1,1)];
% Oxygen_coord = [Oxygen_coord(2:end,:) 8*ones(Oxygen_size(1)-1,1)];
% Sulfur_coord = [Sulfur_coord(2:end,:) 16*ones(Sulfur_size(1)-1,1)];
% Nitrogen_coord = [Nitrogen_coord(2:end,:) 7*ones(Nitrogen_size(1)-1,1)];

Carbon_coord = Carbon_coord(2:end,:);
Oxygen_coord = Oxygen_coord(2:end,:);
Sulfur_coord = Sulfur_coord(2:end,:);
Nitrogen_coord = Nitrogen_coord(2:end,:);

Carbon_B_factor = Carbon_B_factor(2:end,:);
Oxygen_B_factor = Oxygen_B_factor(2:end,:);
Sulfur_B_factor = Sulfur_B_factor(2:end,:);
Nitrogen_B_factor = Nitrogen_B_factor(2:end,:);

% length of each array
Carbon_length = length(Carbon_B_factor);
Oxygen_length = length(Oxygen_B_factor);
Sulfur_length = length(Sulfur_B_factor);
Nitrogen_length = length(Nitrogen_B_factor);

% The standard deviation of each individual atom
Carbon_sigma = sqrt(Carbon_B_factor/(8*pi^2));
Oxygen_sigma = sqrt(Oxygen_B_factor/(8*pi^2));
Sulfur_sigma = sqrt(Sulfur_B_factor/(8*pi^2));
Nitrogen_sigma = sqrt(Nitrogen_B_factor/(8*pi^2));



% In this part of the code, 3D electron density volume array is calculated
% by interpolation.
% try interp3 in future.

% creat a triple cell with dimension 40*40*(3*45)
triple_cell = zeros(41,41,3*46);
trans_vector = [20 -2 5]; % translate the asymmetric unit to the center

% For the follow part, you can write a separate function file.
% It will look much more concise.

% calculate the density projection by Carbon
Carbon_coord_shift = Carbon_coord+ repmat(trans_vector,[Carbon_length,1]);
Carbon_e_density = zeros(41,41,46);
% sigma = 1; % You should change the code a little bit. Read the B factor
                % from pdb file and calculate sigma.
% B = 8*pi^2*sigma^2;
xyz_low = floor(Carbon_coord_shift);
xyz_frac = Carbon_coord_shift-xyz_low;
%xyz_high = xyz_low + 1;
for i=1:Carbon_length
    % distance to surrounding 8 points
    r1 = norm([xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r2 = norm([xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r3 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r4 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r5 = norm([1-xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r6 = norm([1-xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r7 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r8 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    
    sigma = Carbon_sigma(i);
    
    Carbon_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)) = Carbon_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3))+ 6*gaussmf(r1,[sigma,0]);
    Carbon_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) = Carbon_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) + 6*gaussmf(r2,[sigma,0]);
    Carbon_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) = Carbon_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) + 6*gaussmf(r3,[sigma,0]);
    Carbon_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1) = Carbon_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1)+ 6*gaussmf(r4,[sigma,0]);
    Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) = Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) + 6*gaussmf(r5,[sigma,0]);
    Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) = Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) + 6*gaussmf(r6,[sigma,0]);
    Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) = Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) + 6*gaussmf(r7,[sigma,0]);
    Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) =  Carbon_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) + 6*gaussmf(r8,[sigma,0]);
    
end
    
% calculate the density projection by Oxygen
Oxygen_coord_shift = Oxygen_coord+ repmat(trans_vector,[Oxygen_length,1]);
Oxygen_e_density = zeros(41,41,46);
sigma = 1; % You should change the code a little bit. Read the B factor
                % from pdb file and calculate sigma.
xyz_low = floor(Oxygen_coord_shift);
xyz_frac = Oxygen_coord_shift-xyz_low;
%xyz_high = xyz_low + 1;
for i=1:Oxygen_length
    % distance to surrounding 8 points
    r1 = norm([xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r2 = norm([xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r3 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r4 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r5 = norm([1-xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r6 = norm([1-xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r7 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r8 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    
    sigma = Oxygen_sigma(i);
    
    Oxygen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)) = Oxygen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3))+ 8*gaussmf(r1,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) = Oxygen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) + 8*gaussmf(r2,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) = Oxygen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) + 8*gaussmf(r3,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1) = Oxygen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1)+ 8*gaussmf(r4,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) = Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) + 8*gaussmf(r5,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) = Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) + 8*gaussmf(r6,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) = Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) + 8*gaussmf(r7,[sigma,0]);
    Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) =  Oxygen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) + 8*gaussmf(r8,[sigma,0]);
    
end

% calculate the density projection by Nitrogen
Nitrogen_coord_shift = Nitrogen_coord+ repmat(trans_vector,[Nitrogen_length,1]);
Nitrogen_e_density = zeros(41,41,46);
sigma = 1; % You should change the code a little bit. Read the B factor
                % from pdb file and calculate sigma.
xyz_low = floor(Nitrogen_coord_shift);
xyz_frac = Nitrogen_coord_shift-xyz_low;
% xyz_high = xyz_low + 1;
for i=1:Nitrogen_length
    % distance to surrounding 8 points
    r1 = norm([xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r2 = norm([xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r3 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r4 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r5 = norm([1-xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r6 = norm([1-xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r7 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r8 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    
    sigma = Nitrogen_sigma(i);
    
    Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)) = Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3))+ 7*gaussmf(r1,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) = Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) + 7*gaussmf(r2,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) = Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) + 7*gaussmf(r3,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1) = Nitrogen_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1)+ 7*gaussmf(r4,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) = Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) + 7*gaussmf(r5,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) = Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) + 7*gaussmf(r6,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) = Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) + 7*gaussmf(r7,[sigma,0]);
    Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) =  Nitrogen_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) + 7*gaussmf(r8,[sigma,0]);
    
end

% calculate the density projection by Sulfur
Sulfur_coord_shift = Sulfur_coord+ repmat(trans_vector,[Sulfur_length,1]);
Sulfur_e_density = zeros(41,41,46);
sigma = 1; % You should change the code a little bit. Read the B factor
                % from pdb file and calculate sigma.
xyz_low = floor(Sulfur_coord_shift);
xyz_frac = Sulfur_coord_shift-xyz_low;
% xyz_high = xyz_low + 1;
for i=1:Sulfur_length
    % distance to surrounding 8 points
    r1 = norm([xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r2 = norm([xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r3 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r4 = norm([xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r5 = norm([1-xyz_frac(i,1) xyz_frac(i,2) xyz_frac(i,3)]);
    r6 = norm([1-xyz_frac(i,1) xyz_frac(i,2) 1-xyz_frac(i,3)]);
    r7 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) xyz_frac(i,3)]);
    r8 = norm([1-xyz_frac(i,1) 1-xyz_frac(i,2) 1-xyz_frac(i,3)]);
    
    sigma = Sulfur_sigma(i);
    
    Sulfur_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)) = Sulfur_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3))+ 16*gaussmf(r1,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) = Sulfur_e_density(xyz_low(i,1), xyz_low(i,2),xyz_low(i,3)+1) + 16*gaussmf(r2,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) = Sulfur_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)) + 16*gaussmf(r3,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1) = Sulfur_e_density(xyz_low(i,1), xyz_low(i,2)+1,xyz_low(i,3)+1)+ 16*gaussmf(r4,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) = Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)) + 16*gaussmf(r5,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) = Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2),xyz_low(i,3)+1) + 16*gaussmf(r6,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) = Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)) + 16*gaussmf(r7,[sigma,0]);
    Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) =  Sulfur_e_density(xyz_low(i,1)+1, xyz_low(i,2)+1,xyz_low(i,3)+1) + 16*gaussmf(r8,[sigma,0]);
    
end

% Final 3d electron density map for one unit cell with one asymetric unit
unit_cell = Carbon_e_density+Nitrogen_e_density+Oxygen_e_density+Sulfur_e_density;
save('unit_cell.mat','unit_cell');
    
    























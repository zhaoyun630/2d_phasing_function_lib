function [spinelPotential,spinelHKL] = getSpinelCrystal(bFactor,maxHKL)

%*************************************************************************
% This function is used to get spinel crystal charge density distribution
% at given resolution(max hkl).
% Editted by Yun Zhao on May 5th, 2015
%*************************************************************************

% elements in crystal
latticeConstant = 8.080; %8~9 Angstrom
ZMg = 12; % Atom number
ZAl = 13;
ZO = 8;

% Atom coordinate in asymmetric unit?
% x=3/8;
% MgCoord = [0 0 0;1/4 1/4 1/4];
% AlCoord = [5/8 5/8 5/8;5/8 7/8 7/8;7/8 5/8 7/8;7/8 7/8 5/8];
% OCoord = [x x x;x -x -x;-x x -x;-x -x x;1/4-x 1/4-x 1/4-x;1/4-x 1/4+x 1/4+x;1/4+x 1/4-x 1/4+x;1/4+x 1/4+x 1/4-x];

% Non-equivalent Atom coordinate in one unit cell. 
% Here divide the cell to 8 layers [from 0 to 7th].
MgCoord0 = [0 0 0;2 2 0]/4;
MgCoord2 = [1 1 1;3 3 1]/4;
MgCoord4 = [0 2 2;2 0 2]/4;
MgCoord6 = [1 3 3;3 1 3]/4;
MgCoord =[MgCoord0;MgCoord2;MgCoord4;MgCoord6];

AlCoord1 = [1 5 1;3 7 1;5 1 1;7 3 1]/8;
AlCoord3 = [1 7 3;3 5 3;5 3 3;7 1 1]/8;
AlCoord5 = [1 1 1;3 3 5;5 5 5;7 7 5]/8;
AlCoord7 = [1 3 7;3 1 7;5 7 7;7 5 7]/8;
AlCoord = [AlCoord1;AlCoord3;AlCoord5;AlCoord7];

OCoord1 = [1 3 1;1 7 1;3 1 1;3 5 1;5 3 1;5 7 1;7 1 1;7 5 1]/8;
OCoord3 = [1 1 3;1 5 3;3 3 3;3 7 3;5 1 3;5 5 3;7 3 3;7 7 3]/8;
OCoord5 = [1 3 5;1 7 5;3 1 5;3 5 5;5 3 5;5 7 5;7 1 5;7 5 5]/8;
OCoord7 = [1 1 7;1 5 7;3 3 7;3 7 7;5 1 7;5 5 7;7 3 7;7 7 7]/8;
OCoord = [OCoord1;OCoord3;OCoord5;OCoord7];



% Get Gaussian coefficients and make A and B matrices for atomic scattering
% factors.
MgGCoeff = element(ZMg);
AlGCoeff = element(ZAl);
OGCoeff = element(ZO);

% Note each Coefficient has 6 parameters.
MgA = MgGCoeff(1,:);
MgB = MgGCoeff(2,:);
AlA = AlGCoeff(1,:);
AlB = AlGCoeff(2,:);
CA = OGCoeff(1,:);
CB = OGCoeff(2,:);

% Generate matrix with hkl values. Generates an (2N+1)x(2N+1)x(2N+1)
% matrix of Miller indices that will correspond to 3D matrix of Fourier
% coefficients, V(g). At same time, generate 3D structure factor
% matrix, Fhkl.
    h=-maxHKL:maxHKL;
    k=-maxHKL:maxHKL;
    l=-maxHKL:maxHKL;
    Fhkl=zeros(length(h),length(k),length(l));
%     Vg=zeros(length(h),length(k),length(l));
%    sa=zeros(1,numel(Vg));
%    fs=zeros(1,numel(Vg));
	numel=length(h)*length(k)*length(l);
	sa=zeros(1,numel);
	fsMg=zeros(1,numel);
    fsAl=zeros(1,numel);
    fsO=zeros(1,numel);
    count=1;
    for m=1:length(h)
        for n=1:length(k)
            for p=1:length(l)
                % First HKL matrices as vectors
                g=[h(m) k(n) l(p)];
                %hkl(m,n,p)={g};
                % Now Vg matrix
                SS=(g*g')/latticeConstant^2;
                S=0.5*sqrt(SS); % S=1/2*d
                sa(count)=S;
         
                % Doyle-Turner SF for electrons is Sum(i=1,4)A(i)*
                % exp(-B(i)*S^2)
                fMg = 0;
                fAl = 0;
                fO = 0;
                for q=1:4
                    fMg=fMg+MgA(q)*exp(-MgB(q)*S^2);
                    fAl=fAl+AlA(q)*exp(-AlB(q)*S^2);
                    fO=fO+CA(q)*exp(-CB(q)*S^2);
                end
                fsMg(count)=fMg*exp(-bFactor*S^2); % DW factor (or called B factor)
                fsAl(count)=fAl*exp(-bFactor*S^2);
                fsO(count)=fO*exp(-bFactor*S^2);
                % Now form structure factor.
             
                Fhkl(m,n,p)=0;
                for nMg = 1:length(MgCoord)
                    Fhkl(m,n,p)=Fhkl(m,n,p)+fsMg(count)*exp(2*pi*1i*dot(g,MgCoord(nMg,:)));
                end
                
                for nAl=1:length(AlCoord)
                    Fhkl(m,n,p)=Fhkl(m,n,p)+fsAl(count)*exp(2*pi*1i*dot(g,AlCoord(nAl,:)));
                end
                
                for nO=1:length(OCoord)
                    Fhkl(m,n,p)=Fhkl(m,n,p)+fsO(count)*exp(2*pi*1i*dot(g,OCoord(nO,:)));
                end
                
                count=count+1;
            end
        end
    end
% Now form Vg matrix: Fourier coefficients in volts.
    spinelHKL=(47.87*Fhkl)/(latticeConstant^3);
    spinelPotential = real(ifftn(ifftshift(spinelHKL)));
    plot(fsMg,sa);
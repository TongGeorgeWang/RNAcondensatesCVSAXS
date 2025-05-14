
clear
 
%% Create an initial state for a LAMMPS MD simulation 
%   see https://docs.lammps.org/2001/data_format.html for formatting doc
%
%   This is specifically for building polymers (molecules of bonded atoms) 
%       For example: rU30, rA30, PR30
%       Here, Atoms and Bonds are most important to define.
%       'Molecule' and 'polymer' are used interchangeably to define a
%       bonded array of 'atoms'.
%
%   This particular script defines two types of molecule in the box 
%       (eg: ssRNA strand interspersed with PR30 strand) 
%   Both types of molecules will have the same # of molecules (ie same
%   concentration, equimolar); differential stoichiometry may be programmed
%   in later.
%
%   MPiPi force field atom types, all atoms relevant to PR/K-RNA simulations: 
%       3:K|5:R|18:P||41:rA|42:rC|43:rG|44:rU | see doc for rest of atom types
%   Charges: 
%       K & R have charge +0.75; NAs have charge -0.75
%   bond types: 
%       1: NA bond (5.00A) | 2: AA bond (3.81A)
%
%   GW April 2024
%
%   August 2024 - added foldovers in molecule 2 to allow for longer peptides 
%


%% Define 

saveName = 'rU30PR30_250uM_equimolar';
atomMasses = [131.2 57.05 128.2 101.1 156.2 71.08 115.1 129.1 163.2 99.06999999999999 113.2 128.1 186.2 ...
    147.2 87.08 137.1 114.1 97.12 103.1 113.2 131.2 57.05 128.2 101.1 156.2 71.08 115.1 129.1 163.2 99.06999999999999 ...
    113.2 128.1 186.2 147.2 87.08 137.1 114.1 97.12 103.1 113.2 329.2 305.2 345.2 306.2]; % Masses of each unique atom type (in order)
molLen = [30 61]; % Number of atoms in each molecule
cM = 250E-6; % molecule concentration, in mol/L 

moleculeAtomTypeArray1 = ones(1,molLen(1))*44; % molecule 1 sequence (indexed in order that atomMasses is defined) 
nMonomersinFold1 = numel(moleculeAtomTypeArray1); % how many polymers to include in y-series before folding to avoid overlap in space 
moleculeAtomChargeArray1 = ones(1,molLen(1))*-0.75; % molecule 1 charge sequence
moleculeBondTypeArray1 = ones(1,molLen(1)-1);


moleculeAtomTypeArray2 = [5,repmat([18 5],1,30)]; % molecule 2 sequence (indexed in order that atomMasses is defined) 
nMonomersinFold2 = 40; % how many polymers to include in y-series before folding to avoid overlap in space 
moleculeAtomChargeArray2 =  [0.75,repmat([0 0.75],1,30)]; % molecule 2 charge sequence
moleculeBondTypeArray2 = ones(1,molLen(2)-1)*2;

bondLengths = [5.00 3.81]; % Lengths of bonds: NA=5.0, AA=3.81

% Distances from edge atom to edge of box 
x0 = 70; 
y0 = 60; 
z0 = 80; 

Nstacks = 5; % (in z) How many molecules you want to stack before moving to next stack
Narrays = 5; % (in x) How many stack arrays 
NpolymersInARow = 3; % (in y) How many polymers you want to array, ass-to-mouth



%% Compute the rest of the box parameters 

[dx,dy,dz,Lz,Lx,Ly] = setUpBox(cM, numel(molLen),z0,x0,y0,Nstacks,Narrays,NpolymersInARow,[nMonomersinFold1 nMonomersinFold2],bondLengths);
nMolecules = [Nstacks*Narrays*NpolymersInARow Nstacks*Narrays*NpolymersInARow]; % Number of each molecule 


%% Define a polymer and copy-translate it according to the input dimensions

nAtoms1 = molLen(1).*nMolecules(1);
moleculeIDs1 = repelem((1:1:nMolecules(1)),molLen(1));

AtomLines1 = zeros([molLen(1)*nMolecules(1) 7]);
% Write atom info

%atomIDLines = (1:1:molLen(1)) + (molLen(1).*(molNum-1)); % equivalently, atomIDLines=(1:1:nAtoms)

AtomLines1(:,1) = (1:1:nAtoms1); % atom identifier
AtomLines1(:,2) = moleculeIDs1'; % molecule identifier
AtomLines1(:,3) = repmat(moleculeAtomTypeArray1,1,nMolecules(1)); % atom type
AtomLines1(:,4) = repmat(moleculeAtomChargeArray1,1,nMolecules(1)); % atom charge


%% Lay out coordinates, arraying polymers in [Nstacks x Narrays x NPolymersInARow]

% Build unit polymer 
xPolymer(1) = x0; yPolymer(1) = y0(1); zPolymer(1) = z0; % declare position of first atom 
atomIDLines = (1:1:molLen(1));
for N = 1:molLen(1)-1
    xPolymer(N+1) = xPolymer(N);
    yPolymer(N+1) = yPolymer(N) + bondLengths(1);
    zPolymer(N+1) = zPolymer(N);
end

% Translate unit polymer in Cartesian space, separately for each coordinate
% z (stacks occur along z)
zSeries = repelem((zPolymer(1):dz:zPolymer(1)+dz*(Nstacks-1)),molLen(1));
zPolymers = repmat(zSeries,1,(Narrays*NpolymersInARow));
AtomLines1(:,7) = zPolymers'; % write atom z-coordinate

% x (parallel stack arrays occur along x)
xSeries = repelem((xPolymer(1):dx:dx*(Narrays)), molLen(1)*Nstacks);
xPolymers = repmat(xSeries,1,(NpolymersInARow));
AtomLines1(:,5) = xPolymers'; % write atom x-coordinate

% y (the polymer is extended along y)
polymerLen = (molLen(1)-1)*bondLengths(1);
yUnit(1,:) = yPolymer;
for i = 2:NpolymersInARow
    yUnit(i,:) = yUnit(i-1,:)+dy(1)+polymerLen;
end
yUnit = yUnit'; ySeries = yUnit(:); ySeries = ySeries';
yPolymers = repmat(ySeries,1,(Narrays*Nstacks));
AtomLines1(:,6) = yPolymers'; % write atom y-coordinate


%% Write bond matrix

nBondsPerMol = (molLen(1)-1);
nBonds1 = nMolecules(1)*nBondsPerMol;
BondLines1 = zeros([nBonds1 4]);
BondLines1(:,1) = (1:1:nBonds1); % write bond ID 
BondLines1(:,2) = repmat(moleculeBondTypeArray1,1,nMolecules(1));

for j = 1:nMolecules(1)
    BondLines1(1+(j-1)*nBondsPerMol:(j)*nBondsPerMol,3) = (1+(j-1)*molLen(1):nBondsPerMol+(j-1)*molLen(1));
    BondLines1(1+(j-1)*nBondsPerMol:(j)*nBondsPerMol,4) = (2+(j-1)*molLen(1):1+nBondsPerMol+(j-1)*molLen(1));
end


%% Second polymer: 
%% Define a polymer and copy-translate it according to the input dimensions

nAtoms2 = molLen(2).*nMolecules(2);
moleculeIDs = repelem((nMolecules(1)+1:1:nMolecules(2)+nMolecules(1)),molLen(2)); % start where first molecule IDs left off

AtomLines2 = zeros([molLen(2)*nMolecules(2) 7]);
% Write atom info
AtomLines2(:,1) = (nAtoms1+1:1:nAtoms1+nAtoms2); % atom identifier
AtomLines2(:,2) = moleculeIDs'; % molecule identifier
AtomLines2(:,3) = repmat(moleculeAtomTypeArray2,1,nMolecules(2)); % atom type
AtomLines2(:,4) = repmat(moleculeAtomChargeArray2,1,nMolecules(2)); % atom charge


%% Lay out coordinates, arraying polymers in [Nstacks x Narrays x NPolymersInARow]

% Build unit polymer 
clear xPolymer yPolymer zPolymer xPolymers yPolymers zPolymers yUnit

% Predefine array sizes 
xPolymer = zeros([1 molLen(2)]); yPolymer = zeros([1 molLen(2)]); zPolymer = zeros([1 molLen(2)]);

xPolymer(1) = x0; yPolymer(1) = y0(1); zPolymer(1) = z0; % declare position of first atom 
atomIDLines = (1:1:molLen(2));

nFoldovers = ceil(molLen(2)./nMonomersinFold2); % if too long, will need to fold polymer or else polymers will overlap in space

for fold = 1:nFoldovers
    if fold == nFoldovers % if on the last fold, may need to terminate polymer early if length is not an integer multiple of 60
        nPolymersInLastFold = molLen(2)-(fold-1).*nMonomersinFold2 - 1 ;
        for N = nMonomersinFold2*(fold-1) + (1:nPolymersInLastFold)
            if rem(fold, 2) == 1 % odd fold #
                xPolymer(N+1) = xPolymer(N);
                yPolymer(N+1) = yPolymer(N) + bondLengths(2); % advance polymer in y
                zPolymer(N+1) = zPolymer(N);

            else % odd fold #
                xPolymer(N+1) = xPolymer(N);
                yPolymer(N+1) = yPolymer(N) - bondLengths(2); % advance polymer in y in opposite direction (ie, fold)
                zPolymer(N+1) = zPolymer(N);

            end
        end

    else % not on last fold 


        if rem(fold, 2) == 1 % odd fold #
            for N = nMonomersinFold2*(fold-1) + (1:nMonomersinFold2)
                xPolymer(N+1) = xPolymer(N);
                yPolymer(N+1) = yPolymer(N) + bondLengths(2); % advance polymer in y
                zPolymer(N+1) = zPolymer(N);
            end
            yPolymer(N+1) = yPolymer(N); % for last polymer in row before folding, do not add a bond length to avoid skewing
        else % even fold #
            for N = nMonomersinFold2*(fold-1) + (1:nMonomersinFold2)
                xPolymer(N+1) = xPolymer(N);
                yPolymer(N+1) = yPolymer(N) - bondLengths(2); % advance polymer in y in opposite direction (ie, fold)
                zPolymer(N+1) = zPolymer(N);
            end
            yPolymer(N+1) = yPolymer(N); % for last polymer in row before folding, do not add a bond length to avoid skewing
        end
        zPolymer(N+1) = zPolymer(N) + bondLengths(2); % move to next z row (fold)
    end
end


% Translate unit polymer in Cartesian space, separately for each coordinate
% z (stacks occur along z)
%zSeries = repelem((zPolymer(1):dz:zPolymer(1)+dz*(Nstacks-1)),molLen(2));
%zPolymers = repmat(zSeries,1,(Narrays*NpolymersInARow));

zSeries = zPolymer;
for fold = 1:Nstacks-1
    zSeries = [zSeries, (zPolymer+dz.*(fold))];
end
zPolymers = repmat(zSeries,1,(Narrays*NpolymersInARow));
AtomLines2(:,7) = zPolymers'; % write atom z-coordinate

% x (parallel stack arrays occur along x)
xSeries = repelem((xPolymer(1)+(dx/2):dx:dx*(Narrays)+(dx/2)), molLen(2)*Nstacks);
xPolymers = repmat(xSeries,1,(NpolymersInARow));
AtomLines2(:,5) = xPolymers'; % write atom x-coordinate

% y (the polymer is extended along y)
polymerLen = nMonomersinFold2*bondLengths(2);
%polymerLen = (molLen(2)-1)*bondLengths(2);
yUnit(1,:) = yPolymer;
for i = 2:NpolymersInARow
    yUnit(i,:) = yUnit(i-1,:)+dy(2)+polymerLen;
end
yUnit = yUnit'; ySeries = yUnit(:); ySeries = ySeries';
yPolymers = repmat(ySeries,1,(Narrays*Nstacks));
AtomLines2(:,6) = yPolymers'; % write atom y-coordinate


%% Write bond matrix

nBondsPerMol = (molLen(2)-1);
nBonds2 = nMolecules(2)*nBondsPerMol;
BondLines2 = zeros([nBonds2 4]);
BondLines2(:,1) = (nBonds1+1:1:nBonds1+nBonds2); % write bond ID 
BondLines2(:,2) = repmat(moleculeBondTypeArray2,1,nMolecules(2));

for j = 1:nMolecules(2)
    BondLines2(1+(j-1)*nBondsPerMol:(j)*nBondsPerMol,3) = (nAtoms1+1+(j-1)*molLen(2):nAtoms1+nBondsPerMol+(j-1)*molLen(2));
    BondLines2(1+(j-1)*nBondsPerMol:(j)*nBondsPerMol,4) = (nAtoms1+2+(j-1)*molLen(2):nAtoms1+1+nBondsPerMol+(j-1)*molLen(2));
end



%% Define box size and dimensions based on how molecules are packed

%boxDimX = [0, max(xPolymers) + 2.*x0];
%boxDimY = [0, max(yPolymers) + 2.*y0];
%boxDimZ = [0, max(zPolymers) + 2.*z0];

boxDimX = [0,Lx];
boxDimY = [0,Ly];
boxDimZ = [0,Lz];


%% Define header 
nAtomTypes = numel(atomMasses);
nBondTypes = numel(bondLengths); 

hLine = 1; % to keep track of which line of the header is being written to
header{hLine,1} = 'LAMMPS data file written via buildLAMMPSconfig.m'; hLine = hLine+1; % Title 
header{hLine,1} = ' '; hLine = hLine+1; % necessary whitespace
header{hLine,1} = [num2str(nAtoms1+nAtoms2),' atoms']; hLine = hLine+1; % Total number of atoms
header{hLine,1} = [num2str(nBonds1+nBonds2),' bonds']; hLine = hLine+1; % Total number of bonds
header{hLine,1} = ' '; hLine = hLine+1; % necessary whitespace
header{hLine,1} = [num2str(nAtomTypes),' atom types']; hLine = hLine+1; % Number of unique atoms
header{hLine,1} = [num2str(nBondTypes),' bond types']; hLine = hLine+1; % Number of unique bonds 

header{hLine,1} = ' '; hLine = hLine+1; % necessary whitespace
header{hLine,1} = [num2str(boxDimX(1)),' ',num2str(boxDimX(2)),' xlo xhi']; hLine = hLine+1; % box x dimension
header{hLine,1} = [num2str(boxDimY(1)),' ',num2str(boxDimY(2)),' ylo yhi']; hLine = hLine+1; % box y dimension
header{hLine,1} = [num2str(boxDimZ(1)),' ',num2str(boxDimZ(2)),' zlo zhi']; hLine = hLine+1; % box z dimension

header{hLine,1} = ' '; hLine = hLine+1; % necessary whitespace
header{hLine,1} = 'Masses'; hLine = hLine+1; 
header{hLine,1} = ' '; hLine = hLine+1; % necessary whitespace

for a = 1:nAtomTypes % Define atom masses for each atom type
    header{hLine,1} = [num2str(a),' ',num2str(atomMasses(a))]; hLine = hLine+1;
end




%% Write to .dat (text) file) 

AtomLines1 = arrayfun(@num2str, AtomLines1, 'UniformOutput', 0);
AtomLines1cat = strcat(AtomLines1(:,1),{' '},AtomLines1(:,2),{' '},AtomLines1(:,3),{' '},AtomLines1(:,4),{' '},...
                    AtomLines1(:,5),{' '},AtomLines1(:,6),{' '},AtomLines1(:,7),{' '});
BondLines1 = arrayfun(@num2str, BondLines1, 'UniformOutput', 0);
BondLines1cat = strcat(BondLines1(:,1),{' '},BondLines1(:,2),{' '},BondLines1(:,3),{' '},BondLines1(:,4));


AtomLines2 = arrayfun(@num2str, AtomLines2, 'UniformOutput', 0);
AtomLines2cat = strcat(AtomLines2(:,1),{' '},AtomLines2(:,2),{' '},AtomLines2(:,3),{' '},AtomLines2(:,4),{' '},...
                    AtomLines2(:,5),{' '},AtomLines2(:,6),{' '},AtomLines2(:,7),{' '});
BondLines2 = arrayfun(@num2str, BondLines2, 'UniformOutput', 0);
BondLines2cat = strcat(BondLines2(:,1),{' '},BondLines2(:,2),{' '},BondLines2(:,3),{' '},BondLines2(:,4));
%BondLines2cat = BondLines2cat(nBonds1+1:end); % remove emptyness


fullEssay = [header; ' '; 'Atoms'; ' '; AtomLines1cat; AtomLines2cat; ' '; 'Bonds'; ' '; BondLines1cat; BondLines2cat];
writecell(fullEssay, [saveName,'.dat'], 'FileType','text','Delimiter',' ','QuoteStrings',0)




%% functions 
function [dx,dy,dz,Lz,Lx,Ly] = setUpBox(cM, NmoleculeTypes,z0,x0,y0,Ns,Na,Np,lp,lb)
% Determine box parameters (see Computational Notebook A page 140-141 for full formalism)

c = cM*6.023E23/0.001/1E30; % convert to molecules/A^3
N = Ns*Na*Np; % total # of molecules

dz = ((N/c)^(1/3) - 2*z0)/(Ns-1) ;
Lz = 2*z0 + (Ns-1)*dz;

if NmoleculeTypes == 1
    dx = ((N/c)^(1/3) - 2*x0)/(Na-1);
    Lx = 2*x0 + (Na-1)*dx;
    dy = ((N/c)^(1/3) - 2*y0 - (lp-1)*lb*Np)/(Np-1); 
    if dy <= 0 
        error('Not enough spacing between polymers (Np in y-direction is too large)')
    end
    Ly = 2*y0 + (Np-1)*dy+(lp-1)*lb*Np;

elseif NmoleculeTypes == 2
    dx = 2*((N/c)^(1/3) - 2*x0) / (2*Na-1);
    Lx = 2*x0 + (Na-1)*dx + dx/2;
    dy(1) = ((N/c)^(1/3) - 2*y0 - (lp(1)-1)*lb(1)*Np)/(Np-1); 
    dy(2) = ((N/c)^(1/3) - 2*y0 - (lp(2)-1)*lb(2)*Np)/(Np-1); 
    % This error is suppressed since polymer foldovers are incorporated now 
    % if dy(1) <= 0 || dy(2) <= 0 
    %     error('Not enough spacing between polymers (Np in y-direction is too large)')
    % end
    Ly = 2*y0 + (Np-1)*dy(1)+(lp(1)-1)*lb(1)*Np; % Ly should be the same for both 1 and 2 
end


V = Lz*Lx*Ly;
actualC = N/V; % check actual concentration in molecules/A^3
actualCm = actualC*0.001*1E30/6.023E23*1E6; % convert to uM
end


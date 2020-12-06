%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% layers structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the conduction band offset in eV
% second column is the length of the layer in nm
% third column is the n doping volumique of that layer in 1e18cm-3 

% You have to put a resonable amount of doping! Otherwise, it will diverge 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% superlattice

% GaAs=0;
% AlGaAs40=0.36;
% meff   = 0.067;
% Epsi   = 10;
% 
% M=[
% GaAs       20   2    % contact on the left
% GaAs       150  0    % spacer
% 
% AlGaAs40   10  0.5
% GaAs       10  0
% AlGaAs40   10  0.5
% GaAs       10  0
% AlGaAs40   10  0.5
% GaAs       10  0
% AlGaAs40   10  0.5
% GaAs       10  0
% AlGaAs40   10  0.5
% GaAs       10  0
% AlGaAs40   10  0.5
% 
% GaAs       150  0    % spacer
% GaAs       20   1    % contact on the right
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% superlattice with Schottky contact

% GaAs=0;
% AlGaAs40=0.36;
% meff   = 0.067;
% Epsi   = 10;
% 
% M=[
% GaAs       20   0    % contact on the left
% GaAs       150  0    % spacer
% 
% AlGaAs40   10  0.1
% GaAs       10  0
% AlGaAs40   10  0.1
% GaAs       10  0
% AlGaAs40   10  0.1
% GaAs       10  0
% AlGaAs40   10  0.1
% GaAs       10  0
% AlGaAs40   10  0.1
% GaAs       10  0
% AlGaAs40   10  0.1
% 
% GaAs       150  0    % spacer
% GaAs       20   1    % contact on the right
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coupled Quantum wells

% InGaAs=0;
% AlInAs=0.52;
% AlInAs2=0.2;
% meff   = 0.042;
% Epsi   = 10;
% 
% M=[
% InGaAs     20   1    % contact on the left
% InGaAs     50   0    % spacer
% 
% AlInAs      5   0
% AlInAs      1   5
% AlInAs      5   0
% InGaAs      8   0
% AlInAs      3   0
% InGaAs     20   0
% AlInAs      5   0
% AlInAs      1   5
% AlInAs      5   0
% 
% InGaAs     50   0    % spacer
% InGaAs     20   1    % contact on the right
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InGaAs=-0.2;
GaAs=0;
meff   = 0.05;
Epsi   = 10;

M=[
GaAs    50   1    % contact on the left
GaAs    70   0    % spacer

GaAs    10   0
GaAs     1   1
GaAs     5   0
InGaAs   20  0
GaAs     5   0
GaAs     1   1
GaAs    10   0

GaAs    70   0    % spacer
GaAs    50   1    % contact on the right
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GaN  = 0;
% AlN  = 1.8;
% meff = 0.22;
% Epsi = 10;
% DF   = 10;        % Electrical field discontinuity [MV/cm]
% 
% Lb = 1;           % barrier thickness [nm]
% Lw = 1;           % well thickness [nm]
% Ls = 0.05;        % doping spike thickness [nm] (in order to get the E-field)
% 
% Fb = +DF*(Lw+2*Ls)/(Lw+Lb+4*Ls)*1e6*1e2;   %[V/m]
% Fw = -DF*(Lb+2*Ls)/(Lw+Lb+4*Ls)*1e6*1e2;   %[V/m]
% 
% dopS = DF*1e6*1e2*Epsi*Epsi0/e;   % charge/m2 MUST BE added on the interface for GaN/AlN for Wurtzite
% dopV = dopS/(Ls*1e-9);            % charge/m3
% dopV = dopV*1e-6*1e-18;           % charge 1e18cm-3
% 
% 
% M=[
%     
% GaN     10   50    % contact on the left
% GaN     10   0     % spacer
% 
% GaN     Ls  +dopV
% AlN     Lb   0
% GaN     Ls  -dopV
% GaN     Lw   5
% GaN     Ls  +dopV
% AlN     Lb   0
% GaN     Ls  -dopV
% GaN     Lw   5
% GaN     Ls  +dopV
% AlN     Lb   0
% GaN     Ls  -dopV
% GaN     Lw   5
% GaN     Ls  +dopV
% AlN     Lb   0
% GaN     Ls  -dopV
% GaN     Lw   5
% GaN     Ls  +dopV
% AlN     Lb   0
% GaN     Ls  -dopV
% GaN     Lw   5
% GaN     Ls  +dopV
% AlN     Lb   0
% GaN     Ls  -dopV
% 
% GaN     10   0     % spacer
% GaN     10   50    % contact on the right
% ];

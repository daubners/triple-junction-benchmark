%% Define physical input parameters
%  together with refernece values for non-dimensionalization

% Physical length in 3 dimensions
L = 100e-6;  % [m]
B = 200e-6; % [m]

% Interfacial energies
gamma_ab = 1.0; % [J/m2]
gamma_ac = 1.0;
gamma_bc = [0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,1.9,2.0];

% Mobilities
M_ab = 1e-14;   % [m^4/(Js)]

%% Non-dimensional form
% Reference values
t_ref = 100;  % [s]
x_ref = 1e-6; % [m]
O_ref = 1e6;  % [J/m3]

% Non-dimensional length
Lx = L/x_ref;
Ly = B/x_ref;

% Non-dimensional surface energies
gab = gamma_ab/(O_ref*x_ref);
gac = gamma_ac/(O_ref*x_ref);
gbc = gamma_bc/(O_ref*x_ref);
mab = M_ab*O_ref*t_ref/x_ref;
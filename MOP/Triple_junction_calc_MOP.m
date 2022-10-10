% Explicit multi-order parameter solver based on finite difference scheme
%
% Simon Daubner simon.daubner@kit.edu
% Department of Mechanical Engineering
% Karlsruhe Institute of Technology

clear;
close all;

%% Input for Discretisation
addpath('../')
input_parameters_triple;

% Grid is discretised with one row of boundary cells
nx = 102;          % Number of steps in space(x)
ny = 202;          % Number of steps in space(y)
dt = 0.1;          % Width of each time step
t_end = 10000.0;   % Total time (non-dimensional)
frames = 10;       % Plot and output frames

% Creates output as ascii .dat files and .vtk
writeOutput = 0;   % Write output true or false
folder = ['..\Simulations\Benchmark-2\MOP-grad0-wellmoelans-L0\'];
grid   = ['100x400'];
% welltoth wellgarcke wellmoelans

assignmentStrategy = 'L_const'; % 'k_const' 'L_const' 'Om_const'

compare = 0;           % Compare with other simulation true or false
comparepath = '..\Simulations\Benchmark-1\MOP-grad0-wellmoelans-k0\';

plotsAndvalidation = 1;

%% Prepare simulation parameters
ntot=cast(t_end/dt,'int32'); % Total number of time steps
nt_out=ntot/frames;

% Using one layer of ghost cells for boundary conditions
dx = Lx/(nx-2); % Width of space step(x)
dy = Ly/(ny-2); % Width of space step(y)

% Position of the cell centers
x=-dx/2:dx:Lx+dx/2;
y=-dy/2:dy:Ly+dy/2;

eps = 5*dx;
K = 9;

switch assignmentStrategy
    case 'k_const'
        params.k = eps*gab;
        params.Oab = K*gab/eps;
        params.Oac = K*gac/eps;
    case 'L_const'
        params.kab = eps*gab;
        params.kac = eps*gac;
        params.Oab = K*gab/eps;
        params.Oac = K*gac/eps;
    case 'Om_const'
        params.kab = eps*gab;
        params.kac = eps*gac;
        params.O = K*gab/eps;
end

params.nx = nx;
params.ny = ny;
params.dx = dx;
params.dy = dy;
params.Lx = Lx;
params.Ly = Ly;

% Mobility for dual interaction evolution equation
M0  = mab/eps;

for gamma=gbc

params.gratio = gamma/gab;
switch assignmentStrategy
    case 'k_const'
        params.Obc = K*gamma^2/(gab*eps);
    case 'L_const'
        params.kbc = eps*gamma;
        params.Obc = K*gamma/eps;
    case 'Om_const'
        params.kbc = eps*gamma^2/gab;
end

path = [folder,grid,'_gamma_ratio',num2str(gamma/gab)];

%% Initial conditions
phia = zeros(nx,ny); %Initializing phase alpha
phib = zeros(nx,ny); %Initializing phase beta

% Initilize fields for laplacians
laplace_a = zeros(nx,ny);
laplace_b = zeros(nx,ny);
laplace_c = zeros(nx,ny);

% Set initial and boundary conditions
phia(:,ny-20:ny) = ones(nx,21);
phib(1:nx/2,1:ny-21) = ones(nx/2,ny-21);
phic = ones(nx,ny) - phia - phib;

% Initialize figures
if plotsAndvalidation
    fig_phia  = figure;
    fig_ghost = figure;
end

if compare
    fig_diff  = figure;
end

if writeOutput
    save([path,'_phia_F0.dat'],'phia','-ascii');
    save([path,'_phib_F0.dat'],'phib','-ascii');
    save([path,'_phic_F0.dat'],'phic','-ascii');
    write_vtk([path,'_phia_F0.vtk'],'phia',phia,dx,dy,'ascii');
    write_vtk([path,'_phib_F0.vtk'],'phib',phib,dx,dy,'ascii');
    write_vtk([path,'_phic_F0.vtk'],'phic',phic,dx,dy,'ascii');
end
count=1; % Frame number now starts at 1;

%% Computation
i=2:nx-1;
j=2:ny-1;

tic
for it=1:ntot
    %% Multi-order parameter evolution compute all cells
    % [dphia, dphib, dphic] = calcPhi_k0_toth(phia,phib,phic,params);
    % [dphia, dphib, dphic] = calcPhi_L0_toth(phia,phib,phic,params);
    % [dphia, dphib, dphic] = calcPhi_Om0_toth(phia,phib,phic,params);
    % [dphia, dphib, dphic] = calcPhi_k0_moelans(phia,phib,phic,params);
    [dphia, dphib, dphic] = calcPhi_L0_moelans(phia,phib,phic,params);
    % [dphia, dphib, dphic] = calcPhi_Om0_moelans(phia,phib,phic,params);
    phia(i,j) = phia(i,j)-dt.*M0.*dphia;
    phib(i,j) = phib(i,j)-dt.*M0.*dphib;
    phic(i,j) = phic(i,j)-dt.*M0.*dphic;

    %% Boundary conditions
    % No-flux BC at left and right sides
    phia(1,:)  = phia(2,:);
    phia(nx,:) = phia(nx-1,:);
    phib(1,:)  = phib(2,:);
    phib(nx,:) = phib(nx-1,:);
    phic(1,:)  = phic(2,:);
    phic(nx,:) = phic(nx-1,:);
    
    % Periodic BC at left and right sides
    % phia_new(1,:)  = phia_new(nx-1,:);
    % phia_new(nx,:) = phia_new(2,:);
    % phib_new(1,:)  = phib_new(nx-1,:);
    % phib_new(nx,:) = phib_new(2,:);
    
    % No-flux BC at bottom and top
    phia(:,1)  = phia(:,2);
    phia(:,ny) = phia(:,ny-1);
    phib(:,1)  = phib(:,2);
    phib(:,ny) = phib(:,ny-1);
    phic(:,1)  = phic(:,2);
    phic(:,ny) = phic(:,ny-1);
    
    %% Create output every nt_out timestep
    if mod(it,nt_out)==0
        if plotsAndvalidation
            plot_field(phia.*phib+phia.*phic+phib.*phic,fig_phia,cast(it,'double')*dt,dx,dy);
            plot_xy(x,phia(:,2),fig_ghost,cast(it,'double')*dt,0,1);
        end
        % Write output
        if writeOutput
            save([path,'_phia_F',num2str(count),'.dat'],'phia','-ascii');
            save([path,'_phib_F',num2str(count),'.dat'],'phib','-ascii');
            save([path,'_phic_F',num2str(count),'.dat'],'phic','-ascii');
            write_vtk([path,'_phia_F',num2str(count),'.vtk'],'phia',phia,dx,dy,'ascii');
            write_vtk([path,'_phib_F',num2str(count),'.vtk'],'phib',phib,dx,dy,'ascii');
            write_vtk([path,'_phic_F',num2str(count),'.vtk'],'phic',phic,dx,dy,'ascii');
        end
        
        % Comparison to validation solution
        if compare
            phia_ref = load([comparepath,grid,'_gamma_ratio',num2str(gamma/gab),'_phia_F',num2str(count),'.dat'],'-ascii');
            phi_diff = phia - phia_ref;
            plot_field(phi_diff,fig_diff,cast(it,'double')*dt,dx,dy);
        end
        
        count = count+1;
    end
end
toc

end

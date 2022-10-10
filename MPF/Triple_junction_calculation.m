% Explicit multi-phase solver based on finite difference scheme
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
folder = ['..\Simulations\Benchmark-1\MPF-grad2-wellmoelans-chi15\'];
grid   = ['100x200'];

compare = 0;           % Compare with other simulation true or false
comparepath = '..\Simulations\00-First-shot\100x100_grad1_well_toth_assign1';

plotsAndvalidation = 1;

gradient  = 'grad1';
potential = 'obstacle'; % obstacle welltoth wellgarcke wellmoelans

%% Prepare simulation parameters
ntot=cast(t_end/dt,'int32'); % Total number of time steps
nt_out=ntot/frames;

% Using one layer of ghost cells for boundary conditions
dx = Lx/(nx-2); % Width of space step(x)
dy = Ly/(ny-2); % Width of space step(y)

% Position of the cell centers
x=-dx/2:dx:Lx+dx/2;
y=-dy/2:dy:Ly+dy/2;

eps = 5*dx; %4.8634*dx;

switch potential
    case 'welltoth'
        K = 9;
    case 'wellgarcke'
        K = 9;
    case 'wellmoelans'
        K = 9;
    case 'obstacle'
        K = 16/pi^2;
    otherwise
        error('Choose welltoth, wellmoelans, wellgarcke or obstacle for potential term!')
end

params.kab = eps*gab;
params.kac = eps*gac;
params.Oab = K*gab/eps;
params.Oac = K*gac/eps;

switch potential
    case 'wellgarcke'
    params.Oabc = 100.0*params.Oab;
end

params.nx = nx;
params.ny = ny;
params.dx = dx;
params.dy = dy;
params.Lx = Lx;
params.Ly = Ly;

% Mobility for dual interaction evolution equation
M0  = mab/(3*eps);

for gamma=gbc

params.kbc = eps*gamma;
params.Obc = K*gamma/eps;
path = [folder,grid,'_gamma_ratio',num2str(gamma/gab)];

%% Initial conditions
phia = zeros(nx,ny); %Initializing phase alpha
phib = zeros(nx,ny); %Initializing phase beta% Initilize fields for laplacians
laplace_a = zeros(nx,ny);
laplace_b = zeros(nx,ny);

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
    laplace_a(i,j) = (phia(i+1,j)-2*phia(i,j)+phia(i-1,j))/(dx*dx) + (phia(i,j+1)-2*phia(i,j)+phia(i,j-1))/(dy*dy);
    laplace_b(i,j) = (phib(i+1,j)-2*phib(i,j)+phib(i-1,j))/(dx*dx) + (phib(i,j+1)-2*phib(i,j)+phib(i,j-1))/(dy*dy);

    %% Multi-phase field evolution compute all cells
    % [dphia, dphib] = calcPhi_grad1_toth(phia,phib,phic,laplace_a,laplace_b,params);
    % [dphia, dphib] = calcPhi_grad2_moelans(phia,phib,phic,laplace_a,laplace_b,params);
    % [dphia, dphib] = calcPhi_grad2_garcke(phia,phib,phic,laplace_a,laplace_b,params);
    % phia(i,j) = phia(i,j)-dt.*M0.*dphia;
    % phib(i,j) = phib(i,j)-dt.*M0.*dphib;
    
    %% Multi-phase field evolution compute only interface cells
    [active, dphia, dphib] = calcPhi_grad1_obstacle_red(phia,phib,phic,laplace_a,laplace_b,params);
    % [active, dphia, dphib] = calcPhi_grad2_obstacle_red(phia,phib,phic,laplace_a,laplace_b,params);
    phia(active) = phia(active)-dt.*M0.*dphia;
    phib(active) = phib(active)-dt.*M0.*dphib;
    
    phic = ones(nx,ny) - phia - phib;
    
    %% Gibbs simplex
    [phia,phib,phic] = calc_gibbs_simplex(phia,phib,phic);

    %% Boundary conditions
    % No-flux BC at left and right sides. Enable for benchmark 2.
    % phia(1,:)  = phia(2,:);
    % phia(nx,:) = phia(nx-1,:);
    % phib(1,:)  = phib(2,:);
    % phib(nx,:) = phib(nx-1,:);
    % phic(1,:)  = phic(2,:);
    % phic(nx,:) = phic(nx-1,:);
    
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
            plot_field(phia,fig_phia,cast(it,'double')*dt,dx,dy);
            plot_xy(x,phia(:,2),fig_ghost,cast(it,'double')*dt,min(phia(:,2)),1);
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
            phia_ref = load([comparepath,'_phia_F',num2str(count),'.dat'],'-ascii');
            phi_diff = phia - phia_ref;
            plot_field(phi_diff,fig_diff,cast(it,'double')*dt,params);
        end
        
        count = count+1;
    end
end
toc

end

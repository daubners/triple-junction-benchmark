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
t_end = 10000.0;   % Total time (non-dimensional)
frames = 10;       % Plot and output frames

folder = ['..\Simulations\Benchmark-2\MOP-grad0-wellmoelans-L0\'];
grid   = ['100x400'];

assignmentStrategy = 'L_const'; % 'k_const' 'L_const' 'Om_const'
potential = 'wellmoelans'; % obstacle welltoth wellgarcke wellmoelans
benchmark = 'steadyState'; % static steadyState

%% Prepare simulation parameters
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

fileID = fopen([folder,grid,'_overall_metrics.dat'],'w');
switch benchmark
    case 'static'
        fprintf(fileID,'%9s \t\t %9s \t\t %9s \t\t %9s \t\t %9s \t\t %9s \t %9s \t %9s \t %9s \t\t %9s\n','gab/g0','theta','Ftot','TP_x','TP_y', 'ghost_a', 'ghost_b', 'ghost_c', 'eps_f', 'eps_TP');
    case 'steadyState'
        fprintf(fileID,'%9s \t\t %9s \t\t %9s \t\t %9s \t %9s \t\t %9s \t %9s \t %9s \t\t %9s \t\t %9s \t %9s \t %9s \t %9s\n','gab/g0','time','theta','h_GB','v_TP','v_max','L2norm','TP_x','TP_y', 'ghost_a', 'ghost_b', 'ghost_c', 'd_h_GB');
    otherwise
        error('Choose grad1 or grad2 for gradient term!')
end

fclose(fileID);

for gamma=gbc
    fig_TP     = figure;
    fig_phis   = figure;
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
    
    for k=1:frames
        phia = load([path,'_phia_F',num2str(k),'.dat'],'phia','-ascii');
        phib = load([path,'_phib_F',num2str(k),'.dat'],'phib','-ascii');
        phic = load([path,'_phic_F',num2str(k),'.dat'],'phic','-ascii');
        plot_field(phib,fig_phis,t_end/frames*k,dx,dy);
        switch benchmark
            case 'static'
                [Ftot(k),TP_x(k),TP_y(k),ghost(k,:)] = calc_metrics_static(assignmentStrategy,potential,phia,phib,phic,fig_TP,params);
            case 'steadyState'
                [Ftot(k),TP_x(k),TP_y(k),max_GB(k),h_GB(k),L2(k),ghost(k,:)] = calc_metrics_steadyState(assignmentStrategy,potential,phia,phib,phic,fig_TP,params);
        end
        
        time(k) = t_end/frames*k;
    end
    %Individual files
    fileID = fopen([path,'_metrics.dat'],'w');
    switch benchmark
        case 'static'
            fprintf(fileID,'%9s \t\t %9s \t\t %9s \t\t %9s \t\t %9s \t %9s \t %9s\n','time','Ftot','TP_x','TP_y', 'ghost_a', 'ghost_b', 'ghost_c');
            M_out=[time; Ftot; TP_x; TP_y; ghost.'];
            fprintf(fileID,'%9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \n',M_out);
        case 'steadyState'
            v_TP  = [0,-diff(TP_y)./(t_end/frames)];
            v_max = [0,-diff(max_GB)./(t_end/frames)];
            conv  = [1,abs(diff(h_GB))./Lx];
            fprintf(fileID,'%9s \t\t %9s \t\t %9s \t\t %9s \t\t %9s \t %9s \t %9s \t %9s \t %9s \t\t %9s \t\t %9s\n','time','Ftot','TP_x','TP_y', 'h_GB', 'ghost_a', 'ghost_b', 'ghost_c', 'v_TP', 'v_max', 'd_h_GB');
            M_out=[time; Ftot; TP_x; TP_y; h_GB; ghost.'; v_TP; v_max; conv];
            fprintf(fileID,'%9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \n',M_out);
    end
    fprintf(fileID,'\n');
    fclose(fileID);

    % Overall file for simulation study
    fileID = fopen([folder,grid,'_overall_metrics.dat'],'a+');
    switch benchmark
        case 'static'
            if TP_x(frames) == -1
                theta = -1;
            else
                theta = 2*atan(0.5*Lx/(180-TP_y(frames)))*180/pi;
            end
            converge1 = abs(Ftot(frames)-Ftot(frames-1))/Ftot(frames);
            converge2 = abs(TP_y(frames)-TP_y(frames-1))/Ly;
            P_out=[gamma; theta; Ftot(frames); TP_x(frames); TP_y(frames); ghost(frames,:).'; converge1; converge2];
            fprintf(fileID,'%9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f\n',P_out);
        case 'steadyState'
            if isempty(min(find(TP_x==-1)))
                last = frames;
            elseif min(find(TP_x==-1))==1
                last = 2;
            else
                last = min(find(TP_x==-1))-1;
            end
            if TP_x(frames) == -1
                theta = -1;
            else
                theta = (7/5*exp(-h_GB(last)/100))^(3.5);
                dtheta = 1;
                while dtheta>1e-10
                    dtheta = (exp(-h_GB(last)/100)*(sin(0.5*theta))^(1/(theta-pi))-1) / (0.5*(pi-theta)*cot(0.5*theta)+log(sin(0.5*theta))) *(pi-theta)^2;
                    theta = theta + dtheta;
                end
                theta = theta*180/pi;
            end
            
            P_out=[gamma; time(last); theta; h_GB(last); v_TP(last); v_max(last); L2(last); TP_x(last); TP_y(last); ghost(last,:).'; conv(last)];
            fprintf(fileID,'%9.5f \t %9.2f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f \t %9.8f\n',P_out);
    end
    fclose(fileID);
end



%% gradient 1 + obstacle formulation
function [active, dphia, dphib] = calcPhi_grad1_obstacle_red(phia,phib,phic,laplace_a,laplace_b,params)
    nx = params.nx;
    ny = params.ny;
    kab = params.kab;
    kac = params.kac;
    kbc = params.kbc;
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
    % Identify cells with laplacian to reduce actively computed cells
    [row,col] = find(laplace_a~=0 | laplace_b~=0);
    active = row+nx*(col-1);
    
    dF_dphia = kab.*laplace_b(active) - kac.*(laplace_a(active)+laplace_b(active)) + Oab.*phib(active) + Oac.*phic(active);
    dF_dphib = kab.*laplace_a(active) - kbc.*(laplace_a(active)+laplace_b(active)) + Oab.*phia(active) + Obc.*phic(active);
    dF_dphic = kac.*laplace_a(active) + kbc.*laplace_b(active) + Oac.*phia(active) + Obc.*phib(active);

    dphia = 2.*dF_dphia - dF_dphib - dF_dphic;
    dphib = 2.*dF_dphib - dF_dphia - dF_dphic;
end
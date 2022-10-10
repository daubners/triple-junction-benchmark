%% Grad_interpolate + well Garcke formulation
function [dphia,dphib,dphic] = calcPhi_k0_garcke(phia,phib,phic,params)
    dx = params.dx;
    dy = params.dy;
    
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
    Oabc = params.Oabc;
    
    i=2:params.nx-1;
    j=2:params.ny-1;
    
    laplace_a(i,j) = (phia(i+1,j)-2*phia(i,j)+phia(i-1,j))/(dx*dx) + (phia(i,j+1)-2*phia(i,j)+phia(i,j-1))/(dy*dy);
    laplace_b(i,j) = (phib(i+1,j)-2*phib(i,j)+phib(i-1,j))/(dx*dx) + (phib(i,j+1)-2*phib(i,j)+phib(i,j-1))/(dy*dy);
    laplace_c(i,j) = (phic(i+1,j)-2*phic(i,j)+phic(i-1,j))/(dx*dx) + (phic(i,j+1)-2*phic(i,j)+phic(i,j-1))/(dy*dy);
    
    dphia = - params.k.*laplace_a(i,j)...
            + 2.*Oab.*phia(i,j).*phib(i,j).^2 + 2.*Oac.*phia(i,j).*phic(i,j).^2 + 2.*Oabc.*phia(i,j).*phib(i,j).^2.*phic(i,j).^2;
    dphib = - params.k.*laplace_b(i,j)...
            + 2.*Oab.*phib(i,j).*phia(i,j).^2 + 2.*Obc.*phib(i,j).*phic(i,j).^2 + 2.*Oabc.*phib(i,j).*phia(i,j).^2.*phic(i,j).^2;
    dphic = - params.k.*laplace_c(i,j)...
            + 2.*Oac.*phic(i,j).*phia(i,j).^2 + 2.*Obc.*phic(i,j).*phib(i,j).^2 + 2.*Oabc.*phic(i,j).*phia(i,j).^2.*phib(i,j).^2;
end
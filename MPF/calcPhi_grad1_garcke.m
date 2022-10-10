%% gradient 1 + well Garcke formulation
function [dphia, dphib] = calcPhi_grad1_garcke(phia,phib,phic,laplace_a,laplace_b,params)
    kab = params.kab;
    kac = params.kac;
    kbc = params.kbc;
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
    Oabc = params.Oabc;
    
    i=2:params.nx-1;
    j=2:params.ny-1;
    
    dF_dphia = kab.*laplace_b(i,j) - kac.*(laplace_a(i,j)+laplace_b(i,j))...
               + 2.*Oab.*phia(i,j).*phib(i,j).^2 + 2.*Oac.*phia(i,j).*phic(i,j).^2 + 2.*Oabc.*phia(i,j).*phib(i,j).^2.*phic(i,j).^2;
    dF_dphib = kab.*laplace_a(i,j) - kbc.*(laplace_a(i,j)+laplace_b(i,j))...
               + 2.*Oab.*phib(i,j).*phia(i,j).^2 + 2.*Obc.*phib(i,j).*phic(i,j).^2 + 2.*Oabc.*phib(i,j).*phia(i,j).^2.*phic(i,j).^2;
    dF_dphic = kac.*laplace_a(i,j) + kbc.*laplace_b(i,j)...
               + 2.*Oac.*phic(i,j).*phia(i,j).^2 + 2.*Obc.*phic(i,j).*phib(i,j).^2 + 2.*Oabc.*phic(i,j).*phia(i,j).^2.*phib(i,j).^2;
           
    dphia = 2.*dF_dphia - dF_dphib - dF_dphic;
    dphib = 2.*dF_dphib - dF_dphia - dF_dphic;
end
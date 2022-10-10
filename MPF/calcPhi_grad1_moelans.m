%% gradient 1 + well Moelans formulation
function [dphia, dphib] = calcPhi_grad1_moelans(phia,phib,phic,laplace_a,laplace_b,params)
    zero = 1e-15;
    nx = params.nx;
    ny = params.ny;
    kab = params.kab;
    kac = params.kac;
    kbc = params.kbc;
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
    i=2:nx-1;
    j=2:ny-1;

    sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
    [row,col] = find(sumSquares<=zero);
    
    bulk = row+(nx-2)*(col-1);
    sumSquares(bulk)=1.0;
    Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;
    
    dF_dphia = kab.*laplace_b(i,j) - kac.*(laplace_a(i,j)+laplace_b(i,j))...
               + 1.5.*Oab.*phia(i,j).*phib(i,j).^2 + 1.5.*Oac.*phia(i,j).*phic(i,j).^2 + 0.5.*Ozero.*(phia(i,j).^3 - phia(i,j))...
               + phia(i,j).*((Oab-Ozero).*phib(i,j).^2+(Oac-Ozero).*phic(i,j).^2)./sumSquares.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4);
    dF_dphib = kab.*laplace_a(i,j) - kbc.*(laplace_a(i,j)+laplace_b(i,j))...
               + 1.5.*Oab.*phib(i,j).*phia(i,j).^2 + 1.5.*Obc.*phib(i,j).*phic(i,j).^2 + 0.5.*Ozero.*(phib(i,j).^3 - phib(i,j))...
               + phib(i,j).*((Oab-Ozero).*phia(i,j).^2+(Obc-Ozero).*phic(i,j).^2)./sumSquares.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4);
    dF_dphic = kac.*laplace_a(i,j) + kbc.*laplace_b(i,j)...
               + 1.5.*Oac.*phic(i,j).*phia(i,j).^2 + 1.5.*Obc.*phic(i,j).*phib(i,j).^2 + 0.5.*Ozero.*(phic(i,j).^3 - phic(i,j))...
               + phic(i,j).*((Oac-Ozero).*phia(i,j).^2+(Obc-Ozero).*phib(i,j).^2)./sumSquares.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4);
           
    dphia = 2.*dF_dphia - dF_dphib - dF_dphic;
    dphib = 2.*dF_dphib - dF_dphia - dF_dphic;
end
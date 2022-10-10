%% Grad_interpolate + welltoth
function [dphia,dphib,dphic] = calcPhi_k0_toth(phia,phib,phic,params)
    zero = 1e-15;
    nx = params.nx;
    ny = params.ny;
    dx = params.dx;
    dy = params.dy;
    
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
    i=2:nx-1;
    j=2:ny-1;
    
    laplace_a = (phia(i+1,j)-2*phia(i,j)+phia(i-1,j))/(dx*dx) + (phia(i,j+1)-2*phia(i,j)+phia(i,j-1))/(dy*dy);
    laplace_b = (phib(i+1,j)-2*phib(i,j)+phib(i-1,j))/(dx*dx) + (phib(i,j+1)-2*phib(i,j)+phib(i,j-1))/(dy*dy);
    laplace_c = (phic(i+1,j)-2*phic(i,j)+phic(i-1,j))/(dx*dx) + (phic(i,j+1)-2*phic(i,j)+phic(i,j-1))/(dy*dy);

    sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
    [row,col] = find(sumSquares<=zero);
    
    bulk = row+(nx-2)*(col-1);
    sumSquares(bulk)=1.0;
    Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;
    L_int = (phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + params.gratio.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;
    L_int(bulk) = 1.0;
    
    dphia = - params.k.*laplace_a...
            + Oab.*phia(i,j).*phib(i,j).^2 + Oac.*phia(i,j).*phic(i,j).^2 + Ozero.*(phia(i,j).^3 - phia(i,j).^2)...
            + 2.*phia(i,j).*((Oab-Ozero).*phib(i,j).^2+(Oac-Ozero).*phic(i,j).^2)./sumSquares.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12);
    dphib = - params.k.*laplace_b...
            + Oab.*phib(i,j).*phia(i,j).^2 + Obc.*phib(i,j).*phic(i,j).^2 + Ozero.*(phib(i,j).^3 - phib(i,j).^2)...
            + 2.*phib(i,j).*((Oab-Ozero).*phia(i,j).^2+(Obc-Ozero).*phic(i,j).^2)./sumSquares.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12);
    dphic = - params.k.*laplace_c...
            + Oac.*phic(i,j).*phia(i,j).^2 + Obc.*phic(i,j).*phib(i,j).^2 + Ozero.*(phic(i,j).^3 - phic(i,j).^2)...
            + 2.*phic(i,j).*((Oac-Ozero).*phia(i,j).^2+(Obc-Ozero).*phib(i,j).^2)./sumSquares.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12);
        
    dphia = dphia.*L_int;
    dphib = dphib.*L_int;
    dphic = dphic.*L_int;
end
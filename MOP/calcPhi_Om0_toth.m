%% Grad_interpolate + welltoth
function [dphia,dphib,dphic] = calcPhi_Om0_toth(phia,phib,phic,params)
    zero = 1e-15;
    nx = params.nx;
    ny = params.ny;
    dx = params.dx;
    dy = params.dy;
    
    kab = params.kab;
    kac = params.kac;
    kbc = params.kbc;
    
    i=2:nx-1;
    j=2:ny-1;

    laplace_a = (phia(i+1,j)-2*phia(i,j)+phia(i-1,j))/(dx*dx) + (phia(i,j+1)-2*phia(i,j)+phia(i,j-1))/(dy*dy);
    laplace_b = (phib(i+1,j)-2*phib(i,j)+phib(i-1,j))/(dx*dx) + (phib(i,j+1)-2*phib(i,j)+phib(i,j-1))/(dy*dy);
    laplace_c = (phic(i+1,j)-2*phic(i,j)+phic(i-1,j))/(dx*dx) + (phic(i,j+1)-2*phic(i,j)+phic(i,j-1))/(dy*dy);

    sumSquares = phia.^2.*phib.^2 + phia.^2.*phic.^2 + phib.^2.*phic.^2;
    [row,col] = find(sumSquares<=zero);
    
    bulk = row+nx*(col-1);
    sumSquares(bulk)=1.0;
    kzero = (kab.*phia.^2.*phib.^2 + kac.*phia.^2.*phic.^2 + kbc.*phib.^2.*phic.^2)./sumSquares;
    kzero(bulk) = kab;
    sumGradSquare =   (phia(i+1,j)-phia(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx) + (phia(i,j+1)-phia(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
                    + (phib(i+1,j)-phib(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx) + (phib(i,j+1)-phib(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
                    + (phic(i+1,j)-phic(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx) + (phic(i,j+1)-phic(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy);
    L_int = (phia.^2.*phib.^2 + phia.^2.*phic.^2 + phib.^2.*phic.^2./params.gratio)./sumSquares;
    L_int(bulk) = 1.0;
    
    dphia = - kzero(i,j).*laplace_a...
            - (kzero(i+1,j)-kzero(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx) - (kzero(i,j+1)-kzero(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
            + phia(i,j).*((kab-kzero(i,j)).*phib(i,j).^2+(kac-kzero(i,j)).*phic(i,j).^2)./sumSquares(i,j).*sumGradSquare...
            + params.O.*phia(i,j).*phib(i,j).^2 + params.O.*phia(i,j).*phic(i,j).^2 + params.O.*(phia(i,j).^3 - phia(i,j).^2);
    dphib = - kzero(i,j).*laplace_b...
            - (kzero(i+1,j)-kzero(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx) - (kzero(i,j+1)-kzero(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
            + phib(i,j).*((kab-kzero(i,j)).*phia(i,j).^2+(kbc-kzero(i,j)).*phic(i,j).^2)./sumSquares(i,j).*sumGradSquare...
            + params.O.*phib(i,j).*phia(i,j).^2 + params.O.*phib(i,j).*phic(i,j).^2 + params.O.*(phib(i,j).^3 - phib(i,j).^2);
    dphic = - kzero(i,j).*laplace_c...
            - (kzero(i+1,j)-kzero(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx) - (kzero(i,j+1)-kzero(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy)...
            + phic(i,j).*((kac-kzero(i,j)).*phia(i,j).^2+(kbc-kzero(i,j)).*phib(i,j).^2)./sumSquares(i,j).*sumGradSquare...
            + params.O.*phic(i,j).*phia(i,j).^2 + params.O.*phic(i,j).*phib(i,j).^2 + params.O.*(phic(i,j).^3 - phic(i,j).^2);
                          
    dphia = dphia.*L_int(i,j);
    dphib = dphib.*L_int(i,j);
    dphic = dphic.*L_int(i,j);
end
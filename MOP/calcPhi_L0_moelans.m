%% Grad_interpolate + wellmoelans
function [dphia,dphib,dphic] = calcPhi_L0_moelans(phia,phib,phic,params)
    zero = 1e-15;
    nx = params.nx;
    ny = params.ny;
    dx = params.dx;
    dy = params.dy;
    
    kab = params.kab;
    kac = params.kac;
    kbc = params.kbc;
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
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
    Ozero = (Oab.*phia.^2.*phib.^2 + Oac.*phia.^2.*phic.^2 + Obc.*phib.^2.*phic.^2)./sumSquares;
    kzero(bulk) = kab;
    Ozero(bulk) = Oab;
    sumGradSquare =   (phia(i+1,j)-phia(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx) + (phia(i,j+1)-phia(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
                    + (phib(i+1,j)-phib(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx) + (phib(i,j+1)-phib(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
                    + (phic(i+1,j)-phic(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx) + (phic(i,j+1)-phic(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy);
    
    dphia = - kzero(i,j).*laplace_a...
            - (kzero(i+1,j)-kzero(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx) - (kzero(i,j+1)-kzero(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
            + phia(i,j).*((kab-kzero(i,j)).*phib(i,j).^2+(kac-kzero(i,j)).*phic(i,j).^2)./sumSquares(i,j).*sumGradSquare...
            + 1.5.*Oab.*phia(i,j).*phib(i,j).^2 + 1.5.*Oac.*phia(i,j).*phic(i,j).^2 + 0.5.*Ozero(i,j).*(phia(i,j).^3 - phia(i,j))...
            + phia(i,j).*((Oab-Ozero(i,j)).*phib(i,j).^2+(Oac-Ozero(i,j)).*phic(i,j).^2)./sumSquares(i,j).*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4);
    dphib = - kzero(i,j).*laplace_b...
            - (kzero(i+1,j)-kzero(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx) - (kzero(i,j+1)-kzero(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
            + phib(i,j).*((kab-kzero(i,j)).*phia(i,j).^2+(kbc-kzero(i,j)).*phic(i,j).^2)./sumSquares(i,j).*sumGradSquare...
            + 1.5.*Oab.*phib(i,j).*phia(i,j).^2 + 1.5.*Obc.*phib(i,j).*phic(i,j).^2 + 0.5.*Ozero(i,j).*(phib(i,j).^3 - phib(i,j))...
            + phib(i,j).*((Oab-Ozero(i,j)).*phia(i,j).^2+(Obc-Ozero(i,j)).*phic(i,j).^2)./sumSquares(i,j).*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4);
    dphic = - kzero(i,j).*laplace_c...
            - (kzero(i+1,j)-kzero(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx) - (kzero(i,j+1)-kzero(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy)...
            + phic(i,j).*((kac-kzero(i,j)).*phia(i,j).^2+(kbc-kzero(i,j)).*phib(i,j).^2)./sumSquares(i,j).*sumGradSquare...
            + 1.5.*Oac.*phic(i,j).*phia(i,j).^2 + 1.5.*Obc.*phic(i,j).*phib(i,j).^2 + 0.5.*Ozero(i,j).*(phic(i,j).^3 - phic(i,j))...
            + phic(i,j).*((Oac-Ozero(i,j)).*phia(i,j).^2+(Obc-Ozero(i,j)).*phib(i,j).^2)./sumSquares(i,j).*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4);
end
%% Calc integral of total interfacial energy for MOP models
function [Ftot] = calc_Ftot(assignmentStrategy,potential,phia,phib,phic,params)
    i=2:params.nx-1;
    j=2:params.ny-1;
    dx = params.dx;
    dy = params.dy;
    
    switch assignmentStrategy
        case 'k_const'
        Fgrad = 0.5*sum(sum(   params.k.*(phia(i+1,j)-phia(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx)...
                             + params.k.*(phia(i,j+1)-phia(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
                             + params.k.*(phib(i+1,j)-phib(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx)...
                             + params.k.*(phib(i,j+1)-phib(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
                             + params.k.*(phic(i+1,j)-phic(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx)...
                             + params.k.*(phic(i,j+1)-phic(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy)))*dx*dy;
        case 'L_const'
        sumSquares = phia.^2.*phib.^2 + phia.^2.*phic.^2 + phib.^2.*phic.^2;
        [row,col] = find(sumSquares<=1e-15);
    
        bulk = row+params.nx*(col-1);
        sumSquares(bulk)=1.0;
        kzero = (params.kab.*phia.^2.*phib.^2 + params.kac.*phia.^2.*phic.^2 + params.kbc.*phib.^2.*phic.^2)./sumSquares;
        kzero(bulk) = params.kab;
        sumGradSquare =   (phia(i+1,j)-phia(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx) + (phia(i,j+1)-phia(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
                        + (phib(i+1,j)-phib(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx) + (phib(i,j+1)-phib(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
                        + (phic(i+1,j)-phic(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx) + (phic(i,j+1)-phic(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy);
        Fgrad = 0.5*sum(sum( kzero(i,j).*sumGradSquare ))*dx*dy;
        case 'Om_const'
        sumSquares = phia.^2.*phib.^2 + phia.^2.*phic.^2 + phib.^2.*phic.^2;
        [row,col] = find(sumSquares<=1e-15);
    
        bulk = row+params.nx*(col-1);
        sumSquares(bulk)=1.0;
        kzero = (params.kab.*phia.^2.*phib.^2 + params.kac.*phia.^2.*phic.^2 + params.kbc.*phib.^2.*phic.^2)./sumSquares;
        kzero(bulk) = params.kab;
        sumGradSquare =   (phia(i+1,j)-phia(i-1,j)).*(phia(i+1,j)-phia(i-1,j))./(2*dx*2*dx) + (phia(i,j+1)-phia(i,j-1)).*(phia(i,j+1)-phia(i,j-1))./(2*dy*2*dy)...
                        + (phib(i+1,j)-phib(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx) + (phib(i,j+1)-phib(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
                        + (phic(i+1,j)-phic(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx) + (phic(i,j+1)-phic(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy);
        Fgrad = 0.5*sum(sum( kzero(i,j).*sumGradSquare ))*dx*dy;
        otherwise
        error('Gradient term is only defined with const for now!')
    end
    
    switch potential
        case 'welltoth'
        switch assignmentStrategy
            case 'k_const'
            Oab = params.Oab;
            Oac = params.Oac;
            Obc = params.Obc;
            sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
            [row,col] = find(sumSquares<=1e-15);

            bulk = row+(params.nx-2)*(col-1);
            sumSquares(bulk)=1.0;
            Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;

            Fpot  = sum(sum(  0.5.*Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                             +0.5.*Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                             +0.5.*Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                             +Ozero.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12) ))*dx*dy;
            case 'L_const'
            Oab = params.Oab;
            Oac = params.Oac;
            Obc = params.Obc;
            sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
            [row,col] = find(sumSquares<=1e-15);

            bulk = row+(params.nx-2)*(col-1);
            sumSquares(bulk)=1.0;
            Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;

            Fpot  = sum(sum(  0.5.*Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                             +0.5.*Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                             +0.5.*Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                             +Ozero.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12) ))*dx*dy;
            case 'Om_const'
            Fpot  = sum(sum(  0.5.*params.O.*phia(i,j).^2.*phib(i,j).^2 ...
                             +0.5.*params.O.*phia(i,j).^2.*phic(i,j).^2 ...
                             +0.5.*params.O.*phib(i,j).^2.*phic(i,j).^2 ...
                             +params.O.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12) ))*dx*dy;
        end
        case 'wellmoelans'
        switch assignmentStrategy
            case 'k_const'
            Oab = params.Oab;
            Oac = params.Oac;
            Obc = params.Obc;
            sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
            [row,col] = find(sumSquares<=1e-15);

            bulk = row+(params.nx-2)*(col-1);
            sumSquares(bulk)=1.0;
            Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;

            Fpot  = sum(sum(  0.75.*Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                             +0.75.*Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                             +0.75.*Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                             +0.5.*Ozero.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4) ))*dx*dy;
            case 'L_const'
            Oab = params.Oab;
            Oac = params.Oac;
            Obc = params.Obc;
            sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
            [row,col] = find(sumSquares<=1e-15);

            bulk = row+(params.nx-2)*(col-1);
            sumSquares(bulk)=1.0;
            Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;

            Fpot  = sum(sum(  0.75.*Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                             +0.75.*Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                             +0.75.*Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                             +0.5.*Ozero.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4) ))*dx*dy;
            case 'Om_const'
            Fpot  = sum(sum(  0.75.*params.O.*phia(i,j).^2.*phib(i,j).^2 ...
                             +0.75.*params.O.*phia(i,j).^2.*phic(i,j).^2 ...
                             +0.75.*params.O.*phib(i,j).^2.*phic(i,j).^2 ...
                             +0.5.*params.O.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4) ))*dx*dy;
        end
        otherwise
        error('Choose welltoth or wellmoelans for potential term!')
    end
    
    Ftot = Fgrad + Fpot;
end
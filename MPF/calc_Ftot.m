%% Calc integral of total interfacial energy for MPF models
function [Ftot] = calc_Ftot(gradient,potential,phia,phib,phic,params)
    i=2:params.nx-1;
    j=2:params.ny-1;
    dx = params.dx;
    dy = params.dy;
    Oab = params.Oab;
    Oac = params.Oac;
    Obc = params.Obc;
    
    switch gradient
        case 'grad1'
        Fgrad = sum(sum( -params.kab.*(phia(i+1,j)-phia(i-1,j)).*(phib(i+1,j)-phib(i-1,j))./(2*dx*2*dx)...
                         -params.kab.*(phia(i,j+1)-phia(i,j-1)).*(phib(i,j+1)-phib(i,j-1))./(2*dy*2*dy)...
                         -params.kac.*(phia(i+1,j)-phia(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx)...
                         -params.kac.*(phia(i,j+1)-phia(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy)...
                         -params.kbc.*(phib(i+1,j)-phib(i-1,j)).*(phic(i+1,j)-phic(i-1,j))./(2*dx*2*dx)...
                         -params.kbc.*(phib(i,j+1)-phib(i,j-1)).*(phic(i,j+1)-phic(i,j-1))./(2*dy*2*dy)))*dx*dy;
        case 'grad2'
        Fgrad = sum(sum(  params.kab.*( phib(i,j).*(phia(i+1,j)-phia(i-1,j))./(2*dx) - phia(i,j).*(phib(i+1,j)-phib(i-1,j))./(2*dx) ).^2 ...
                         +params.kab.*( phib(i,j).*(phia(i,j+1)-phia(i,j-1))./(2*dy) - phia(i,j).*(phib(i,j+1)-phib(i,j-1))./(2*dy) ).^2 ...
                         +params.kac.*( phic(i,j).*(phia(i+1,j)-phia(i-1,j))./(2*dx) - phia(i,j).*(phic(i+1,j)-phic(i-1,j))./(2*dx) ).^2 ...
                         +params.kac.*( phic(i,j).*(phia(i,j+1)-phia(i,j-1))./(2*dy) - phia(i,j).*(phic(i,j+1)-phic(i,j-1))./(2*dy) ).^2 ...
                         +params.kbc.*( phic(i,j).*(phib(i+1,j)-phib(i-1,j))./(2*dx) - phib(i,j).*(phic(i+1,j)-phic(i-1,j))./(2*dx) ).^2 ...
                         +params.kbc.*( phic(i,j).*(phib(i,j+1)-phib(i,j-1))./(2*dy) - phib(i,j).*(phic(i,j+1)-phic(i,j-1))./(2*dy) ).^2 ))*dx*dy;
        otherwise
        error('Choose grad1 or grad2 for gradient term!')
    end
    
    switch potential
        case 'welltoth'
        sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
        [row,col] = find(sumSquares<=1e-15);

        bulk = row+(params.nx-2)*(col-1);
        sumSquares(bulk)=1.0;
        Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;

        Fpot  = sum(sum(  0.5.*Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                         +0.5.*Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                         +0.5.*Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                         +Ozero.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^3 +phib(i,j).^3 +phic(i,j).^3)./3 + 1/12) ))*dx*dy;
        case 'wellmoelans'
        sumSquares = phia(i,j).^2.*phib(i,j).^2 + phia(i,j).^2.*phic(i,j).^2 + phib(i,j).^2.*phic(i,j).^2;
        [row,col] = find(sumSquares<=1e-15);

        bulk = row+(params.nx-2)*(col-1);
        sumSquares(bulk)=1.0;
        Ozero = (Oab.*phia(i,j).^2.*phib(i,j).^2 + Oac.*phia(i,j).^2.*phic(i,j).^2 + Obc.*phib(i,j).^2.*phic(i,j).^2)./sumSquares;

        Fpot  = sum(sum(  0.75.*Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                         +0.75.*Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                         +0.75.*Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                         +0.5.*Ozero.*((phia(i,j).^4 + phib(i,j).^4 + phic(i,j).^4)./4 - (phia(i,j).^2 +phib(i,j).^2 +phic(i,j).^2)./2 + 1/4) ))*dx*dy;
        case 'wellgarcke'
        Oabc = params.Oabc;
        Fpot  = sum(sum(  Oab.*phia(i,j).^2.*phib(i,j).^2 ...
                         +Oac.*phia(i,j).^2.*phic(i,j).^2 ...
                         +Obc.*phib(i,j).^2.*phic(i,j).^2 ...
                         +Oabc.*phia(i,j).^2.*phib(i,j).^2.*phic(i,j).^2 ))*dx*dy;
        case 'obstacle'
        Fpot  = sum(sum(  Oab.*phia(i,j).*phib(i,j)...
                         +Oac.*phia(i,j).*phic(i,j)...
                         +Obc.*phib(i,j).*phic(i,j) ))*dx*dy;
        otherwise
        error('Choose welltoth, wellgarcke or obstacle for potential term!')
    end
    
    Ftot = Fgrad + Fpot;
end
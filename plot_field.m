%% Plot function for fields
function plot_field(field,fig_name,time,dx,dy)
    [nx,ny] = size(field);
    Lx = (nx-2)*dx;
    Ly = (ny-2)*dy;
    x = -dx/2:dx:Lx+dx/2;
    y = -dy/2:dy:Ly+dy/2;

    figure(fig_name);
    g=pcolor(x,y,(field).');
    pbaspect([Lx Ly 1]);
    axis([0 Lx 0 Ly]);
    shading interp;
    colorbar;
    colormap turbo; %parula; %turbo; hot; spring; summer; gray;
    title({['Field solution at time t = ',num2str(time)]})
    xlabel('x-axis')
    ylabel('y-axis')
    drawnow;
end
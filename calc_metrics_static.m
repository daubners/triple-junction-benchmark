%% Calc metrics for static triple junction
% Integral of total energy should be interfacial energy
% Position of triple junction defines angles
% Check spurious phases
function [Ftot,TP_x,TP_y,ghost] = calc_metrics_static(gradient,potential,phia,phib,phic,fig_name,params)
    i=2:params.nx-1;
    j=2:params.ny-1;
    dx = params.dx;
    dy = params.dy;
    
    %% Calculate overall energy
    Ftot = calc_Ftot(gradient,potential,phia,phib,phic,params);
    
    %% Triple point position
    x = -dx/2:dx:params.Lx+dx/2;
    y = -dy/2:dy:params.Ly+dy/2;
    
    iso_line = contourc(x(i),y(j),(phia(i,j)-phib(i,j)).',[0.0 0.0]);
    iso_ab = iso_line(:,2:size(iso_line,2));
    iso_line = contourc(x(i),y(j),(phia(i,j)-phic(i,j)).',[0.0 0.0]);
    iso_ac = iso_line(:,2:size(iso_line,2));
    iso_line = contourc(x(i),y(j),(phib(i,j)-phic(i,j)).',[0.0 0.0]);
    iso_bc = iso_line(:,2:size(iso_line,2));
    
    [TP_x,TP_y] = calc_intersection(iso_ab(1,:),iso_ab(2,:),iso_ac(1,:),iso_ac(2,:));
    % These two are for validation and actually never used
    [TP_x2,TP_y2] = calc_intersection(iso_ab(1,:),iso_ab(2,:),iso_bc(1,:),iso_bc(2,:));
    [TP_x3,TP_y3] = calc_intersection(iso_ac(1,:),iso_ac(2,:),iso_bc(1,:),iso_bc(2,:));
    
    if isempty(TP_x)
        if isempty(TP_y)
            TP_x = -1;
            TP_y = -1;
        else
            error('very strange things happening...')
        end
    end
    
    if size(TP_x,1)>1
        [val,idx]=min(abs(TP_x-50));
        TP_x = TP_x(idx);
        TP_y = TP_y(idx);
    end
    
    % Figure is for instant validation
    figure(fig_name);
    g=pcolor(x,y,(phia).');
    pbaspect([params.Lx params.Ly 1]);
    axis([0 params.Lx 0 params.Ly]);
    shading interp;
    colorbar;
    colormap turbo;
    title({['Triple junction']});
    xlabel('x-axis')
    ylabel('y-axis')
    hold on
        h1=contour(x,y,(phia-phib).',[0.0,0.0],'blue',LineWidth=2);
        h2=contour(x,y,(phia-phic).',[0.0,0.0],'red',LineWidth=2);
        h3=contour(x,y,(phib-phic).',[0.0,0.0],'green',LineWidth=2);
        plot(TP_x,TP_y,'v','MarkerSize',10,'MarkerEdgeColor','black');
    hold off
    drawnow;
    
    %% Max value of Ghost phases in binary interfaces
    % phia evaluated at lower boundary
    % phib at vertical line at 0.75*L and
    % phic at vertical line at 0.25*L.
    ghost = [max(phia(:,2)) max(phib(77,:)) max(phic(26,:))];
end
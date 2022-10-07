%% Calc metrics for calc_metrics_steadyState triple junction
% Integral of total energy should be interfacial energy
% Position of triple junction defines angles
% Check spurious phases
function [Ftot,TP_x,TP_y,max_GB,h_GB,L2,ghost] = calc_metrics_steadyState(gradient,potential,phia,phib,phic,fig_name,params)
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
    
    % If multiple intersections are computed due to numerical artefacts
    % select the one in the middle (others usually occur at the domain
    % boundaries).
    if size(TP_x,1)>1
        [val,idx]=min(abs(TP_x-50));
        TP_x = TP_x(idx);
        TP_y = TP_y(idx);
    end
    
    %% Maximum of curved profile
    max_GB = max(iso_ab(2,:));
    if isempty(max_GB)
        max_GB = -1;
    end
    % As the problem is symmetric this value should be the same
    % max(iso_ac(2,:))
    
    h_GB = max(iso_ab(2,:)) - TP_y;
    if isempty(h_GB)
        h_GB = -1;
    end
    
    %% L2-norm of curved profile
    shape_analytic = 100/(pi-2*acos(0.5*params.gratio))*log( cos((pi-2*acos(0.5*params.gratio))/100*x(2:params.nx-1)) );
    shape_analytic(params.nx/2:end) = shape_analytic((params.nx-2)/2:-1:1);
    
    if TP_x == -1
        % If triple point is ill-defined, assume TP_x=50.0.
        shape_numeric = [iso_ab(:,find(iso_ab(1,:)<50-0.25*dx)),iso_ac(:,find(iso_ac(1,:)>50+0.25*dx))];
    else
        % If triple point is well-defined, use actual TP_x.
        shape_numeric = [iso_ab(:,find(iso_ab(1,:)<TP_x-0.25*dx)),iso_ac(:,find(iso_ac(1,:)>TP_x+0.25*dx))];
    end
    shape_numeric = sortrows(shape_numeric')';
    y_numeric = interp1(shape_numeric(1,:),shape_numeric(2,:)-max_GB,x(2:params.nx-1));
    L2 = sqrt( sum((y_numeric-shape_analytic).^2)/size(shape_analytic,2) )/params.Lx;

    % Figure for debugging errors with L2 norm
    % figure;
    % h=plot(x(2:params.nx-1),shape_analytic,'r');
    % xlim([x(1) x(end)])
    % title(['Boundary profile'])
    % hold on
    %     h1=plot(iso_ab(1,:),iso_ab(2,:)- max_GB,'blue');
    %     h2=plot(iso_ac(1,:),iso_ac(2,:)- max_GB,'green');
    %     h3=plot(shape_numeric(1,:),shape_numeric(2,:)-max_GB,'black');
    %     plot(x(2:params.nx-1), y_numeric, '.', 'markersize', 8);
    % hold off
    % drawnow;
    
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
    % phib at vertical line at L and
    % phic at vertical line at 0.
    ghost = [max(phia(:,2)) max(phib(101,:)) max(phic(2,:))];
end
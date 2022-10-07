%% Plot function for scalar metrics
function plot_xy(x,y,fig_name,time,low,high)
    figure(fig_name);
    h=plot(x,y,'r');
    xlim([x(1) x(end)])
    ylim([low high])
    title(['xy data at time t = ',num2str(time)])
    refreshdata(h)
    drawnow;
end
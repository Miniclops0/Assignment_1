function [] = plot_traj(colour, xp, x, yp, y, box)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:10       
    x_o = xp(i,:);
    y_o = yp(i,:);
    x_n = x(i,:);
    y_n = y(i,:);
    figure(3)
    if box == 1
        rectangle('position', [0.8e-7 0.6e-7 0.4e-7  0.4e-7])
        hold on
        rectangle('position', [0.8e-7 0 0.4e-7  0.4e-7])
        hold on
    end
    plot([x_o x_n], [y_o y_n], 'color', colour(i,:), 'linewidth', 2)
    xlim([0, 200e-9])
    ylim([0, 100e-9])
%     title(['Semiconductor Temp ', num2str(avg_temp), ' K'])
    title('Electron Trajectories')
    xlabel('x (m)')
    ylabel('y (m)')
    hold on
%     pause(0.01)
end
%
end


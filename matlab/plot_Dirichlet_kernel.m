% Dirichlet kernel example
N = 16;
% basic setting
% clear
close all;
clc
% basic setting
% bueatiful color scheme
colors = ...
 1/256*[ 31,119,180; % 1 默认蓝色
        255,127,14;  % 2 橘色
         44,160,44;  % 3 绿色
        214,39,40;   % 4 红色
        148,103,189; % 5 紫色
        140,86,75;   % 6 棕色 
       227,119,194;  % 7 粉色
       127,127,127;  % 8 灰色
       188,189,34;   % 9 青棕
       23,190,207;   % 10 淡蓝
       26,85,255;    % 11 鲜蓝色
       ];

% marker list
markers = {'o','s','d','^','v','>','<','p','h','+','x','*'};

% line style
linestyles = {'-','--',':','-.'};
% linewidth
linewidth = 1.;

% marker size
markersize = 8;

% set global font and fontsize
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',10);
% set latex as default text interpreter 
set(0,'DefaultTextInterpreter','latex');

% generate the data
x = linspace(-pi,pi,1000);
% Dirichlet kernel , consider the case x = 0


% plot the Dirichlet kernel
figure(1)
Dirichlet_kernel = plot(x,D(x,N),'k','LineWidth',1, 'Color', colors(1,:));
hold on;
% vertical line stand for the x-axis
plot([-pi pi],[0 0],'k','LineWidth',1);
hold on;
% plot points with markers(1) and the corresponding line vertical to x-axis
% plot the points
points_x = linspace(-pi,pi,17);
points_y = D(points_x,N);

% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(3,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points
Ideal_sampling_without_offset = plot(points_x,points_y,'o','MarkerSize',5,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:));
hold on;


% add offset
offset = 0.1;
points_x_offset = points_x + offset;
points_y_offset = D(points_x_offset,N);

% plot the vertical line
for i = 1:17
    plot([points_x_offset(i) points_x_offset(i)],[0 points_y_offset(i)],'LineWidth',0.5,'Color',colors(4,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points
Discrete_samoling_error = plot(points_x_offset, points_y_offset,'o','MarkerSize',5,'MarkerFaceColor',colors(4,:),'MarkerEdgeColor',colors(4,:));

xlim([-pi pi]);
ylim([-0.4 1.1]);
xlabel('$x$ (rad)');
ylabel('$D_N(x)$');

% legend
legend([Dirichlet_kernel,Ideal_sampling_without_offset,Discrete_samoling_error],{'Dirichlet kernel ($N = 16$)','Ideal sampling','Sampling with offset'},'Location','northwest', 'Interpreter','latex');
title('(d) The cause of power tail','Interpreter','latex');
grid on;

% D = @(x,N) Dirichlet(x,N);
function y = D(x, N)
    if abs(x) < eps
        y = 1;
    else
        y = 1/N * sin(N*x/2)./sin(x/2);
    end
    % turn nan to 1
    y(isnan(y)) = 1;
end




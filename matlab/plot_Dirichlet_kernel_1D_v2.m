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
subplot(2,1,1)
Dirichlet_kernel = plot(x,D(x,16, 0),'k','LineWidth',1, 'Color', colors(1,:));
hold on;
% 给从 -0.25 到 0.25 的 Dirichlet_kernel 的积分区域上色
% 定义要上色区域的边界
x_fill = linspace(-pi/16, pi/16, 500); % 从 -0.25 到 0.25，你可以根据需要调整精度
y_fill = D(x_fill,16,0); % 计算对应的函数值

% 在指定区域内填充颜色
fill_zone = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(3,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;
% vertical line stand for the x-axis
plot([-pi pi],[0 0],'k','LineWidth',1);
hold on;
% plot points with markers(1) and the corresponding line vertical to x-axis
% plot the points
points_x = linspace(-pi,pi,17);
points_y = D(points_x,16, 0);

% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(3,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points
Ideal_sampling_without_offset = plot(points_x,points_y,'o','MarkerSize',5,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:));
hold on;

% twice stepping
points_x_2_step = linspace(-pi,pi,33);
points_y_2_step = D(points_x_2_step,16, 0);

% plot the vertical line
for i = 1:33
    plot([points_x_2_step(i) points_x_2_step(i)],[0 points_y_2_step(i)],'LineWidth',0.5,'Color',colors(4,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points no face only edge
Twice_stepping = plot(points_x_2_step,points_y_2_step,'o','MarkerSize',10,'MarkerEdgeColor',colors(4,:));


xlim([-0.5 1.5]);
ylim([-0.3 1.4]);
xlabel('$x$ (rad)');
ylabel('$D_N(x)$');

% legend
legend([Dirichlet_kernel,Ideal_sampling_without_offset,Twice_stepping, fill_zone],{'Dirichlet kernel ($N = 16$)','Ideal sampling','Angular-domain sampling with $\delta = {\lambda}/{4}$', 'Integral region corresponding to $(\sigma_l^f)^2$'},'Location','northeast', 'Interpreter','latex');
title('(a) Cause of the power leakage: oversampling.');

% add text at [0, 1], with latex interpreter, color set to 3
text(-0.25,1.2,'$(\sigma_l^f)^2 = \int_{\Omega_l (\cdot)} A^2 (\cdot)$','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);

text(0.02,0.5,'$(\sigma_l^f)^2$','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);

text(-0.03, - 0.12,'$\Omega_l$','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);


% arrow with x:0.33,0.33, y:0.772,0.846, set transparency to 0.8
annotation(figure(1),'arrow',[0.34 0.34],[0.772 0.849],'LineWidth',0.8, 'Color', colors(8,:));
annotation(figure(1),'arrow',[0.34 0.34],[0.7149,0.6492],'LineWidth',0.8, 'Color', colors(8,:));

annotation(figure(1),'arrow',[0.303,0.2494],[0.627,0.627],'LineWidth',0.8, 'Color', colors(8,:));
annotation(figure(1),'arrow',[0.3458,0.3982],[0.627,0.627],'LineWidth',0.8, 'Color', colors(8,:));


grid on;

subplot(2,1,2)
offset_value = 0.1;
Dirichlet_kernel = plot(x,D(x,16, 0),'k','LineWidth',1, 'Color', colors(1,:));
hold on;
Dirichlet_kernel_with_offset = plot(x,D(x,16, offset_value),'k','LineWidth',1, 'Color', colors(2,:), 'LineStyle',linestyles{2});
hold on;

% vertical line stand for the x-axis
plot([-pi pi],[0 0],'k','LineWidth',1);
hold on;

% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(3,:), 'LineStyle',linestyles{2});
    hold on;
    plot([points_x_with_offset(i) points_x_with_offset(i)],[0 points_y_with_offset(i)],'LineWidth',0.5,'Color',colors(4,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points
Ideal_sampling_without_offset = plot(points_x,points_y,'o','MarkerSize',5,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:));
hold on;
Ideal_sampling_with_offset = plot(points_x_with_offset,points_y_with_offset,'o','MarkerSize',5,'MarkerFaceColor',colors(4,:),'MarkerEdgeColor',colors(4,:));

xlim([-0.5 1.5]);
ylim([-0.3 1.4]);
xlabel('$x$ (rad)');
ylabel('$D_N(x)$ \& $D_N(x - {\rm offset})$');
grid on;
legend([Dirichlet_kernel,Ideal_sampling_without_offset,Dirichlet_kernel_with_offset,Ideal_sampling_with_offset],{'Dirichlet kernel ($N = 16$)','Ideal sampling','Dirichlet kernel ($N = 16$) with offset','Ideal sampling with offset'},'Location','northeast', 'Interpreter','latex');
title('(b) Cause of the power leakage: offset.');
% D = @(x,N) Dirichlet(x,N);
function y = D(x, N, offset)
    x = x - offset;
    if abs(x) < eps
        y = 1;
    else
        y = 1/N * sin(N*x/2)./sin(x/2);
    end
    % turn nan to 1
    y(isnan(y)) = 1;
end




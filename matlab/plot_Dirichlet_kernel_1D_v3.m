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
subplot(2,1,2)
points_x = linspace(-pi,pi,17);
Dirichlet_kernel_1 = plot(x,D(x,16, 0),'k','LineWidth',1, 'Color', colors(1,:), 'LineStyle', linestyles{2});
hold on;
% vertical line stand for the x-axis
plot([-pi pi],[0 0],'k','LineWidth',1);
hold on;
x_fill = linspace(-pi/16, pi/16, 500); 
y_fill = D(x_fill,16,0); 
fill_zone_1 = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(1,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;

% plot sampling points
points_y = D(points_x,16, 0);
% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(1,:), 'LineStyle',linestyles{2});
    hold on;
end


Dirichlet_kernel_2 = plot(x,D(x,16, pi/8 * 1),'k','LineWidth',1, 'Color', colors(2,:), 'LineStyle', linestyles{2});
hold on;
x_fill = linspace(-pi/16, pi/16, 500)+ pi/8 * 1; 
y_fill = D(x_fill,16,pi/8 * 1); 
fill_zone_2 = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(2,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;

% plot sampling points
points_y = D(points_x,16, pi/8 * 1);
% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(2,:), 'LineStyle',linestyles{2});
    hold on;
end

Dirichlet_kernel_3 = plot(x,D(x,16, pi/8 * 2),'k','LineWidth',1, 'Color', colors(3,:), 'LineStyle', linestyles{2});
hold on;
x_fill = linspace(-pi/16, pi/16, 500)+ pi/8 * 2; 
y_fill = D(x_fill,16,pi/8 * 2); 
fill_zone_3 = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(3,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;

% plot sampling points
points_y = D(points_x,16, pi/8 * 2);
% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(3,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points
points_y = D(points_x,16, 0);
Sampling_1 = plot(points_x,points_y,'o','MarkerSize',5,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
hold on;

points_y = D(points_x,16, pi/8 * 1);
Sampling_2 = plot(points_x,points_y,'o','MarkerSize',8,'MarkerEdgeColor',colors(2,:));

points_y = D(points_x,16, pi/8 * 2);
Sampling_3 = plot(points_x,points_y,'o','MarkerSize',11,'MarkerEdgeColor',colors(3,:));

legend([Dirichlet_kernel_1, Dirichlet_kernel_2, Dirichlet_kernel_3, Sampling_1, Sampling_2, Sampling_3], ...
    {'$l_x = 0$', '$l_x = 1$', '$l_x = 2$', '$n_x^\prime = 0$', '$n_x^\prime = 1$', '$n_x^\prime = 2$'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

xlim([-0.5 2]);
ylim([-0.3 1.4]);
xlabel('$\gamma = \frac{2\pi}{N_x} (l_x - n_x^\prime)$');
ylabel('$D_{N_x}(\gamma)$');
grid on;
title('(b) Equavalent Dirichlet kernel sampling in angular domain with $\delta = 1/2 \lambda$', 'Interpreter', 'latex'); 

figure(1);
text(-0.161235408560312,1.15,0,'$(\sigma_0^f)^2 = \int_{\Omega_0 } A^2 $','Color',colors(1,:),'Interpreter','latex', 'FontSize',12);
% annotation(figure(1),'arrow',[0.2381,0.2751],[0.8552,0.8484],'LineWidth',0.8, 'Color', colors(1,:));

text(0.227,1.15,0,'$(\sigma_1^f)^2 = \int_{\Omega_1 } A^2 $','Color',colors(2,:),'Interpreter','latex', 'FontSize',12);
% annotation(figure(1),'arrow',[0.3173,0.3504],[0.8614,0.8523],'LineWidth',0.8, 'Color', colors(2,:));

text(0.614542801556419,1.15,'$(\sigma_2^f)^2 = \int_{\Omega_2 } A^2 $','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);
% annotation(figure(1),'arrow',[0.4574,0.4152],[0.8653,0.851],'LineWidth',0.8, 'Color', colors(3,:));

text(-0.027480544747082,-0.18,0,'$\Omega_0$','Color',colors(1,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'doublearrow',[0.3451,0.2275],[0.16,0.16],'LineWidth',0.8, 'Color', colors(1,:));

text(0.361624513618677,-0.18,0,'$\Omega_1$','Color',colors(2,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'doublearrow',[0.4664,0.3488],[0.16,0.16],'LineWidth',0.8, 'Color', colors(2,:));

text(0.76045719844358,-0.18,0,'$\Omega_2$','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'doublearrow',[0.5878,0.4701],[0.16,0.16],'LineWidth',0.8, 'Color', colors(3,:));


subplot(2,1,1)
points_x = linspace(-pi,pi,17);
Dirichlet_kernel_1 = plot(x,D(x,16, 0),'k','LineWidth',1, 'Color', colors(1,:), 'LineStyle', linestyles{2});
hold on;
% vertical line stand for the x-axis
plot([-pi pi],[0 0],'k','LineWidth',1);
hold on;
x_fill = linspace(-pi/32, pi/32, 500);
y_fill = D(x_fill,16,0);
fill_zone_1 = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(1,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;

% plot sampling points
points_y = D(points_x,16, 0);
% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(1,:), 'LineStyle',linestyles{2});
    hold on;
end


Dirichlet_kernel_2 = plot(x,D(x,16, pi/8 * 0.5),'k','LineWidth',1, 'Color', colors(2,:), 'LineStyle', linestyles{2});
hold on;
x_fill = linspace(-pi/32, pi/32, 500)+ pi/8 * 0.5;   
y_fill = D(x_fill,16,pi/8 * 0.5);
fill_zone_2 = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(2,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;

% plot sampling points
points_y = D(points_x,16, pi/8 * 0.5);
% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(2,:), 'LineStyle',linestyles{2});
    hold on;
end


Dirichlet_kernel_3 = plot(x,D(x,16, pi/8 * 1),'k','LineWidth',1, 'Color', colors(3,:), 'LineStyle', linestyles{2});
hold on;
x_fill = linspace(-pi/32, pi/32, 500)+ pi/8 * 1;
y_fill = D(x_fill,16,pi/8 * 1);
fill_zone_3 = fill([x_fill, fliplr(x_fill)], [zeros(size(y_fill)), fliplr(y_fill)], colors(3,:), 'FaceAlpha', 0.2,'EdgeAlpha', 0);
hold on;

% plot sampling points
points_y = D(points_x,16, pi/8 * 1);
% plot the vertical line
for i = 1:17
    plot([points_x(i) points_x(i)],[0 points_y(i)],'LineWidth',0.5,'Color',colors(3,:), 'LineStyle',linestyles{2});
    hold on;
end

% plot the points
points_y = D(points_x,16, 0);
Sampling_1 = plot(points_x,points_y,'o','MarkerSize',5,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
hold on;

points_y = D(points_x,16, pi/8 * 0.5);
Sampling_2 = plot(points_x,points_y,'o','MarkerSize',8,'MarkerEdgeColor',colors(2,:));

points_y = D(points_x,16, pi/8 * 1);
Sampling_3 = plot(points_x,points_y,'o','MarkerSize',11,'MarkerEdgeColor',colors(3,:));

legend([Dirichlet_kernel_1, Dirichlet_kernel_2, Dirichlet_kernel_3, Sampling_1, Sampling_2, Sampling_3], ...
    {'Dirichlet kernel at wavenumber-domain $l_x = 0$', 'Dirichlet kernel at wavenumber-domain $l_x = 1$', 'Dirichlet kernel at wavenumber-domain $l_x = 2$', 'Angular-domain sampling result at $n_x^\prime = 0$', 'Angular-domain sampling result at $n_x^\prime = 1$', 'Angular-domain sampling result at $n_x^\prime = 2$'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10,'Position', [0.539763633202325,0.684938747351909,0.361904394384066,0.234375000000001]);

xlim([-0.5 2]);
ylim([-0.3 1.4]);
xlabel('$\gamma = \frac{2\pi}{N_x} (\frac{1}{2}l_x - n_x^\prime)$');
ylabel('$D_{N_x}(\gamma)$');
grid on;
title('(a) Equavalent Dirichlet kernel sampling in angular domain with $\delta = 1/4 \lambda$', 'Interpreter', 'latex'); 

text(-0.377675097276265,1.2,'$(\sigma_0^f)^2 = \int_{\Omega_0 } A^2 $','Color',colors(1,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'arrow',[0.2381,0.2751],[0.8552,0.8484],'LineWidth',0.8, 'Color', colors(1,:));

text(0.016293774319065,1.2,'$(\sigma_1^f)^2 = \int_{\Omega_1 } A^2 $','Color',colors(2,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'arrow',[0.3173,0.3504],[0.8614,0.8523],'LineWidth',0.8, 'Color', colors(2,:));

text(0.439445525291828,1.2,'$(\sigma_2^f)^2 = \int_{\Omega_2 } A^2 $','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'arrow',[0.4574,0.4152],[0.8653,0.851],'LineWidth',0.8, 'Color', colors(3,:));

text(-0.027480544747082,-0.166675860438544,0,'$\Omega_0$','Color',colors(1,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'doublearrow',[0.3143,0.2546],[0.6333,0.6333],'LineWidth',0.8, 'Color', colors(1,:));

text(0.162208171206225,-0.166675860438544,0,'$\Omega_1$','Color',colors(2,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'doublearrow',[0.3761,0.3164],[0.6333,0.6333],'LineWidth',0.8, 'Color', colors(2,:));

text(0.368920233463034,-0.166675860438544,0,'$\Omega_2$','Color',colors(3,:),'Interpreter','latex', 'FontSize',12);
annotation(figure(1),'doublearrow',[0.4387,0.379],[0.6333,0.6333],'LineWidth',0.8, 'Color', colors(3,:));

% set window position : 1000,800.3333333333333,784.3333333333333,437.3333333333333
set(gcf,'Position',[436.3333333333333,406.3333333333333,884.6666666666667,514.6666666666667]);
% save .fig and .pdf to ../../manuscript/figure/AD_compare5.fig and .pdf
saveas(gcf,'../../manuscript/figure/AD_compare5.fig');
saveas(gcf,'../../manuscript/figure/AD_compare5.pdf');

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




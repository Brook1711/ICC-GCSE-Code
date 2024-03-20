clc
clear all;
close all;
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

perfect_color = colors(4, :);
offset_color = colors(11, :);

linewidth = 1;
markersize = 15;
AD_xlabel = '$\sin(\theta)\cos(\phi)$';
AD_ylabel = '$\sin(\theta)\sin(\phi)$';
zlabel_str = '$D(x;N=10)\otimes D(y;N=10)$';
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontName', 'Times New Roman')


% Parameters
N = 10;  % Order of the kernel
x_range = -pi:0.02:pi;  % Range for x values
y_range = -pi:0.02:pi;  % Range for y values

% define sampling points
x_perfect_sam = linspace(-pi,pi,N+1);
y_perfect_sam = linspace(-pi,pi,N+1);

x_offset_sam = linspace(-pi,pi,N+1) + 0.2;
y_offset_sam = linspace(-pi,pi,N+1) + 0.2;

% 2D sampling matrix
[x_perfect_sam,y_perfect_sam] = meshgrid(x_perfect_sam,y_perfect_sam);
Z_sam = D2(x_perfect_sam, y_perfect_sam, N);
[x_offset_sam,y_offset_sam] = meshgrid(x_offset_sam,y_offset_sam);
Z_offset_sam = D2(x_offset_sam, y_offset_sam, N);


% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);


% Evaluate the 2D Dirichlet kernel
Z = D2(X, Y, N);

fig = figure(1);
subplot(1,2,1);
% Plot the surface, no edge line
surf(X/pi, Y/pi, Z, "EdgeAlpha", 0, "FaceAlpha", 0.6);
hold on;  % This allows us to overlay multiple plots
scatter3(x_perfect_sam(:)/pi, y_perfect_sam(:)/pi, Z_sam(:), markersize,perfect_color, 'filled');  % Plot sampled points as red dots

% Plot vertical lines from scatter points to the xy-plane
for i = 1:numel(x_perfect_sam)
    line([x_perfect_sam(i) x_perfect_sam(i)]/pi, [y_perfect_sam(i) y_perfect_sam(i)]/pi, [Z_sam(i) 0], 'Color', perfect_color, 'LineWidth', linewidth, 'LineStyle', '--');
end

xlabel(AD_xlabel);
ylabel(AD_ylabel);
xlim([-1 1]);
ylim([-1 1]);
zlim([-0.3 1]);
view([120, 20]);
zlabel(zlabel_str);
% colorbar("alpha", 0.5);
% title('2D Dirichlet Kernel Surface');

subplot(1,2,2);
% Plot the surface, no edge line
surf(X/pi, Y/pi, Z, "EdgeAlpha", 0, "FaceAlpha", 0.6);
hold on;  % This allows us to overlay multiple plots
scatter3(x_offset_sam(:)/pi, y_offset_sam(:)/pi, Z_offset_sam(:), markersize,  offset_color,  'filled');  % Plot sampled points as blue dots
for i = 1:numel(x_offset_sam)
    line([x_offset_sam(i) x_offset_sam(i)]/pi, [y_offset_sam(i) y_offset_sam(i)]/pi, [Z_offset_sam(i) 0], 'Color', offset_color, 'LineWidth', linewidth, 'LineStyle', '--');
end


xlabel(AD_xlabel);
ylabel(AD_ylabel);
xlim([-1 1]);
ylim([-1 1]);
zlim([-0.3 1]);
view([120, 20]);

zlabel(zlabel_str);
% title('2D Dirichlet Kernel Surface');
h = axes(fig, 'visible', 'off');
c = colorbar(h,'Position', [0.92 0.11 0.02 0.8150]);
colormap(c, 'Parula');

% colorbar;

% Define the 2D Dirichlet kernel function
function y = D2(x, y, N)
    Dx = D1(x, N);
    Dy = D1(y, N);
    y = Dx .* Dy;
end

function y = D1(x, N)
    if abs(x) < eps
        y = 1;
    else
        y = 1/N * sin(N*x/2)./sin(x/2);
    end
    y(isnan(y)) = 1;
end
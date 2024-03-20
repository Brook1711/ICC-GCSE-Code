% clear
close all;
clc
% basic setting
% set global font and fontsize
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',6);
% plot label using latex
set(0,'DefaultTextInterpreter','latex');

% load nosie-free data
load ../data/generate_channel_data.mat
real_WND_2D_mat = channel.WND_2D;
real_AD_2D_mat = channel.AD_2D;

% data loading under low SNR (5dB)
load ../data/SNR_5_Nx_129_RF_1000_spacing_4/alg_GCSE_WND.mat
low_GCSE_WND_2D_mat = abs(vec_H_a_recovered);

load ../data/SNR_5_Nx_129_RF_1000_spacing_4/alg_OMP_WND.mat
low_OMP_WND_2D_mat = abs(vec_H_a_recovered);

load ../data/SNR_5_Nx_129_RF_1000_spacing_4/alg_GCSE_AD.mat
low_GCSE_AD_2D_mat = abs(vec_H_AD_recovered);

load ../data/SNR_5_Nx_129_RF_1000_spacing_4/alg_OMP_AD.mat
low_OMP_AD_2D_mat = abs(vec_H_AD_recovered);

% data loading under high SNR (50dB)
load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_GCSE_WND.mat
high_GCSE_WND_2D_mat = abs(vec_H_a_recovered);

load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_OMP_WND.mat
high_OMP_WND_2D_mat = abs(vec_H_a_recovered);

load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_GCSE_AD.mat
high_GCSE_AD_2D_mat = abs(vec_H_AD_recovered);

load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_OMP_AD.mat
high_OMP_AD_2D_mat = abs(vec_H_AD_recovered);

% plot 2 rows and 5 columns
% 1st row: WND
% 2nd row: AD
% 1st column: noise-free
% 2nd column: high SNR, GCSE
% 3rd column: high SNR, OMP
% 4th column: low SNR, GCSE
% 5th column: low SNR, OMP

% WND xlabel: L_x
% WND ylabel: L_y
WND_xlabel = 'L_x';
WND_ylabel = 'L_y';

% AD xlabel: sin(phi)
% AD ylabel: sin(theta)
AD_xlabel = 'sin(phi)';
AD_ylabel = 'sin(theta)';

% WND colormap: hot
% AD colormap: jet

% set xtick and ytick
% xtick: -(N-1)/2:(N-1)/2
% ytick: -(N-1)/2:(N-1)/2
WND_xtick = -(size(real_WND_2D_mat,2)-1)/2:(size(real_WND_2D_mat,2)-1)/2;
WND_ytick = -(size(real_WND_2D_mat,1)-1)/2:(size(real_WND_2D_mat,1)-1)/2;
real_WND_2D_mat(isnan(real_WND_2D_mat)) = max(max(real_WND_2D_mat));



% WND
figure(1);
% WD from (2,5,1) to (2,5,5)
WND_data_list = {real_WND_2D_mat, high_GCSE_WND_2D_mat, high_OMP_WND_2D_mat, low_GCSE_WND_2D_mat, low_OMP_WND_2D_mat};
WND_title_list = {'WD, noise-free', 'WD, high SNR, GCSE', 'WD, high SNR, OMP', 'WD, low SNR, GCSE', 'WD, low SNR, OMP'};
for i = 1:5
    plot_WND(WND_data_list{i}, WND_title_list{i}, i);
end


% AD
% AD from (2,5,6) to (2,5,10)
AD_data_list = {real_AD_2D_mat, high_GCSE_AD_2D_mat, high_OMP_AD_2D_mat, low_GCSE_AD_2D_mat, low_OMP_AD_2D_mat};
AD_title_list = {'AD, noise-free', 'AD, high SNR, GCSE', 'AD, high SNR, OMP', 'AD, low SNR, GCSE', 'AD, low SNR, OMP'};
for i = 1:5
    plot_AD(AD_data_list{i}, AD_title_list{i}, i+5, 2, 5, 'jet');
end



% new figure, 2cols, 1row
% 1st col: AD, non-trimmed
% 2nd col: AD, trimmed as 1/2 x and 1/2 y
filtered_AD_2D_mat = filter_AD(abs(real_AD_2D_mat), 0.005);

figure(2);
% AD from (1,2,1) to (1,2,2)
AD_data_list = {filtered_AD_2D_mat, real_AD_2D_mat, real_AD_2D_mat};
% trim to reserving the central part, from 1/4 to 3/4
AD_data_list{2} = AD_data_list{2}(size(AD_data_list{2},1)/4+1:3*size(AD_data_list{2},1)/4, size(AD_data_list{2},2)/4+1:3*size(AD_data_list{2},2)/4);
AD_title_list = {'(c) Sparse element indicator, $\delta = \lambda/4$', '(a) ${\bf h}^a$, $\delta = \lambda/2$', '(b) ${\bf h}^a$, $\delta = \lambda/4$'};
subplot_index = [3, 1, 2];

for i = 1:3
    % trim to reserving the central part, from 1/4 to 3/4
    % if i == 1 use hot else use jet
    sub_i = subplot_index(i);
    figure(10+i);
    if i == 1
        plot_AD(AD_data_list{i}, AD_title_list{i}, sub_i, 1, 3, 'hot');
    else
        plot_AD(AD_data_list{i}, AD_title_list{i}, sub_i, 1, 3, 'jet');
        % add colorbar
        
    end
end
% new figure
figure(3);
% plot channel variance in wavenumebr domain
plot_WND(abs(channel.variance), 'Channel variance', 0)

% WD support indicator
figure(4);
WND_mat = channel.variance;
WND_xtick = -(size(WND_mat,2)-1)/2:(size(WND_mat,2)-1)/2;
WND_ytick = -(size(WND_mat,1)-1)/2:(size(WND_mat,1)-1)/2;
% set nan to max
WND_mat(isnan(WND_mat)) = max(max(WND_mat));
WND_xlabel = '$l_x$';
WND_ylabel = '$l_y$';
offset = (size(WND_mat,1) ) / 2;
for i = 1:size(WND_mat,1)
    for j = 1:size(WND_mat,2)
        if sqrt((i - offset)^2 + (j - offset)^2) >= (size(WND_mat,1) ) / 2 -1.5
            WND_mat(i,j) = max(max(WND_mat));
        end
        if WND_mat(i,j) > 0.1
            WND_mat(i,j) = max(max(WND_mat));
        end
    end
end
fig = figure(4);
fig = subplot(1,3,1);
imagesc(abs(WND_mat)...,
    ,'XData',WND_xtick,'YData',WND_ytick);
colormap(fig, hot);
title('Sparse element indicator, $\delta = \lambda / 4$');
xlabel(WND_xlabel);
ylabel(WND_ylabel);
axis equal;
% x, y lim: -32 to 32
xlim([-32, 32]);
ylim([-32, 32]);

% function plot WND, (2,5,i), i from 1 to 5
function plot_WND(WND_mat, title_str, i)
    WND_xtick = -(size(WND_mat,2)-1)/2:(size(WND_mat,2)-1)/2;
    WND_ytick = -(size(WND_mat,1)-1)/2:(size(WND_mat,1)-1)/2;
    % set nan to max
    WND_mat(isnan(WND_mat)) = max(max(WND_mat));
    WND_xlabel = '$l_x$';
    WND_ylabel = '$l_y$';
    if i == 0
        % plot in 3D
        offset = (size(WND_mat,1) ) / 2;
        for i = 1:size(WND_mat,1)
            for j = 1:size(WND_mat,2)
                if sqrt((i - offset)^2 + (j - offset)^2) >= (size(WND_mat,1) ) / 2 -1.5
                    WND_mat(i,j) = nan;
                end
            end
        end
        fig = figure(3);
        surf(WND_mat,'EdgeColor', 'none');
        title(title_str);
        xlabel('$l_x + (|\xi|+1)/2$');
        ylabel('$l_y + (|\xi|+1)/2$');
        xlim([1, size(WND_mat,1)]);
        ylim([1, size(WND_mat,2)]);
        % set surf edge transparent
        % set(fig, 'EdgeColor', 'none');
        % set az and el
        view(25,35);
    else
        fig = subplot(2,5,i);
        imagesc(abs(WND_mat)...,
            ,'XData',WND_xtick,'YData',WND_ytick);
        colormap(fig, hot);
        title(title_str);
        xlabel(WND_xlabel);
        ylabel(WND_ylabel);
        axis equal;
    end

end

% function plot AD, (2,5,i), i from 6 to 10
function plot_AD(AD_mat, title_str, i, row_num, col_num, color_map)
    % flip AD_mat, both x and y
    AD_mat = flip(AD_mat, 1);
    AD_mat = flip(AD_mat, 2);

    AD_xtick = -(size(AD_mat,2)-1)/2:(size(AD_mat,2)-1)/2;
    AD_ytick = -(size(AD_mat,1)-1)/2:(size(AD_mat,1)-1)/2;
    % ticks normalization
    AD_xtick = AD_xtick / max(abs(AD_xtick));
    AD_ytick = AD_ytick / max(abs(AD_ytick));
    % set nan to max
    AD_mat(isnan(AD_mat)) = max(max(AD_mat));
    % xlabel: latex : $\sin(\theta)$
    AD_xlabel = '$\sin(\theta)\cos(\phi)$';
    AD_ylabel = '$\sin(\theta)\sin(\phi)$';
    fig = subplot(row_num,col_num,i);
    imagesc(abs(AD_mat)...,
        ,'XData',AD_xtick,'YData',AD_ytick);
    colormap(fig, color_map);
    title(title_str);
    xlabel(AD_xlabel);
    ylabel(AD_ylabel);
    % set x, y eq
    axis equal;
    % x, y lim: -1 to 1
    xlim([-1, 1]);
    ylim([-1, 1]);
    % set size bigger
    % set(fig, 'Position', [fig.Position(1), fig.Position(2), fig.Position(3)*1.5, fig.Position(4)*1.5]);
    % add colorbar
    % colorbar;
end

% function filter AD_mat, turn the element to 1 if > threshold, 0 if < threshold
function AD_mat = filter_AD(AD_mat, threshold)
    AD_mat(AD_mat > threshold) = 1;
    AD_mat(AD_mat < threshold) = 0;
end
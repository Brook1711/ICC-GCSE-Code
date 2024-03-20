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

title_list = {  '(b) Angular-domain power distribution, ${\bf h}^a$, $\delta = \lambda/4$', ...
                '(a) Angular-domain power distribution, ${\bf h}^a$, $\delta = \lambda/2$,',...
                '(b) Angular-domain sparse element indicator, $\delta = \lambda/4$', ...
                '(a) Wavenumber-domain sparse element indicator, $\delta = \lambda/4$'};


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

% AD xlabel: sin(phi)
% AD ylabel: sin(theta)
AD_xlabel = '$\sin(\theta)\cos(\phi)$';
AD_ylabel = '$\sin(\theta)\sin(\phi)$';


% WND colormap: hot
% AD colormap: jet

% set xtick and ytick
% xtick: -(N-1)/2:(N-1)/2
% ytick: -(N-1)/2:(N-1)/2
WND_mat = channel.variance;
WND_xtick = -(size(WND_mat,2)-1)/2:(size(WND_mat,2)-1)/2;
WND_ytick = -(size(WND_mat,1)-1)/2:(size(WND_mat,1)-1)/2;
% set nan to max
WND_mat(isnan(WND_mat)) = max(max(WND_mat));
WND_xlabel = 'Horizontal wavenumber-domain index, $l_x$';
WND_ylabel = 'Vertical wavenumber-domain index, $l_y$';
offset = (size(WND_mat,1) ) / 2;
for i = 1:size(WND_mat,1)
    for j = 1:size(WND_mat,2)
        if sqrt((i - offset)^2 + (j - offset)^2) >= (size(WND_mat,1) ) / 2 -1.5
            WND_mat(i,j) = max(max(WND_mat));
        end
        if WND_mat(i,j) > 0.1
            WND_mat(i,j) = max(max(WND_mat));
        end
        if WND_mat(i,j) <= 0.1
            WND_mat(i,j) = 0;
        end
    end
end


% new figure, 2cols, 1row
% 1st col: AD, non-trimmed
% 2nd col: AD, trimmed as 1/2 x and 1/2 y
filtered_AD_2D_mat = filter_AD(abs(real_AD_2D_mat), 0.005);

% AD from (1,2,1) to (1,2,2)
AD_data_list = {filtered_AD_2D_mat, real_AD_2D_mat, real_AD_2D_mat};
% trim to reserving the central part, from 1/4 to 3/4
AD_data_list{2} = AD_data_list{2}(size(AD_data_list{2},1)/4+1:3*size(AD_data_list{2},1)/4, size(AD_data_list{2},2)/4+1:3*size(AD_data_list{2},2)/4);

% plot 1 row 2 col according to `AD_data_list(3)` and `AD_data_list(2)`
figure(1);
fig = subplot(1,2,2);
AD_xtick = -(size(AD_data_list{3},2)-1)/2:(size(AD_data_list{3},2)-1)/2;
AD_ytick = -(size(AD_data_list{3},1)-1)/2:(size(AD_data_list{3},1)-1)/2;
AD_xtick = AD_xtick / max(abs(AD_xtick));
AD_ytick = AD_ytick / max(abs(AD_ytick));
imagesc(abs(AD_data_list{3})...,
        ,'XData',AD_xtick,'YData',AD_ytick);colormap(fig, jet);
colormap(fig, jet);
axis equal;
xlim([-1,1]);
ylim([-1,1]);
% set label
xlabel(AD_xlabel);
ylabel(AD_ylabel);
title(title_list{1});
% move tile slightly to the top
title_handle = get(gca,'Title');
title_position = get(title_handle,'Position');
title_position(2) = title_position(2) + 2.5;


fig = subplot(1,2,1);
AD_xtick = -(size(AD_data_list{2},2)-1)/2:(size(AD_data_list{2},2)-1)/2;
AD_ytick = -(size(AD_data_list{2},1)-1)/2:(size(AD_data_list{2},1)-1)/2;
AD_xtick = AD_xtick / max(abs(AD_xtick));
AD_ytick = AD_ytick / max(abs(AD_ytick));
imagesc(abs(AD_data_list{2})...,
        ,'XData',AD_xtick,'YData',AD_ytick);colormap(fig, jet);
axis equal;
xlim([-1,1]);
ylim([-1,1]);
% set label
xlabel(AD_xlabel);
ylabel(AD_ylabel);
title(title_list{2});
% move tile slightly to the top
title_handle = get(gca,'Title');
title_position = get(title_handle,'Position');
title_position(2) = title_position(2) + 2.5;



% plot 1 row 2 col according to `AD_data_list(1)` and `WD_mat`
figure(2);
fig = subplot(1,2,2);
AD_xtick = -(size(AD_data_list{1},2)-1)/2:(size(AD_data_list{1},2)-1)/2;
AD_ytick = -(size(AD_data_list{1},1)-1)/2:(size(AD_data_list{1},1)-1)/2;
AD_xtick = AD_xtick / max(abs(AD_xtick));
AD_ytick = AD_ytick / max(abs(AD_ytick));
imagesc(abs(AD_data_list{1})...,
        ,'XData',AD_xtick,'YData',AD_ytick);colormap(fig, hot);
axis equal;
xlim([-1,1]);
ylim([-1,1]);
% set label
xlabel(AD_xlabel);
ylabel(AD_ylabel);
title(title_list{3});
% move tile slightly to the top
title_handle = get(gca,'Title');
title_position = get(title_handle,'Position');
title_position(2) = title_position(2) + 2.5;


fig = subplot(1,2,1);
% WND_mat flip y
WND_mat = flipud(WND_mat);
% WND_mat flip x 
WND_mat = fliplr(WND_mat);


imagesc(abs(WND_mat)...,
        ,'XData',WND_xtick,'YData',WND_ytick);colormap(fig, hot);
axis equal;
xlim([-32,32]);
ylim([-32,32]);
% set label
xlabel(WND_xlabel);
ylabel(WND_ylabel);
title(title_list{4});
% move tile slightly to the top
title_handle = get(gca,'Title');
title_position = get(title_handle,'Position');
title_position(2) = title_position(2) + 2.5;





% function filter AD_mat, turn the element to 1 if > threshold, 0 if < threshold
function AD_mat = filter_AD(AD_mat, threshold)
    AD_mat(AD_mat > threshold) = 1;
    AD_mat(AD_mat < threshold) = 0;
end
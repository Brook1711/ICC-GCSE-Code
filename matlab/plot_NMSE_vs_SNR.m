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

GCSE_WD_color = colors(4,:);
OMP_WD_color = colors(1,:);
GCSE_AD_color = colors(4,:);
OMP_AD_color = colors(1,:);

U_GCSE_WD_color = colors(7,:);
U_OMP_WD_color = colors(10,:);
U_GCSE_AD_color = colors(7,:);
U_OMP_AD_color = colors(10,:);

GCSE_WD_marker = markers{1};
OMP_WD_marker = markers{2};
GCSE_AD_marker = markers{1};
OMP_AD_marker = markers{2};

GCSE_WD_linestyle = linestyles{1};
OMP_WD_linestyle = linestyles{1};
GCSE_AD_linestyle = linestyles{2};
OMP_AD_linestyle = linestyles{2};

% linewidth
linewidth = 1.5;

% marker size
markersize = 8;

% set global font and fontsize
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',14);
% data load
load('NMSE_vs_SNR.mat');
% 1st col, 2nd col, 3rd col, 4th col, 5th col
% GCSE_WD, OMP_WD,  GCSE_AD, OMP_AD,  SNR

GCSE_WD_data = NMSE_vs_SNR(:,1);
OMP_WD_data = NMSE_vs_SNR(:,2);
GCSE_AD_data = NMSE_vs_SNR(:,3);
OMP_AD_data = NMSE_vs_SNR(:,4);
SNR_list = NMSE_vs_SNR(:,5);

% plot
figure(1);
plot(SNR_list,...
        GCSE_WD_data,...
        'Color',GCSE_WD_color,...
        'Marker',GCSE_WD_marker,...
        'LineStyle',GCSE_WD_linestyle,...
        'LineWidth',linewidth,...
        'MarkerSize',markersize);
hold on;

plot(SNR_list,...
        OMP_WD_data,...
        'Color',OMP_WD_color,...
        'Marker',OMP_WD_marker,...
        'LineStyle',OMP_WD_linestyle,...
        'LineWidth',linewidth,...
        'MarkerSize',markersize);
hold on;

plot(SNR_list,...
        GCSE_AD_data,...
        'Color',GCSE_AD_color,...
        'Marker',GCSE_AD_marker,...
        'LineStyle',GCSE_AD_linestyle,...
        'LineWidth',linewidth,...
        'MarkerSize',markersize);
hold on;

plot(SNR_list,...
        OMP_AD_data,...
        'Color',OMP_AD_color,...
        'Marker',OMP_AD_marker,...
        'LineStyle',OMP_AD_linestyle,...
        'LineWidth',linewidth,...
        'MarkerSize',markersize);

grid on;
xlabel('SNR (dB)');
ylabel('NMSE');
% set log scale
set(gca,'YScale','log');

legend('GCSE, WD','OMP, WD','GCSE, AD','OMP, AD','Location','NorthEast');

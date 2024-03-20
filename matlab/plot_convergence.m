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
markersize = 1;

% set global font and fontsize
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',14);


% data loading
% load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_GCSE_WND.mat
load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_GCSE_WND.mat;
NMSE_GCSE_WND = NMSE_list;
UNMSE_GCSE_WND = NMSE_list_v2;

% load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_OMP_WND.mat
load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_OMP_WND.mat;
NMSE_OMP_WND = NMSE_list;
UNMSE_OMP_WND = NMSE_list_v2;

% load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_GCSE_AD.mat
load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_GCSE_AD.mat;
NMSE_GCSE_AD = NMSE_list;
UNMSE_GCSE_AD = NMSE_list_v2;

% load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_OMP_AD.mat
load ../data/SNR_50_Nx_129_RF_1000_spacing_4/alg_OMP_AD.mat;
NMSE_OMP_AD = NMSE_list;
UNMSE_OMP_AD = NMSE_list_v2;

% start plot
figure(1);
plot(1:length(NMSE_GCSE_WND),...
    NMSE_GCSE_WND,...
    'Color',GCSE_WD_color,...
    'LineStyle',GCSE_WD_linestyle,...
    'LineWidth',linewidth,...
    'MarkerSize',markersize);
hold on;

plot(1:length(NMSE_OMP_WND),...
    NMSE_OMP_WND,...
    'Color',OMP_WD_color,...
    'LineStyle',OMP_WD_linestyle,...
    'LineWidth',linewidth,...
    'MarkerSize',markersize);
hold on;

plot(1:length(NMSE_GCSE_AD),...
    NMSE_GCSE_AD,...
    'Color',GCSE_AD_color,...
    'LineStyle',GCSE_AD_linestyle,...
    'LineWidth',linewidth,...
    'MarkerSize',markersize);
hold on;

plot(1:length(NMSE_OMP_AD),...
    NMSE_OMP_AD,...
    'Color',OMP_AD_color,...
    'LineStyle',OMP_AD_linestyle,...
    'LineWidth',linewidth,...
    'MarkerSize',markersize);
hold on;

% plot(1:length(UNMSE_GCSE_WND),...
%     UNMSE_GCSE_WND,...
%     'Color',U_GCSE_WD_color,...
%     'LineStyle',GCSE_WD_linestyle,...
%     'LineWidth',linewidth,...
%     'MarkerSize',markersize);
% hold on;

% plot(1:length(UNMSE_OMP_WND),...
%     UNMSE_OMP_WND,...
%     'Color',U_OMP_WD_color,...
%     'LineStyle',OMP_WD_linestyle,...
%     'LineWidth',linewidth,...
%     'MarkerSize',markersize);
% hold on;

% plot(1:length(UNMSE_GCSE_AD),...
%     UNMSE_GCSE_AD,...
%     'Color',U_GCSE_AD_color,...
%     'LineStyle',GCSE_AD_linestyle,...
%     'LineWidth',linewidth,...
%     'MarkerSize',markersize);
% hold on;

% plot(1:length(UNMSE_OMP_AD),...
%     UNMSE_OMP_AD,...
%     'Color',U_OMP_AD_color,...
%     'LineStyle',OMP_AD_linestyle,...
%     'LineWidth',linewidth,...
%     'MarkerSize',markersize);
% hold on;

grid on;
% set axis
xlabel('Iteration index');
ylabel('NMSE');
xlim([1,1128]);
ylim([10^-2,10^1.5]);

% set log scale
set(gca,'YScale','log');
legend('GCSE, WD','OMP, WD','GCSE, AD','OMP, AD','Location','NorthEast');
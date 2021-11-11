%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	    Author: Shahab Golshan, 2021 	     %%
%% A sample Matlab code for processing solids holdup %%
%% measurements.				     %%
%%						     %%
%% Before running the code, data from two probe      %%
%% bundles should be saved as experiment1.mat and    %%
%% experiment2.mat in the same folder as the main    %%
%% code.					     %%
%%						     %%
%% Inputs: data (time-series) from two probe bundles %%
%% (experiment1.mat and experiment2.mat)	     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

% Loading experimental data. Here it is assumed that data are saved in .mat format and data from bundle 1 and 2 are stored as experiment1
% and experiment2
load('experiment1.mat');
load('experiment2.mat');

% Segmentation of signals
n_segments = input('Please enter number of data segments ','s');
n_segments = str2num(n_segments);

% Standard deviation elimination threshold
std_elimination_threshold = input('Please enter standard deviation elimination threshold ','s');
std_elimination_threshold = str2num(std_elimination_threshold);

data_number_in_segment = floor(length(voltage_signal_bundle1) / n_segments);
for i = 1 : n_segments
    bundle1_segment(:,i) = voltage_signal_bundle1(((i - 1) * data_number_in_segment + 1) : i * data_number_in_segment);
    bundle2_segment(:,i) = voltage_signal_bundle2(((i - 1) * data_number_in_segment + 1) : i * data_number_in_segment);
end

% Average voltages of bundle1 and bundle2
for i = 1 : n_segments
    average_voltage_bundle1(i) = mean(bundle1_segment(:,i));
    average_voltage_bundle2(i) = mean(bundle2_segment(:,i));
end

% Using calibration equataion to obtain local holdups from bundle1 and bundle2
for i = 1 : n_segments
    local_holdup_bundle1(i) = 0.133 * average_voltage_bundle1(i) + 0.177;
    local_holdup_bundle2(i) = 0.079 * average_voltage_bundle2(i) + 0.457;
end

% Elimination of outliers
local_holdup_bundle1 = local_holdup_bundle1(find(abs(local_holdup_bundle1 - mean(local_holdup_bundle1)) < std_elimination_threshold * std(local_holdup_bundle1)));
local_holdup_bundle2 = local_holdup_bundle2(find(abs(local_holdup_bundle2 - mean(local_holdup_bundle2)) < std_elimination_threshold * std(local_holdup_bundle2)));

% Local holdup is the average of these two values (bundles 1 and 2)
bundle_local_holdup1 = mean(local_holdup_bundle1);
bundle_local_holdup2 = mean(local_holdup_bundle2);
local_holdup = (bundle_local_holdup1 + bundle_local_holdup2) / 2;

disp('The local holdup is: ');
disp(local_holdup);


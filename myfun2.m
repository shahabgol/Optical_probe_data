%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	    Author: Shahab Golshan, 2021 	     %%
%% A sample Matlab code for processing particles'    %%
%% velocity measurements.			     %%
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

% Sampling freuqnecy in experiment
sampling_frequency = input('Please enter sampling frequency ','s');
sampling_frequency = str2num(sampling_frequency);

% Probe equivalent bundle distance
probe_equivalent_length = input('Please enter probe equivalent length ','s');
probe_equivalent_length = str2num(probe_equivalent_length);

% Cross-correlation elimination threshold
cross_correlation_threshold = input('Please enter acceptable cross correlation threshold value ','s');
cross_correlation_threshold = str2num(cross_correlation_threshold);

% Standard deviation elimination threshold
std_elimination_threshold = input('Please enter standard deviation elimination threshold ','s');
std_elimination_threshold = str2num(std_elimination_threshold);

% Number of data segments
n_segments = input('Please enter number of data segments ','s');
n_segments = str2num(n_segments);

% Loading experimental data. Here it is assumed that data are saved in .mat format and data from bundle 1 and 2 are stored as experiment1
% and experiment2
load('experiment1.mat');
load('experiment2.mat');

% Segmentation of signals
data_number_in_segment = floor(length(voltage_signal_bundle1) / n_segments);
for i = 1 : n_segments
    bundle1_segment(:,i) = voltage_signal_bundle1(((i - 1) * data_number_in_segment + 1) : i * data_number_in_segment);
    bundle2_segment(:,i) = voltage_signal_bundle2(((i - 1) * data_number_in_segment + 1) : i * data_number_in_segment);
end

% Obtaining cross correlation curves for each data segment
for i = 1 : n_segments
    [c,lags] = xcorr(bundle2_segment(:,i) , bundle1_segment(:,i));
    maximum_number(i) = find(c == max(c(:)));
    maximum_time(i) = lags(maximum_number(i)) ./ sampling_frequency;
    cross_correlation_value(i) = c(maximum_number(i));

    if maximum_time(i)  == 0
        error('Maximum cross correlation occured at time 0, inappropriate data');
    end
end

% Converting time lag to local average velocity of particles
for i = 1 : n_segments
    raw_particle_velocity(i) = probe_equivalent_length ./ maximum_time(i);
end

% Elimination of inappropriate data because of data segment size
j = 1;
for i = 1 : length(raw_particle_velocity)
    if abs(lags(maximum_number(i))) == length(bundle1_segment(:,i))
        fprintf('data segment size is too short, you can consider decreasing the number of data segments \n');
    else
        raw_particle_velocity_elimination1(j) =  raw_particle_velocity(i);
        j = j + 1;
    end
end

% Elimination of weak cross-correlation values
j = 1;
for i = 1 : length(raw_particle_velocity_elimination1)
    if cross_correlation_value(i) > cross_correlation_threshold
        raw_particle_velocity_elimination2(j) = raw_particle_velocity_elimination1(i);
        j = j + 1;
    end
end

% Elimination of inappropriate data because of sampling frequency
j = 1;
for i = 1 : length(raw_particle_velocity_elimination2)
    if abs(lags(maximum_number(i))) == 1
        fprintf('sampling frequency is not acceptable, you can consider increasing the sampling frequency \n')
    else
        raw_particle_velocity_elimination3(j) = raw_particle_velocity_elimination2(i);
        j = j + 1;
    end
end

% Elimination of outliers
raw_particle_velocity_elimination3 = raw_particle_velocity_elimination3(find(abs(raw_particle_velocity_elimination3 - mean(raw_particle_velocity_elimination3)) < std_elimination_threshold * std(raw_particle_velocity_elimination3)));

% Local particle velocity is the average of different segment values
particle_velocity = mean(raw_particle_velocity_elimination3);

disp('The local velocity is: ');
disp(particle_velocity);
disp(' m/s');

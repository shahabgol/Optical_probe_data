%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	    Author: Shahab Golshan, 2021 	     %%
%% A sample Matlab code for calculation of bubble    %%
%% size and rise velocity.			     %%
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

% Loading experimental data. Here it is assumed that data are saved in .mat
% format and data from bundles 1 and 2 are stored in columns 1 and 2 of the
% experiment.mat file, respectively
mydata = load('experiment.mat');

% Probe equivalent bundle distance
probe_equivalent_length = input('Please enter probe equivalent length  ','s');
probe_equivalent_length = str2num(probe_equivalent_length);

% Minimum acceptable bubble size
minimum_bubble_size = input('Please enter minimum acceptable bubble size ' , 's' );
minimum_bubble_size = str2num(minimum_bubble_size);

% Sampling freuqnecy in experiment
sampling_frequency = input('Please enter sampling frequency ' , 's' );
sampling_frequency = str2num(sampling_frequency);

% Maximum number of data for searching maximum cross-correlation
bubble_velocity_search_domain = input('Please enter maximum number of lags for cross-correlation function ', 's' );
bubble_velocity_search_domain = str2num(bubble_velocity_search_domain);

% Signal threshold for detection of bubbles for bundle1
bubble_signal_threshold1 = input('Please enter bubble signal threshold for bundle1. Signals below this threshold are considered as a bubble ', 's' );
bubble_signal_threshold1 = str2num(bubble_signal_threshold1);

% Signal threshold for detection of bubbles for bundle2
bubble_signal_threshold2 = input('Please enter bubble signal threshold for bundle2. Signals below this threshold are considered as a bubble ', 's' );
bubble_signal_threshold2 = str2num(bubble_signal_threshold2);

% Defining minimum acceptable cross-correlation value for elimination
cross_correlation_threshold = input('Please enter acceptable cross-correlation threshold value ' , 's' );
cross_correlation_threshold = str2num(cross_correlation_threshold);

% Converting minimum acceptable bubble size to data number using sampling frequency
bubble_length_minimum_data_number = minimum_bubble_size * sampling_frequency;

% Defining a few counters for calculations
j = 1;
n = 0;
l = 1;
p = 0;

% Finding bubbles
for i = 1 : length(mydata)
    if mydata(i,1) < bubble_signal_threshold1
        n = n + 1;
    else 
        rawdata1(j) = n;
        n = 0;
        j = j + 1;
    end
    
    if mydata(i,2) < bubble_signal_threshold2
        p = p + 1;
    else 
        rawdata2(l) = p;
        p = 0;
        l = l + 1;
    end
end

% Finding bubbles larger than minimum acceptable bubble size
k = 1;
for i = 1 : length(rawdata1)
    if rawdata1(i) >= bubble_length_minimum_data_number
        bubble1(k) = rawdata1(i);
        k = k + 1;
    end 
end

k = 1;
for i = 1 : length(rawdata2)
    if rawdata2(i) >= bubble_length_minimum_data_number
        bubble2(k) = rawdata2(i);
        k = k + 1;
    end 
end

% Converting to time
t1 = bubble1 ./ sampling_frequency;
t2 = bubble2 ./ sampling_frequency;

% Finding data positions attributed to detected bubbles
accepted_data_locations1 = find(rawdata1 >= bubble_length_minimum_data_number);
accepted_data_locations2 = find(rawdata2 >= bubble_length_minimum_data_number);

if length(t1) < length(t2)
    time = t1;
    accepted_data_locations = accepted_data_locations1;
    signal1 = mydata(:,1);
    signal2 = mydata(:,2);
    rawdata = rawdata1 + 1;
    cumulative_rawdata = cumsum(rawdata);
else 
    time = t2;
    accepted_data_locations = accepted_data_locations2;
    signal1 = mydata(:,1);
    signal2 = mydata(:,2);
    rawdata = rawdata1 + 1;
    cumulative_rawdata = cumsum(rawdata);
end

% Calculation of cross-correlation for bublle rise velocity
for i = 1 : length(accepted_data_locations)
    signalp1 = signal1(cumulative_rawdata(accepted_data_locations(i)-1)-bubble_velocity_search_domain:cumulative_rawdata(accepted_data_locations(i))+bubble_velocity_search_domain);
    signalp2 = signal2(cumulative_rawdata(accepted_data_locations(i)-1)-bubble_velocity_search_domain:cumulative_rawdata(accepted_data_locations(i))+bubble_velocity_search_domain);
    
    
    [c, lags] = xcorr(signalp2 , signalp1);
    maximum_number(i) = find( c == max (c (:)));
    maximum_time(i) = lags(maximum_number( i ) ) ./ sampling_frequency;
    cross_correlation_value(i) = c (maximum_number(i));
end

data_before_elimination = [maximum_time', cross_correlation_value'];

% Elimination of weak cross-correlation values
    j = 1;
    for i = 1 : length(data_before_elimination)
        if data_before_elimination(i,2) > cross_correlation_threshold
            data_after_elimination1(j,:) = data_before_elimination(i,:);
            time_after_elimination1(j) = time(i);
            j = j + 1;
        end
    end
    
    % Elimination of inappropriate data because of sampling frequency
    k = 1;
    for i = 1 : length(data_after_elimination1) 
        if data_after_elimination1(i,1) > 0 
            data_after_elimination2(k,:) = data_after_elimination1(i,:);
            time_after_elimination2(k) = time_after_elimination1(i);
            k = k + 1;
        end
    end
    
maximum_time   = data_after_elimination2(:,1);

% Calculation of bubble rise velocity and bubble size
bubble_velocity  = probe_equivalent_length ./ maximum_time;
bubble_size = bubble_velocity .* time_after_elimination2';

% Visualization and printing information
hist(bubble_size)
xlabel('Bubble size (m)');
ylabel('Number of bubbles');
mean(bubble_size)
std(bubble_size)

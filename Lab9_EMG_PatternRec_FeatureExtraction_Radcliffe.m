%% Lab 9: Pattern Recognition, EMG Feature Extraction - Erin M. Radcliffe
%  BIOE 5073: Neural Interfaces and Bionic Limbs - Spring 2024
%  Instuctors: Weir, R.F.ff; TAs: Huddle, S.A. and Lorelli, S.

% Lab objectives / deliverables:
% Take one posture and extract features from your training data; do not use validation data.
% 1) Graph of Filtered Raw EMG with 'off' sections
% 2) Graph of EMG with quiescent data removed
% 3) Graph of Mean Absolute Value for each bin size (50ms, 100ms, 150ms, 200ms) for 1 EMG sensor in one posture
% 4) Graph of Mean Absolute Value, Slope Sign Change, Zero Crossings, Waveform Length value for 50ms bin size and 2 EMG sensors, one posture

% Notes:
% A) Eventually you will run all data through your code to later do the LDA
%  pattern recognition, so make the loops as flexible as possible. I would
%  suggest using Matlab structures to keep your data organized with names.
% B) To extract 'off' sections automatically you will either need to: 
%    i) set a threshold based on the "off" data or 
%   ii) find a channel that has the best signal to noise ratio within 
%       a posture and then chop out the 'off' sections across all channels 
%       for those 'off' time indices of the one channel.

%% MatLab directory organization

clear all; close all; clc;

main_dir = 'C:\MATLAB\GitHub\EMG_LDA';
raw_data_dir = 'C:\MATLAB\GitHub\EMG_LDA\Group1_2024_quatro_sensors'; % Group 1 EMG
filt_data_dir = 'C:\MATLAB\GitHub\EMG_LDA\Filtered_Group1_2024_quatro_sensors'; % Filtered Group 1 EMG


%% Load & Save each Raw EMG (Group 1) mat file struct (posture) into a single cell array 

cd(raw_data_dir)
mat_files = dir('*.mat');

% each mat file is a 1x1 struct (named posture) with 4 fields:
% Channels (8x25 char), Fs (2.222222167968750e+03),
% Time (8x167761 double), Data (8x167761 double)

posture_names = {'HC', 'HO', 'off', 'RD', 'Tab', 'Tad', 'UD', 'WE', 'WF', 'WP', 'WS'};

n_postures = 11;

% Raw EMG training posture data
for i = 1:n_postures
    load([posture_names{1,i},'_train.mat'])

    emg_train_struct.(posture_names{1,i})= posture.Data;

    clear posture
end
% emg_train_struct: 1x1 struct with 11 fields, each field (posture.Data) is a 8x167761 double

clear i;

n_files = length(mat_files);

allData_cellArray = cell(cell(n_files, 2));

for i = 1:n_files
    % 1st col: name of file struct (posture)
    allData_cellArray{i, 1} = mat_files(i).name;
    % 2nd col: content of file struct (posture)
    allData_cellArray{i, 2} = load(mat_files(i).name);
end

clear i;

%% Filter Design w/ cutoff frequencies

Fs = allData_cellArray{1, 2}.posture.Fs(1); % 2.2222e+03 Hz

% High Pass Filter - 30Hz to remove motion artefact - choose a multiple of 60Hz
[b_hp, a_hp] = butter(4, 30/(Fs/2), 'high'); % 4th order butterworth

% Low Pass Filter - 480Hz to bandlimit to around 500Hz - choose a multiple of 60Hz
[b_lp, a_lp] = butter(4, 480/(Fs/2), 'low'); % 4th order butterworth

% Notch Filter centered at 60Hz
f0 = 60;  % Notch frequency
Q = 35;   % Quality factor, q = ω0/bw where ω0 is the frequency to remove from the signal
[b_notch, a_notch] = iirnotch(f0/(Fs/2), f0/(Fs/2)/Q); % 2nd order infinite impulse response (IIR) notch filter

% Comb Filter centered at 60Hz and its harmonics (120Hz, 180Hz, 240Hz, 300Hz, 360Hz, 420Hz) up to the passband (480 Hz)
f0 = 60;  % Notch frequency
bw = f0 / Q;  % Bandwidth, bw = (fo/(fs/2))/q;
n_harmonics = floor((Fs / 2) / f0);  % Number of harmonics within Nyquist frequency
[b_comb, a_comb] = iircomb(floor(Fs/f0), bw/(Fs/2), 'notch');  % Design comb filter


%% Apply filters sequentially to each posture data file

% Define notch frequencies at 60Hz and harmonics (120Hz, 180Hz, 240Hz, 300Hz, 360Hz, 420Hz) up to the passband (480 Hz)
notch_freqs = 60:60:480;

for i = 1:length(allData_cellArray)
    data_raw = allData_cellArray{i, 2}.posture.Data;  % extract raw data from current file
    Fs = allData_cellArray{i, 2}.posture.Fs(1);  % extract sampling frequency from current file
    data_filtered = zeros(size(data_raw));  % initialize filtered data array

    for k = 1:size(data_raw, 1)  % iterate over channels
        % High-pass filter
        data_filtered(k, :) = filtfilt(b_hp, a_hp, data_raw(k, :));
        % Low-pass filter
        data_filtered(k, :) = filtfilt(b_lp, a_lp, data_filtered(k, :));

        % Notch filters at harmonics
        for f = notch_freqs
            [b_notch, a_notch] = iirnotch(f/(Fs/2), f/(Fs/2)/Q);
            data_filtered(k, :) = filtfilt(b_notch, a_notch, data_filtered(k, :));
        end
    end

    % Store filtered data back in allData_cellArray
    allData_cellArray{i, 2}.posture.DataFiltered = data_filtered;
end

clear i;

%% save allData_cellArray in filt_data_dir

% Change to the filtered data directory
cd(filt_data_dir);

% Iterate over all posture files in allData_cellArray and save them
for i = 1:length(allData_cellArray)
    % Extract file name without the '.mat' extension
    [~, filename, ~] = fileparts(allData_cellArray{i, 1});

    % Create a new file name for the filtered data
    new_filename = [filename '_Filtered.mat'];

    % Extract the posture data
    posture_data = allData_cellArray{i, 2}.posture;
    
    % Add the filtered data to the posture struct
    posture_data.DataFiltered = allData_cellArray{i, 2}.posture.DataFiltered;

    % Save the updated posture struct to a new file in the filtered data directory
    save(new_filename, 'posture_data');
end


%% Plot filtered EMG data along time axis for each channel of HO_train

% Select the HO_train posture data
ho_train_idx = find(strcmp({allData_cellArray{:,1}}, 'HO_train.mat'));
ho_train_data = allData_cellArray{ho_train_idx, 2}.posture;

figure;
for j = 1:size(ho_train_data.DataFiltered, 1)
    subplot(4, 2, j);
    plot(ho_train_data.Time(j,:), ho_train_data.DataFiltered(j,:));
    title(['Filtered EMG Data, Channel ' num2str(j)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
end
sgtitle('HO_train: Filtered EMG Along Time Axis', 'Interpreter', 'none');


%% Determine method for determining "OFF" sections - notes

% 11 postures, 8 EMGs, 4 features per EMG 
% 6 flexions ~ 60sec, ~30sec after OFF sections are removed

% mean and std dev (or variance) of row - mean per posture 
% mean and std dev (or variance) of col - mean per channel 

% root mean squared of absolute value

% aggressively smooth the filtered EMG - smoothdata or low pass to extreme degree 
% determine indices of rise/fall threshold using median, mean, standard deviation

% determine indices of ON/OFF threshold crossing for strong posture signal 
% superimpose indices on all other postures 

% Remove "OFF" sections, concatanate "On" sections
% plot original signal compared to concatanated signal (with Off sections removed)

% Error-prone but valid method:
% For all 6 postures: 10 sec lead in, flex ON for 5 sec, OFF for 5 sec, 10 sec lead out


%% 1) Graph Filtered EMG training data with 'off' sections for 'HO_train.mat' - sensor 1

% Fs = allData_cellArray{1, 2}.posture.Fs(1); % 2.2222e+03 Hz

upThresh = mean(ho_train_data.DataFiltered(1,:))+(1.96*std(ho_train_data.DataFiltered(1,:))); % Threshold based on deviation, 95th based
downThresh = -upThresh;

% Empty array to store 'on' data
on_data = [];

% 1.5 ms loop spacing 
lower_lim_eval_window = 1;                                          
upper_lim_eval_window = 333;                                        
window_increment = 333;                                      
loop_limit = floor(length(ho_train_data.DataFiltered(1,:))/window_increment);  

% Loop to determine and store 'on' data
for k = 1:loop_limit-1
    % Extract the segment of the data
    segment = ho_train_data.DataFiltered(1, lower_lim_eval_window:upper_lim_eval_window);
    
    % Find the maximum and minimum values in the segment
    checkPtmax = max(segment);     % Looks for window max
    checkPtmin = min(segment);     % Looks for window min
    
    % Check if the max or min exceed the threshold
    if checkPtmax > upThresh || checkPtmin < downThresh
        on_data = [on_data, segment];
        % disp("Yes");
    else 
        % disp("No");
    end

    % Increment the window
    lower_lim_eval_window = 1 + window_increment*k;                   
    upper_lim_eval_window = window_increment + window_increment*k;
end

% Generate the time vector for the 'on' data
on_data_Time = ho_train_data.Time(1,1:length(on_data));

% Plot the original and the 'on' data
figure(50);
subplot(2,1,1);
plot(ho_train_data.Time, ho_train_data.DataFiltered(1,:), 'b-');
title('Original Data');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-4e-4 4e-4]);
xlim([0 ho_train_data.Time(end)]);

subplot(2,1,2);
plot(on_data_Time, on_data, 'r-');
title('Data with "Off" Sections Removed');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-4e-4 4e-4]);
xlim([0 ho_train_data.Time(end)]);

% Adjust the plot settings
sgtitle('HO Train, Sensor 1');


%% 2) Graph Filtered EMG training data with quiescent data removed for 'HO_train.mat' - sensor 1

% Define quiescent threshold
quietThresh = mean(on_data(1,:)) + (0.5*std(on_data(1,:))); % Adjust multiplier as needed

% Initialize array to store 'active' data, where the muscle is not quiescent
active_data = [];
active_time = [];

% Define the bin size for iteration, could be the same as used for 'off' detection
binSize = 1.5e-3; % 1.5 ms in seconds
binSamples = round(binSize * Fs); % Convert bin size to number of samples

% Loop to identify and store 'active' data
for i = 1:binSamples:length(on_data(1,:)) - binSamples
    segment = on_data(1, i:i+binSamples-1);
    % Check if the segment's activity exceeds the quiescent threshold
    if max(abs(segment)) > quietThresh
        active_data = [active_data, segment];
        active_time = [active_time, on_data_Time(1, i:i+binSamples-1)];
    end
end

% Plot the 'active' data
figure;
plot(active_time, active_data, 'r-');
title('HO train, Sensor 1: EMG Data with Quiescent Data Removed');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-4e-4 4e-4]);
xlim([0, on_data_Time(1,end)]);

%% 3) Graph Mean Absolute Value for each bin size (50ms, 100ms, 150ms, 200ms) 
% for 1 EMG sensor in one posture ('HO_train.mat') 

% Define the bin sizes
binSizes = [0.050, 0.1, 0.15, 0.2];  % in seconds
binSamples = round(binSizes * Fs);  % Convert bin sizes to number of samples

% Initialize arrays to store MAV values for each bin size
MAV_values = cell(length(binSizes), 1);

% Loop over each bin size and calculate MAVs
for i = 1:length(binSizes)
    numBins = floor(length(active_data) / binSamples(i));
    MAVs = zeros(1, numBins);
    
    for bin = 1:numBins
        startIdx = (bin - 1) * binSamples(i) + 1;
        endIdx = bin * binSamples(i);
        MAVs(bin) = mean(abs(active_data(startIdx:endIdx)));
    end
    
    MAV_values{i} = MAVs; % MAV_values cell array
end


% Plotting MAV for each bin size
figure;
titles = {'50ms bin size', '100ms bin size', '150ms bin size', '200ms bin size'};
for i = 1:length(binSizes)
    subplot(4, 1, i);
    stem((1:length(MAV_values{i})) * binSizes(i), MAV_values{i}, 'Marker', 'none');
    title(titles{i});
    xlabel('Time (s)');
    ylabel('MAV (V)');
end
sgtitle('Mean Absolute Value for Different Bin Sizes - HO train, Sensor 1');


%% %% %% future step: make nested loop for all EMG sensors
% for now, see clunky/repetative but functional process below

%% 1) Graph Filtered EMG training data with 'off' sections for 'HO_train.mat' - sensor 2

% Fs = allData_cellArray{1, 2}.posture.Fs(1); % 2.2222e+03 Hz

upThresh = mean(ho_train_data.DataFiltered(2,:))+(1.96*std(ho_train_data.DataFiltered(2,:))); % Threshold based on deviation, 95th based
downThresh = -upThresh;

% Empty array to store 'on' data
on_data_2 = [];

% 1.5 ms loop spacing 
lower_lim_eval_window = 1;                                          
upper_lim_eval_window = 333;                                        
window_increment = 333;                                      
loop_limit = floor(length(ho_train_data.DataFiltered(2,:))/window_increment);  

% Loop to determine and store 'on' data
for k = 1:loop_limit-1
    % Extract the segment of the data
    segment = ho_train_data.DataFiltered(2, lower_lim_eval_window:upper_lim_eval_window);
    
    % Find the maximum and minimum values in the segment
    checkPtmax = max(segment);     % Looks for window max
    checkPtmin = min(segment);     % Looks for window min
    
    % Check if the max or min exceed the threshold
    if checkPtmax > upThresh || checkPtmin < downThresh
        on_data_2 = [on_data_2, segment];
        % disp("Yes");
    else 
        % disp("No");
    end

    % Increment the window
    lower_lim_eval_window = 1 + window_increment*k;                   
    upper_lim_eval_window = window_increment + window_increment*k;
end

% Generate the time vector for the 'on' data
on_data_2_Time = ho_train_data.Time(2,1:length(on_data_2));

% Plot the original and the 'on' data
figure(50);
subplot(2,1,1);
plot(ho_train_data.Time, ho_train_data.DataFiltered(2,:), 'b-');
title('Original Data');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-4e-4 4e-4]);
xlim([0 ho_train_data.Time(end)]);

subplot(2,1,2);
plot(on_data_2_Time, on_data_2, 'r-');
title('Data with "Off" Sections Removed');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-4e-4 4e-4]);
xlim([0 ho_train_data.Time(end)]);

% Adjust the plot settings
sgtitle('HO train, Sensor 2');


%% 2) Graph Filtered EMG training data with quiescent data removed for 'HO_train.mat' - sensor 2

% Define quiescent threshold
quietThresh_2 = mean(on_data_2(1,:)) + (0.5*std(on_data_2(1,:))); % Adjust multiplier as needed

% Initialize array to store 'active' data, where the muscle is not quiescent
active_data_2 = [];
active_time_2 = [];

% Define the bin size for iteration, could be the same as used for 'off' detection
binSize = 1.5e-3; % 1.5 ms in seconds
binSamples = round(binSize * Fs); % Convert bin size to number of samples

% Loop to identify and store 'active' data
for i = 1:binSamples:length(on_data_2(1,:)) - binSamples
    segment = on_data_2(1, i:i+binSamples-1);
    % Check if the segment's activity exceeds the quiescent threshold
    if max(abs(segment)) > quietThresh_2
        active_data_2 = [active_data_2, segment];
        active_time_2 = [active_time_2, on_data_2_Time(1, i:i+binSamples-1)];
    end
end

% Plot the 'active' data
figure;
plot(active_time_2, active_data_2, 'r-');
title('HO train, Sensor 2: EMG Data with Quiescent Data Removed');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-4e-4 4e-4]);
xlim([0, on_data_2_Time(1,end)]);


%% 3) Graph Mean Absolute Value for each bin size (50ms, 100ms, 150ms, 200ms) for 'HO_train.mat' - sensor 2
% for 1 EMG sensor in one posture ('HO_train.mat') 

% Define the bin sizes
binSizes = [0.050, 0.1, 0.15, 0.2];  % in seconds
binSamples = round(binSizes * Fs);  % Convert bin sizes to number of samples

% Initialize arrays to store MAV values for each bin size
MAV_values = cell(length(binSizes), 1);

% Loop over each bin size and calculate MAVs
for i = 1:length(binSizes)
    numBins = floor(length(active_data_2) / binSamples(i));
    MAVs = zeros(1, numBins);
    
    for bin = 1:numBins
        startIdx = (bin - 1) * binSamples(i) + 1;
        endIdx = bin * binSamples(i);
        MAVs(bin) = mean(abs(active_data_2(startIdx:endIdx)));
    end
    
    MAV_values{i} = MAVs; % MAV_values cell array
end

% Plotting MAV for each bin size
figure;
titles = {'50ms bin size', '100ms bin size', '150ms bin size', '200ms bin size'};
for i = 1:length(binSizes)
    subplot(4, 1, i);
    stem((1:length(MAV_values{i})) * binSizes(i), MAV_values{i}, 'Marker', 'none');
    title(titles{i});
    xlabel('Time (s)');
    ylabel('MAV (V)');
end
sgtitle('Mean Absolute Value for Different Bin Sizes - HO train, Sensor 2');


%% 4) Graph Mean Absolute Value, Slope Sign Change, Zero Crossings, Waveform Length value for the 50ms bin size and 2 EMG sensors in one posture ('HO_train.mat')

% Prepare the data for feature extraction

% Set the bin size to 50ms
binSize = 0.050; % 50 ms in seconds
binSamples = round(binSize * Fs); % Convert bin size to number of samples per bin

% Initialize feature matrices
MAV = cell(1, 2);
SSC = cell(1, 2);
ZC = cell(1, 2);
WL = cell(1, 2);

% Define active data for HO train, sensor 1 and sensor 2
active_data_s1s2 = {active_data, active_data_2};

% Loop over the 'active' data in bins to calculate features for each sensor
for sensor_idx = 1:2
    active_data_sensor = active_data_s1s2{sensor_idx};
    
    % Determine the number of full bins we can calculate
    numBins = floor(length(active_data_sensor) / binSamples);
    
    % Initialize the feature arrays for each sensor based on the number of bins
    MAV{sensor_idx} = zeros(1, numBins);
    SSC{sensor_idx} = zeros(1, numBins);
    ZC{sensor_idx} = zeros(1, numBins);
    WL{sensor_idx} = zeros(1, numBins);
    
    for bin = 1:numBins
        startIdx = (bin - 1) * binSamples + 1;
        endIdx = bin * binSamples;
        binData = active_data_sensor(startIdx:endIdx);
        
        % Calculate features for each bin
        MAV{sensor_idx}(bin) = mean(abs(binData));
        SSC{sensor_idx}(bin) = sum(diff(sign(diff(binData)))~=0);
        ZC{sensor_idx}(bin) = sum(diff(sign(binData))~=0);
        WL{sensor_idx}(bin) = sum(abs(diff(binData)));
    end
end

% Calculate feature_time based on the actual number of bins calculated
feature_time = ((1:numBins) - 0.5) * binSize;  % Center time in each bin


% Verify the number of bins and adjust feature_time accordingly - test for single feature (MAV)
numBins = min([length(MAV{1}), length(MAV{2})]);  % Use the smaller number of bins if different
feature_time = ((1:numBins) - 0.5) * binSize;

% Trim the feature vectors if necessary
MAV{1} = MAV{1}(1:numBins);
MAV{2} = MAV{2}(1:numBins);

% Debugging prints
disp(['Length of feature_time: ', num2str(length(feature_time))]);
disp(['Length of MAV{1}: ', num2str(length(MAV{1}))]);
disp(['Length of MAV{2}: ', num2str(length(MAV{2}))]);

% Now plot the features with the corrected lengths
figure;
title('Feature Extraction for HO_train');

% Plot Mean Absolute Value for Sensor 1 and Sensor 2
hold on;
plot(feature_time, MAV{1}, 'r');
plot(feature_time, MAV{2}, 'b');
title('Mean Absolute Value (MAV)');
xlabel('Time (s)');
ylabel('MAV');
legend('Sensor 1', 'Sensor 2');
hold off;


%%  Verify the number of bins and adjust feature_time accordingly for all features

% Determine the smallest number of bins across all sensors and features
numBins = min([
    length(MAV{1}), length(MAV{2}),
    length(SSC{1}), length(SSC{2}),
    length(ZC{1}), length(ZC{2}),
    length(WL{1}), length(WL{2})
]);

% Adjust feature vectors to have the same number of elements
MAV{1} = MAV{1}(1:numBins);
MAV{2} = MAV{2}(1:numBins);
SSC{1} = SSC{1}(1:numBins);
SSC{2} = SSC{2}(1:numBins);
ZC{1} = ZC{1}(1:numBins);
ZC{2} = ZC{2}(1:numBins);
WL{1} = WL{1}(1:numBins);
WL{2} = WL{2}(1:numBins);

% Adjust feature_time
feature_time = ((1:numBins) - 0.5) * binSize;

%% Plot extracted features 

figure;
t = tiledlayout(4,1);
title(t,'Feature Extraction for HO_train');

% Plot Mean Absolute Value
nexttile;
hold on;
plot(feature_time, MAV{1}, 'r', 'DisplayName', 'Sensor 1 MAV');
plot(feature_time, MAV{2}, 'b', 'DisplayName', 'Sensor 2 MAV');
title('Mean Absolute Value (MAV)');
xlabel('Time (s)');
ylabel('MAV');
legend;
hold off;


% Plot Slope Sign Changes
nexttile;
hold on;
plot(feature_time, SSC{1}, 'r', 'DisplayName', 'Sensor 1 SSC');
plot(feature_time, SSC{2}, 'b', 'DisplayName', 'Sensor 2 SSC');
title('Slope Sign Changes (SSC)');
xlabel('Time (s)');
ylabel('SSC');
legend;
hold off;

% Plot Zero Crossings
nexttile;
hold on;
plot(feature_time, ZC{1}, 'r', 'DisplayName', 'Sensor 1 ZC');
plot(feature_time, ZC{2}, 'b', 'DisplayName', 'Sensor 2 ZC');
title('Zero Crossings (ZC)');
xlabel('Time (s)');
ylabel('ZC');
legend;
hold off;

% Plot Waveform Length
nexttile;
hold on;
plot(feature_time, WL{1}, 'r', 'DisplayName', 'Sensor 1 WL');
plot(feature_time, WL{2}, 'b', 'DisplayName', 'Sensor 2 WL');
title('Waveform Length (WL)');
xlabel('Time (s)');
ylabel('WL');
legend;
hold off;



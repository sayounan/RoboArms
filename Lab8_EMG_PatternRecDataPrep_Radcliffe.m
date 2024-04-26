%% Lab 8: Pattern Recognition Lab: Data Preparation - Erin M. Radcliffe
%  BIOE 5073: Neural Interfaces and Bionic Limbs - Spring 2024
%  Instuctors: Weir, R.F.ff; TAs: Huddle, S.A. and Lorelli, S.

% Lab objectives / deliverables:
% 1) Graph of Raw EMG
% 2) Filter Design plot w/ cutoff frequencies
%    a) High Pass - 30Hz to remove motion artefact - choose a multiple of 60Hz
%    b) Low Pass - 480Hz to bandlimit to around 500Hz - choose a multiple of 60Hz
%    c) Notch at 60Hz and its harmonics 120Hz, 180Hz, 240Hz, 300Hz, 360Hz, 420Hz up to the Passband
% 3) Frequency vs. Power plot after all filters have been implemented
% 4) Filtered EMG along time axis


%% MatLab directory organization

clear all; close all; clc;

main_dir = 'C:\MATLAB\GitHub\EMG_LDA';
raw_data_dir = 'C:\MATLAB\GitHub\EMG_LDA\Group1_2024_quatro_sensors'; % Group 1 EMG
filt_data_dir = 'C:\MATLAB\GitHub\EMG_LDA\Filtered_Group1_2024_quatro_sensors'; % Filtered Group 1 EMG

cd(raw_data_dir)


%% Data exploration: load Raw EMG (Group 1) and save training and validation posture data into separate structs

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

% Raw EMG validation posture data
for i = 1:n_postures
    load([posture_names{1,i},'_val.mat'])

    emg_val_struct.(posture_names{1,i})= posture.Data;

    clear posture
end
% emg_val_struct: 1x1 struct with 11 fields, each field (posture.Data) is a 8x167761 double

% EMG training data and validation data are equivalent in dimension -->
% save all data in a single storage container


%% Save each Raw EMG mat file struct (posture) in a single cell array

n_files = length(mat_files);

allData_cellArray = cell(cell(n_files, 2));

for i = 1:n_files
    % 1st col: name of file struct (posture)
    allData_cellArray{i, 1} = mat_files(i).name;
    % 2nd col: content of file struct (posture)
    allData_cellArray{i, 2} = load(mat_files(i).name);
end

clear i;


%% Plot all channels of Raw EMG per posture

n_channels = height(allData_cellArray{1,2}.posture.Channels);

for i = 1:length(allData_cellArray)
    figure(i)
    for j = 1:n_channels
        subplot(4, 2, j);
        plot(allData_cellArray{i,2}.posture.Time(j,:), allData_cellArray{i,2}.posture.Data(j,:));
        title("Raw Data, Channel "+ j);
        sgtitle(allData_cellArray{i,1}, 'Interpreter', 'none');
        ylim([-2e-3, 2e-3]);
    end
end


%% Filter Design plotting w/ cutoff frequencies

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

% Define x-axis ticks at every 60 Hz
x_ticks = 0:60:Fs/2;

% Plot frequency response of high-pass filter
[h_hp, w_hp] = freqz(b_hp, a_hp);
f_hp = w_hp * Fs / (2 * pi); % Convert from rad/sample to Hz
mag_hp = 20 * log10(abs(h_hp)); % Convert magnitude to dB
figure;
plot(f_hp, mag_hp);
title('High-pass Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-100, 10]);
xticks(x_ticks);

% Plot frequency response of low-pass filter
[h_lp, w_lp] = freqz(b_lp, a_lp);
f_lp = w_lp * Fs / (2 * pi); % Convert from rad/sample to Hz
mag_lp = 20 * log10(abs(h_lp)); % Convert magnitude to dB
figure;
plot(f_lp, mag_lp);
title('Low-pass Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-200, 10]);
xticks(x_ticks);

% Plot frequency response of notch filter
[h_notch, w_notch] = freqz(b_notch, a_notch);
f_notch = w_notch * Fs / (2 * pi); % Convert from rad/sample to Hz
mag_notch = 20 * log10(abs(h_notch)); % Convert magnitude to dB
figure;
plot(f_notch, mag_notch);
title('Notch Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-5, 1]);
xticks(x_ticks);

% Plot frequency response of comb filter
[h_comb, w_comb] = freqz(b_comb, a_comb);
f_comb = w_comb * Fs / (2 * pi);  % Convert from rad/sample to Hz
mag_comb = 20 * log10(abs(h_comb));  % Convert magnitude to dB
figure;
plot(f_comb, mag_comb);
title('Comb Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-25, 1]);
xticks(x_ticks);


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


%% save allData_cellArray with added DataFiltered field per posture struct in filt_data_dir

% Change to the filtered data directory
cd(filt_data_dir);

% Iterate over all posture files in allData_cellArray and save each modified structure back into a .mat file
for i = 1:length(allData_cellArray)
    % Extract file name without the '.mat' extension
    [~, filename, ~] = fileparts(allData_cellArray{i, 1});

    % Create a new file name for the filtered data
    new_filename = [filename '_Filtered.mat'];

    % Extract the posture data
    posture_data = allData_cellArray{i, 2}.posture;
    
    % Add the filtered data to the posture struct
    posture_data.DataFiltered = allData_cellArray{i, 2}.posture.DataFiltered;

    % Save the updated posture struct to a new .mat file in the filtered data directory
    save(new_filename, 'posture_data');
end

% Change back to the original directory if needed
cd(raw_data_dir);


%% Plot raw and filtered data for comparison (for a single posture)

% Select the HO_train posture data
ho_train_idx = find(strcmp({allData_cellArray{:,1}}, 'HO_train.mat'));
ho_train_data = allData_cellArray{ho_train_idx, 2}.posture;

% Plot raw and filtered data for comparison
figure;
for j = 1:size(ho_train_data.Data, 1)
    subplot(size(ho_train_data.Data, 1), 1, j);
    plot(ho_train_data.Time(j,:), ho_train_data.Data(j,:), 'b');  % Raw data in blue
    hold on;
    plot(ho_train_data.Time(j,:), ho_train_data.DataFiltered(j,:), 'r');  % Filtered data in red
    title(['Channel ' num2str(j)]);
    legend('Raw', 'Filtered');
end
sgtitle('HO_train: Raw vs Filtered EMG Data', 'Interpreter', 'none');


%% Frequency vs. Power plot after all filters have been implemented (for a single posture)

% Fourier transform on filtered data for HO_train
figure;
for h = 1:size(ho_train_data.DataFiltered, 1)
    data_fft = fft(ho_train_data.DataFiltered(h, :));
    f = ho_train_data.Fs * (0:(length(data_fft)/2)) / length(data_fft);
    P2 = abs(data_fft / length(data_fft));
    P1 = P2(1:length(data_fft)/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    % Plot frequency vs. power for each channel in a single figure
    subplot(4, 2, h);
    plot(f, P1);
    title(['Channel ' num2str(h)]);
    xlabel('Frequency (Hz)');
    ylabel('Power');
end
sgtitle('HO_train: Frequency vs. Power', 'Interpreter', 'none');

% Visualize frequency vs. power using pwelch for HO_train
figure;
for j = 1:size(ho_train_data.DataFiltered, 1)
    [P, F] = pwelch(ho_train_data.DataFiltered(j, :), ones(8192, 1), 8192/2, 8192, ho_train_data.Fs(1));
    subplot(4, 2, j);
    plot(F, 10*log10(P));  % Convert power to dB
    title(['Channel ' num2str(j)]);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
end
sgtitle('HO_train: Power Spectral Density (pwelch)', 'Interpreter', 'none');


%% Plot filtered EMG data along time axis for each channel (for a single posture)

% Plot filtered EMG data along time axis for each channel for HO_train
figure;
for j = 1:size(ho_train_data.DataFiltered, 1)
    subplot(4, 2, j);
    plot(ho_train_data.Time(j,:), ho_train_data.DataFiltered(j,:));
    title(['Filtered EMG Data, Channel ' num2str(j)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
end
sgtitle('HO_train: Filtered EMG Along Time Axis', 'Interpreter', 'none');



%% Supplemental (filters + plotting code applied to all postures)

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

    % Plot raw and filtered data for comparison
    figure;
    for j = 1:size(data_raw, 1)
        subplot(size(data_raw, 1), 1, j);
        plot(data_raw(j, :), 'b');  % Raw data in blue
        hold on;
        plot(data_filtered(j, :), 'r');  % Filtered data in red
        title(['Channel ' num2str(j)]);
        legend('Raw', 'Filtered');
    end
    sgtitle(['Raw vs Filtered EMG Data - ' allData_cellArray{i, 1}], 'Interpreter', 'none');


    % Frequency vs. Power plot after all filters have been implemented
    % Fourier transform on filtered data
    figure;
    for h = 1:size(data_filtered, 1)
        data_fft = fft(data_filtered(h, :));
        f = Fs*(0:(length(data_fft)/2))/length(data_fft);
        P2 = abs(data_fft/length(data_fft));
        P1 = P2(1:length(data_fft)/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % Plot frequency vs. power for each channel in a single figure
        subplot(4, 2, h);
        plot(f, P1);
        title(['Channel ' num2str(h)]);
        xlabel('Frequency (Hz)');
        ylabel('Power');
    end
    sgtitle(['Frequency vs. Power - ' allData_cellArray{i, 1}], 'Interpreter', 'none');

    % Visualize frequency vs. power using pwelch for each channel per posture
    figure;
    for j = 1:size(data_filtered, 1)
        [P, F] = pwelch(data_filtered(j, :), ones(8192, 1), 8192/2, 8192, Fs(1));
        subplot(4, 2, j);
        plot(F, 10*log10(P));  % Convert power to dB
        title(['Channel ' num2str(j)]);
        xlabel('Frequency (Hz)');
        ylabel('Power/Frequency (dB/Hz)');
    end
    sgtitle([allData_cellArray{i, 1} ': Power Spectral Density (pwelch)'], 'Interpreter', 'none');


    % Plot filtered EMG data along time axis for each channel per posture
    figure;
    for j = 1:size(data_filtered, 1)
        subplot(4, 2, j);
        plot(allData_cellArray{i,2}.posture.Time(j,:), data_filtered(j,:));
        title(['Filtered EMG Data, Channel ' num2str(j)]);
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    sgtitle(['Filtered EMG along Time axis - ' allData_cellArray{i, 1}], 'Interpreter', 'none');

end


%% Setup of the System
filt = load("LPF.mat");
data = load("my_experiment_data.mat");
LPF = filt.LPF;

fs = data.fs;
samples = data.t;
sample_size = length(samples);

% Import temperature data from every sensor.
temp_air1_ = data.temp_ch3; % Air 1
temp_air2_ = data.temp_ch1; % Air 2
temp_surface1_ = data.temp_ch2; % Surface 1
temp_surface2_ = data.temp_ch0; % Surface 2

%% Time Domain Plots

figure("Name","Time Domain Representation of the Signals")
% Top: Surface 1 & Surface 2
subplot(2,2,1); plot(samples, temp_surface1_, "LineWidth", 1.5); title("Surface 1"); xlabel("Sample"); ylabel("Temperature (C)");
subplot(2,2,2); plot(samples, temp_surface2_, "LineWidth", 1.5); title("Surface 2"); xlabel("Sample"); ylabel("Temperature (C)");

% Bottom: Air 1 & Air 2
subplot(2,2,3); plot(samples, temp_air1_, "LineWidth", 1.5); title("Air 1"); xlabel("Sample"); ylabel("Temperature (C)");
subplot(2,2,4); plot(samples, temp_air2_, "LineWidth", 1.5); title("Air 2"); xlabel("Sample"); ylabel("Temperature (C)");

%% Frequency Domain Plots

% Apply fft
TEMP_AIR1_ = fft(temp_air1_, sample_size);
TEMP_AIR2_ = fft(temp_air2_, sample_size);
TEMP_SURFACE1_ = fft(temp_surface1_, sample_size);
TEMP_SURFACE2_ = fft(temp_surface2_, sample_size);

f = linspace(-fs/2, fs/2, sample_size);

figure("Name","Frequency Domain Representation of the Signals")
% Top: Surface 1 & Surface 2
subplot(2,2,1); stem(f, abs(fftshift(TEMP_SURFACE1_/sample_size))); title("Surface 1 Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);
subplot(2,2,2); stem(f, abs(fftshift(TEMP_SURFACE2_/sample_size))); title("Surface 2 Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);

% Bottom: Air 1 & Air 2
subplot(2,2,3); stem(f, abs(fftshift(TEMP_AIR1_/sample_size))); title("Air 1 Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);
subplot(2,2,4); stem(f, abs(fftshift(TEMP_AIR2_/sample_size))); title("Air 2 Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);

%% Apply Filter

temp_air1_filtered = filtfilt(LPF, 1, temp_air1_);
temp_air2_filtered = filtfilt(LPF, 1, temp_air2_);
temp_surface1_filtered = filtfilt(LPF, 1, temp_surface1_);
temp_surface2_filtered = filtfilt(LPF, 1, temp_surface2_);

figure("Name","Filtered Time Domain Representation of the Signals")
% Top: Surface 1 & Surface 2
subplot(2,2,1); plot(samples, temp_surface1_filtered, "LineWidth", 1.5); title("Surface 1 Filtered");
subplot(2,2,2); plot(samples, temp_surface2_filtered, "LineWidth", 1.5); title("Surface 2 Filtered");

% Bottom: Air 1 & Air 2
subplot(2,2,3); plot(samples, temp_air1_filtered, "LineWidth", 1.5); title("Air 1 Filtered");
subplot(2,2,4); plot(samples, temp_air2_filtered, "LineWidth", 1.5); title("Air 2 Filtered");

%% Frequency Domain of Filtered Signals

% Apply fft to filtered signals
TEMP_AIR1_FILTERED_ = fft(temp_air1_filtered, sample_size);
TEMP_AIR2_FILTERED_ = fft(temp_air2_filtered, sample_size);
TEMP_SURFACE1_FILTERED_ = fft(temp_surface1_filtered, sample_size);
TEMP_SURFACE2_FILTERED_ = fft(temp_surface2_filtered, sample_size);

figure("Name","Frequency Domain Representation of the Filtered Signals")
% Top: Surface 1 & Surface 2
subplot(2,2,1); stem(f, abs(fftshift(TEMP_SURFACE1_FILTERED_/sample_size))); title("Surface 1 Filtered Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);
subplot(2,2,2); stem(f, abs(fftshift(TEMP_SURFACE2_FILTERED_/sample_size))); title("Surface 2 Filtered Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);

% Bottom: Air 1 & Air 2
subplot(2,2,3); stem(f, abs(fftshift(TEMP_AIR1_FILTERED_/sample_size))); title("Air 1 Filtered Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);
subplot(2,2,4); stem(f, abs(fftshift(TEMP_AIR2_FILTERED_/sample_size))); title("Air 2 Filtered Spectrum"); xlabel("Frequency (Hz)"); ylabel("Magnitude"); ylim([0 2]); xlim([-0.05 0.05]);

%% Downsample and Scatter Plot

samples_per_min = round(fs * 60);
minute_indices = 1:samples_per_min:sample_size;
time_in_minutes = samples(minute_indices) / 60; 

% Extract data at indices
air1_min = temp_air1_filtered(minute_indices); % Air 1
air2_min = temp_air2_filtered(minute_indices); % Air 2
surface1_min = temp_surface1_filtered(minute_indices); % Surface 1
surface2_min = temp_surface2_filtered(minute_indices); % Surface 2

figure("Name","Temperature Data at Every Minute")
% Top: Surface 1 & Surface 2
subplot(2,2,1); 
scatter(time_in_minutes, surface1_min, 'filled'); 
title("Surface 1"); xlabel("Time (Minutes)"); ylabel("Temperature (C)"); grid on;

subplot(2,2,2); 
scatter(time_in_minutes, surface2_min, 'filled'); 
title("Surface 2"); xlabel("Time (Minutes)"); ylabel("Temperature (C)"); grid on;

% Bottom: Air 1 & Air 2
subplot(2,2,3); 
scatter(time_in_minutes, air1_min, 'filled'); 
title("Air 1"); xlabel("Time (Minutes)"); ylabel("Temperature (C)"); grid on;

subplot(2,2,4); 
scatter(time_in_minutes, air2_min, 'filled'); 
title("Air 2"); xlabel("Time (Minutes)"); ylabel("Temperature (C)"); grid on;

%% Export to Excel / CSV

% Create a table with named columns
T = table(time_in_minutes(:), ...
          air1_min(:), ... % Air 1
          air2_min(:), ... % Air 2
          surface1_min(:), ... % Surface 1
          surface2_min(:), ... % Surface 2
          'VariableNames', {'Time_Minutes', 'Air_1', 'Air_2', 'Surface_1', 'Surface_2'});

% Export
writetable(T, 'temperature_results.xlsx');
writetable(T, 'temperature_results.csv');

%% Error Calculation (Comparison with Chemical Engineering Data)

% Load reference data from chemical engineering students
ref_data = readtable('experiment_data.csv');

% Extract reference values (47 data points: minutes 0-46)
ref_time = ref_data.Time_min;
ref_air1 = ref_data.Air_1;
ref_air2 = ref_data.Air_2;
ref_surface1 = ref_data.Surface_1;
ref_surface2 = ref_data.Surface_2;

% Ensure we have matching data points (first 47 minutes: 0-46)
num_points = min(length(ref_time), length(time_in_minutes));

% Truncate our data to match reference data length (ensure column vectors)
our_air1 = air1_min(1:num_points); our_air1 = our_air1(:);
our_air2 = air2_min(1:num_points); our_air2 = our_air2(:);
our_surface1 = surface1_min(1:num_points); our_surface1 = our_surface1(:);
our_surface2 = surface2_min(1:num_points); our_surface2 = our_surface2(:);
our_time = time_in_minutes(1:num_points); our_time = our_time(:);

% Truncate reference data to match (ensure column vectors)
ref_air1 = ref_air1(1:num_points); ref_air1 = ref_air1(:);
ref_air2 = ref_air2(1:num_points); ref_air2 = ref_air2(:);
ref_surface1 = ref_surface1(1:num_points); ref_surface1 = ref_surface1(:);
ref_surface2 = ref_surface2(1:num_points); ref_surface2 = ref_surface2(:);

% Calculate absolute errors
error_air1 = abs(our_air1 - ref_air1);
error_air2 = abs(our_air2 - ref_air2);
error_surface1 = abs(our_surface1 - ref_surface1);
error_surface2 = abs(our_surface2 - ref_surface2);

% Calculate percentage errors (avoid division by zero)
pct_error_air1 = 100 * error_air1 ./ max(abs(ref_air1), 1e-6);
pct_error_air2 = 100 * error_air2 ./ max(abs(ref_air2), 1e-6);
pct_error_surface1 = 100 * error_surface1 ./ max(abs(ref_surface1), 1e-6);
pct_error_surface2 = 100 * error_surface2 ./ max(abs(ref_surface2), 1e-6);

% Ensure all error arrays are column vectors
error_air1 = error_air1(:);
error_air2 = error_air2(:);
error_surface1 = error_surface1(:);
error_surface2 = error_surface2(:);
pct_error_air1 = pct_error_air1(:);
pct_error_air2 = pct_error_air2(:);
pct_error_surface1 = pct_error_surface1(:);
pct_error_surface2 = pct_error_surface2(:);

% Calculate Mean Absolute Error (MAE) and Mean Percentage Error (MPE)
MAE_air1 = mean(error_air1);
MAE_air2 = mean(error_air2);
MAE_surface1 = mean(error_surface1);
MAE_surface2 = mean(error_surface2);

MPE_air1 = mean(pct_error_air1);
MPE_air2 = mean(pct_error_air2);
MPE_surface1 = mean(pct_error_surface1);
MPE_surface2 = mean(pct_error_surface2);

% Display error summary
fprintf('\n===== Error Analysis Summary =====\n');
fprintf('Sensor\t\t\tMAE (°C)\tMPE (%%)\n');
fprintf('Air 1\t\t\t%.2f\t\t%.2f\n', MAE_air1, MPE_air1);
fprintf('Air 2\t\t\t%.2f\t\t%.2f\n', MAE_air2, MPE_air2);
fprintf('Surface 1\t\t%.2f\t\t%.2f\n', MAE_surface1, MPE_surface1);
fprintf('Surface 2\t\t%.2f\t\t%.2f\n', MAE_surface2, MPE_surface2);
fprintf('==================================\n');

% Plot comparison: Our data vs Reference data
figure("Name","Comparison: Our Data vs Chemical Engineering Data")
% Top: Surface 1 & Surface 2
subplot(2,2,1);
plot(our_time, our_surface1, 'b-', 'LineWidth', 1.5); hold on;
plot(ref_time(1:num_points), ref_surface1, 'r--', 'LineWidth', 1.5);
title("Surface 1 Comparison"); xlabel("Time (Minutes)"); ylabel("Temperature (C)");
legend("Our Data", "Reference"); grid on;

subplot(2,2,2);
plot(our_time, our_surface2, 'b-', 'LineWidth', 1.5); hold on;
plot(ref_time(1:num_points), ref_surface2, 'r--', 'LineWidth', 1.5);
title("Surface 2 Comparison"); xlabel("Time (Minutes)"); ylabel("Temperature (C)");
legend("Our Data", "Reference"); grid on;

% Bottom: Air 1 & Air 2
subplot(2,2,3);
plot(our_time, our_air1, 'b-', 'LineWidth', 1.5); hold on;
plot(ref_time(1:num_points), ref_air1, 'r--', 'LineWidth', 1.5);
title("Air 1 Comparison"); xlabel("Time (Minutes)"); ylabel("Temperature (C)");
legend("Our Data", "Reference"); grid on;

subplot(2,2,4);
plot(our_time, our_air2, 'b-', 'LineWidth', 1.5); hold on;
plot(ref_time(1:num_points), ref_air2, 'r--', 'LineWidth', 1.5);
title("Air 2 Comparison"); xlabel("Time (Minutes)"); ylabel("Temperature (C)");
legend("Our Data", "Reference"); grid on;

% Plot percentage error over time
figure("Name","Percentage Error Over Time")
subplot(2,2,1);
plot(our_time, pct_error_surface1, 'LineWidth', 1.5);
title("Surface 1 Error"); xlabel("Time (Minutes)"); ylabel("Error (%)"); grid on;

subplot(2,2,2);
plot(our_time, pct_error_surface2, 'LineWidth', 1.5);
title("Surface 2 Error"); xlabel("Time (Minutes)"); ylabel("Error (%)"); grid on;

subplot(2,2,3);
plot(our_time, pct_error_air1, 'LineWidth', 1.5);
title("Air 1 Error"); xlabel("Time (Minutes)"); ylabel("Error (%)"); grid on;

subplot(2,2,4);
plot(our_time, pct_error_air2, 'LineWidth', 1.5);
title("Air 2 Error"); xlabel("Time (Minutes)"); ylabel("Error (%)"); grid on;

% Export error data to CSV
error_table = table(our_time, ...
    our_air1, ref_air1, error_air1, pct_error_air1, ...
    our_air2, ref_air2, error_air2, pct_error_air2, ...
    our_surface1, ref_surface1, error_surface1, pct_error_surface1, ...
    our_surface2, ref_surface2, error_surface2, pct_error_surface2, ...
    'VariableNames', {'Time_Min', ...
    'Our_Air1', 'Ref_Air1', 'AbsErr_Air1', 'PctErr_Air1', ...
    'Our_Air2', 'Ref_Air2', 'AbsErr_Air2', 'PctErr_Air2', ...
    'Our_Surface1', 'Ref_Surface1', 'AbsErr_Surface1', 'PctErr_Surface1', ...
    'Our_Surface2', 'Ref_Surface2', 'AbsErr_Surface2', 'PctErr_Surface2'});

writetable(error_table, 'error_analysis.csv');
writetable(error_table, 'error_analysis.xlsx');

%% Save All Figures as PNG

% Get all open figures
all_figs = findall(0, 'Type', 'figure');

% Get the directory where this script is located
script_dir = fileparts(mfilename('fullpath'));
figures_dir = fullfile(script_dir, 'figures');

% Create a folder for the figures if it doesn't exist
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% Save each figure
for i = 1:length(all_figs)
    fig = all_figs(i);
    fig_name = get(fig, 'Name');
    
    % Create a valid filename from the figure name
    if isempty(fig_name)
        filename = fullfile(figures_dir, sprintf('figure_%d.png', fig.Number));
    else
        % Replace spaces and special characters with underscores
        safe_name = regexprep(fig_name, '[^a-zA-Z0-9]', '_');
        filename = fullfile(figures_dir, sprintf('%s.png', safe_name));
    end
    
    % Save the figure
    saveas(fig, filename);
    fprintf('Saved: %s\n', filename);
end

fprintf('\nAll figures saved to the "figures" folder.\n');
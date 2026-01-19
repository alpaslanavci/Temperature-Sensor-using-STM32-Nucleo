clear all
close all
clc

% Time constants.
sample_size = 6001; 
fs = 2;
% t = (0:sample_size-1) / fs;

% Constants
V_SUPPLY = 3.3;       
R_FIXED = 20000.0;
T_KELVIN_0 = 273.15;

% Steinhart-Hart Coefficients  
A_COEFF = 0.4440809144e-3;
B_COEFF = 2.121308493e-4;
C_COEFF = 2.973400211e-7;

try
    % Create the serial port object
    mc = serialport("/dev/tty.usbmodem1303", 115200);
    mc.Timeout = 10;

    % Create empty matrices for each value of voltage, resistance and temperature.
    voltage_ch0 = zeros(sample_size, 1);
    voltage_ch1 = zeros(sample_size, 1);
    voltage_ch2 = zeros(sample_size, 1);
    voltage_ch3 = zeros(sample_size, 1);

    resistance_ch0 = zeros(sample_size, 1);
    resistance_ch1 = zeros(sample_size, 1);
    resistance_ch2 = zeros(sample_size, 1);
    resistance_ch3 = zeros(sample_size, 1);

    temp_ch0 = zeros(sample_size, 1);
    temp_ch1 = zeros(sample_size, 1);
    temp_ch2 = zeros(sample_size, 1);
    temp_ch3 = zeros(sample_size, 1);

    % Create empty plots first so we can update the plots on the go. For
    % convenience, I will only create plots for temperature.
    figure("Name","Temperature Plots");
    pos_temp_ch0 = subplot(4,1,1); plot_temp_ch0 = plot(pos_temp_ch0, NaN, NaN);
    pos_temp_ch1 = subplot(4,1,2); plot_temp_ch1 = plot(pos_temp_ch1, NaN, NaN);
    pos_temp_ch2 = subplot(4,1,3); plot_temp_ch2 = plot(pos_temp_ch2, NaN, NaN);
    pos_temp_ch3 = subplot(4,1,4); plot_temp_ch3 = plot(pos_temp_ch3, NaN, NaN);

    % Send 's' to start sampling.
    write(mc, "s", "uint8");
    disp("Started sampling")

    for i = 1:sample_size

        data = read(mc,4,"single"); % Read 4 8bit data. 1--> (Channel 0) 2--> (Channel 1) 3--> (Channel) 2 4--> (Channel 3)
        % Time vector for current iteration
        t = (0:i-1) / fs;

        % Update voltage values
        voltage_ch0(i) = data(1) * V_SUPPLY;
        voltage_ch1(i) = data(2) * V_SUPPLY;
        voltage_ch2(i) = data(3) * V_SUPPLY;
        voltage_ch3(i) = data(4) * V_SUPPLY;

        % Update resistance values
        resistance_ch0(i) = (voltage_ch0(i) * R_FIXED) / (V_SUPPLY - voltage_ch0(i));
        if resistance_ch0(i) <= 0
            resistance_ch0(i) = NaN;
        end
        resistance_ch1(i) = (voltage_ch1(i) * R_FIXED) / (V_SUPPLY - voltage_ch1(i));
        if resistance_ch1(i) <= 0
            resistance_ch1(i) = NaN;
        end
        resistance_ch2(i) = (voltage_ch2(i) * R_FIXED) / (V_SUPPLY - voltage_ch2(i));
        if resistance_ch2(i) <= 0
            resistance_ch2(i) = NaN;
        end
        resistance_ch3(i) = (voltage_ch3(i) * R_FIXED) / (V_SUPPLY - voltage_ch3(i));
        if resistance_ch3(i) <= 0
            resistance_ch3(i) = NaN;
        end

        % Update temperature values using Steinhart-Hart Equation.
        if ~isnan(resistance_ch0(i))
            ln_r0 = log(resistance_ch0(i));
            inv_T0 = A_COEFF + (B_COEFF * ln_r0) + (C_COEFF * ln_r0^3);
            temp_ch0(i) = (1.0 / inv_T0) - T_KELVIN_0;
        else
            temp_ch0(i) = NaN;
        end
        
        if ~isnan(resistance_ch1(i))
            ln_r1 = log(resistance_ch1(i));
            inv_T1 = A_COEFF + (B_COEFF * ln_r1) + (C_COEFF * ln_r1^3);
            temp_ch1(i) = (1.0 / inv_T1) - T_KELVIN_0;
        else
            temp_ch1(i) = NaN;
        end
        
        if ~isnan(resistance_ch2(i))
            ln_r2 = log(resistance_ch2(i));
            inv_T2 = A_COEFF + (B_COEFF * ln_r2) + (C_COEFF * ln_r2^3);
            temp_ch2(i) = (1.0 / inv_T2) - T_KELVIN_0;
        else
            temp_ch2(i) = NaN;
        end
        
        if ~isnan(resistance_ch3(i))
            ln_r3 = log(resistance_ch3(i));
            inv_T3 = A_COEFF + (B_COEFF * ln_r3) + (C_COEFF * ln_r3^3);
            temp_ch3(i) = (1.0 / inv_T3) - T_KELVIN_0;
        else
            temp_ch3(i) = NaN;
        end

        % Update temperature plots for this data
        set(plot_temp_ch0, "XData", t, "YData", temp_ch0(1:i));
        set(plot_temp_ch1, "XData", t, "YData", temp_ch1(1:i));
        set(plot_temp_ch2, "XData", t, "YData", temp_ch2(1:i));
        set(plot_temp_ch3, "XData", t, "YData", temp_ch3(1:i));

        fprintf("Sample %d, PA_0: %.1f -- PA_5: %.1f -- PA_3: %.1f -- PA_4:%.1f\n", i, temp_ch0(i), temp_ch1(i), temp_ch2(i), temp_ch3(i));
        
        drawnow
    end
catch ME
    disp(ME.message);
end
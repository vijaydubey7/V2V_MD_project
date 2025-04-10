clc; close all; restoredefaultpath; rehash toolboxcache;

%% Load XML Data
xmlFile = './data/test_highway_rate_0.5_trace.xml';
xmlData = xmlread(xmlFile);
timesteps = xmlData.getElementsByTagName('timestep');
vehicle_data = [];

%% Parameters
max_power = 1;                % Maximum transmission power
injectionRate = 0.5;          % Vehicle injection rate
num_iterations = 100;         % Number of iterations for optimization
f = 28e9;                     % Frequency in Hz (28 GHz mmWave)
threshold_distance = 1500;    % Max distance for beamforming (meters)
hysteresis_threshold = 5;     % Beam switching threshold

% Fixed Heights
RSU_height = 5;        % Height of RSU (meters)
vehicle_height = 1.5;  % Height of vehicles (meters)

%% Define Antenna Array
num_elements_az = 8;  % Number of elements in azimuth
num_elements_el = 8;  % Number of elements in elevation

antenna_array_az = phased.ULA('NumElements', num_elements_az, 'ElementSpacing', 0.5);
antenna_array_el = phased.ULA('NumElements', num_elements_el, 'ElementSpacing', 0.5);

% Extract number of elements
num_elements_az = antenna_array_az.NumElements;
num_elements_el = antenna_array_el.NumElements;

%% Parse Vehicle Data from XML
for i = 0:timesteps.getLength-1
    timestep = timesteps.item(i);
    vehicles = timestep.getElementsByTagName('vehicle');
    
    for j = 0:vehicles.getLength-1
        vehicle = vehicles.item(j);
        x = str2double(char(vehicle.getAttribute('x')));
        y = str2double(char(vehicle.getAttribute('y')));
        speed = str2double(char(vehicle.getAttribute('speed')));
        
        vehicle_data = [vehicle_data; x, y, vehicle_height, speed];
    end
end

%% Normalize Data
if ~isempty(vehicle_data)
    xPos = vehicle_data(:, 1);
    yPos = vehicle_data(:, 2);
    zPos = vehicle_data(:, 3); % Vehicle height (fixed at 1.5m)
    speeds = vehicle_data(:, 4);
    
    if max(speeds) > min(speeds) % Avoid division by zero
        speeds = (speeds - min(speeds)) / (max(speeds) - min(speeds));
    else
        speeds = zeros(size(speeds));
    end
else
    error('No vehicle data found in the XML file.');
end

%% Compute Network Parameters
numVehicles = size(vehicle_data, 1);
road_length = 1000; % Road length (meters)
road_width = 20;    % Road width (meters)
road_area = road_length * road_width;
vehicle_density = numVehicles / road_area;
speed_factor = mean(speeds) * 1.2; % Adjusted speed factor

fprintf('Number of vehicles for injection rate %.1f: %d\n', injectionRate, numVehicles);
% Initialize beam parameters
azimuth_angle = 0; elevation_angle = 0;
previous_angle = [azimuth_angle, elevation_angle];
all_SINR = []; all_throughput = []; all_outage = []; all_energy = [];
%% Simulation Loop
for iter = 1:min(num_iterations, numVehicles)
    % Random vehicle & target
    vehicle_idx = randi([1, numVehicles]);
    vehicle_pos = [xPos(vehicle_idx), yPos(vehicle_idx), zPos(vehicle_idx)];
    target_idx = randi([1, numVehicles]);
    while target_idx == vehicle_idx
        target_idx = randi([1, numVehicles]);
    end
    target_pos = [xPos(target_idx), yPos(target_idx), zPos(target_idx)];
    % Compute horizontal distance and elevation angle
    horizontal_distance = norm(vehicle_pos(1:2) - target_pos(1:2));
    height_difference = RSU_height - vehicle_pos(3);

    if horizontal_distance == 0
    elevation_angle = 0; % Prevent NaN issues
    else
        elevation_angle = atan2d(height_difference, horizontal_distance);
    end
    disp(['Iteration: ', num2str(iter)]); disp('vehicle_pos:'); disp(vehicle_pos);
    disp('target_pos:');  disp(target_pos);
    distance = norm(vehicle_pos - target_pos) + 1e-6;
    disp(['Distance: ', num2str(distance), ' m']); % Add this line
    % Calculate initial gain
    azimuth_angle = previous_angle(1);
    beam_gain = calculate_beam_gain(azimuth_angle, elevation_angle, num_elements_az, num_elements_el);
    environment = determine_environment(vehicle_density);
    
    % --- GWO-based Beamforming Optimization ---
    [optimized_beam_angles, power_allocation, beam_weights, SINR_values, throughput_values, outage_values, energy_values] = ...
        beamforming_with_GWO(previous_angle, true, vehicle_pos, target_pos, speed_factor, vehicle_density, num_elements_az, num_elements_el);
    
    % --- Hysteresis Beam Switching ---
    beam_switch = hysteresis_beam_switch(previous_angle, optimized_beam_angles, hysteresis_threshold);
    previous_angle = beam_switch;
    % --- Display Iteration Result ---
    fprintf('\nIteration %d:\n', iter);
    fprintf('Optimized Beam Angles: Azimuth = %.2f°, Elevation = %.2f°\n', optimized_beam_angles(1), optimized_beam_angles(2));
    fprintf('Power Allocation: %.2f W, Beam Switching Applied: [%.2f°, %.2f°]\n', power_allocation, beam_switch(1), beam_switch(2));
    fprintf('SINR: %.2f dB, Throughput: %.2f Mbps, Outage: %.2f%%, Energy: %.2f\n', SINR_values, throughput_values, outage_values*100, energy_values);
    % --- Mobility Update ---
    xPos(vehicle_idx) = vehicle_pos(1) + speeds(vehicle_idx) * 0.1;
    yPos(vehicle_idx) = vehicle_pos(2) + speeds(vehicle_idx) * 0.1;
end
%% Function Definitions
function SINR = calculate_SINR(beam_angles, power_allocation, vehicle_pos, target_pos, is_LoS, num_elements_az, num_elements_el)
    % Parameters
    noise_power = 5e-10;  % Noise Power in W
    distance = norm(vehicle_pos - target_pos);  % Distance between Tx and Rx
    is_LoS = distance < 1000;  % Line-of-Sight condition
    azimuth_angle = beam_angles(1);
    elevation_angle = beam_angles(2);
    beam_gain = calculate_beam_gain(azimuth_angle, elevation_angle, num_elements_az, num_elements_el);
    % Path Loss Calculation
    path_loss_dB = calculate_path_loss(distance, 28e9, is_LoS, beam_gain);
    path_loss_linear = 10^(path_loss_dB / 10);
     % Signal Power
    power_allocation = max(power_allocation, 0.05);  % Avoid zero
     signal_power = (power_allocation * beam_gain * 1.2) / path_loss_linear; % Increased gain factor
    % SINR Calculation with Minimum Threshold
    SINR_linear = signal_power / noise_power;
    SINR = 10 * log10(max(SINR_linear, 1e-5)); % Avoid near-zero values
    SINR = max(0, SINR); % Clamp minimum SINR to 0 dB (improves reliability)
    % Display
    fprintf('Distance: %.2f m, Path Loss: %.2f dB, Beam Gain: %.2f, Signal Power: %.6f W, SINR: %.2f dB\n', ...
        distance, path_loss_dB, beam_gain, signal_power, SINR);
end
function path_loss_dB = calculate_path_loss(distance, frequency, is_LoS, beam_gain)
    distance = max(distance, 1);  % Avoid zero
    f_GHz = frequency * 1e-9;  % GHz
    if nargin < 4
    beam_gain = 10; % Default beam gain in dB
    end
       % Convert beam gain from dB to linear scale
    gain_linear = 10^(beam_gain / 10);

    % Path Loss Model
    if is_LoS
        path_loss_dB = 32.0 + 18 * log10(f_GHz) + 18 * log10(distance);
    else
        path_loss_dB = 35.0 + 19 * log10(f_GHz) + 25 * log10(distance); % Adjusted NLoS formula
    end

    % Apply Beam Gain Correction Properly
    path_loss_dB = path_loss_dB - 10 * log10(gain_linear);

    % Ensure path loss is within a reasonable range (avoid extreme values)
    path_loss_dB = max(path_loss_dB, 30);
end
% Function for beam switching delay
function beam_switching_delay = calculate_beam_switching_delay(previous_angle, current_angle)
    threshold = 8; % Reduced threshold for smoother switching
    angle_change = norm(current_angle - previous_angle);
    delay_factor = 0.0015; % Reduced delay factor for faster switching

    if angle_change > threshold
        beam_switching_delay = delay_factor * angle_change;  
    else
        beam_switching_delay = (delay_factor / 2) * angle_change;  
    end

    beam_switching_delay = max(1e-6, min(beam_switching_delay, 1.0));
end
function throughput = calculate_throughput(SINR)
    bandwidth = 100e6; % 100 MHz
    SINR_linear = 10^(SINR / 10);
    throughput = bandwidth * log2(1 + SINR_linear) * 1e-6; % Convert to Mbps
    throughput = max(throughput, 0.1); % Ensure a minimum threshold
end

function outage_probability = calculate_outage_probability(SINR)
    outage_threshold = 10;
    slope = 0.5; % Adjusted for realistic transition
    outage_probability = 1 / (1 + exp(-slope * (SINR - outage_threshold)));
    outage_probability = max(0, min(outage_probability, 1)); % Constrain between 0 and 1
end

function energy_consumption = calculate_energy_consumption(power_allocation)
    T = 1; % Time duration
    efficiency = 0.5; % Improved efficiency factor
    energy_consumption = (1 / efficiency) * sum(power_allocation) * T;
    energy_consumption = max(0.05, energy_consumption); % Ensure nonzero energy
end

function packet_loss = calculate_packet_loss(SINR, throughput)
    if nargin < 2
        throughput = 10;
    end
    packet_loss = 0.08 * exp(-SINR / 6) * (1 - exp(-throughput / 25)); % Adjusted decay rates
    packet_loss = max(0, min(packet_loss, 0.1)); % Constrain packet loss to 10% max
end

function latency = calculate_latency(SINR, speed_factor)
    processing_delay = 0.03; % Reduced processing delay
    queuing_delay = 0.8 / (SINR + 3); % More gradual queuing response
    latency = processing_delay + queuing_delay / max(speed_factor, 0.5);
    latency = min(max(latency, 0.01), 0.2); % Bound latency between 10ms and 200ms
end

function jitter = calculate_jitter(latency)
    persistent prev_latency;
    if isempty(prev_latency)
        prev_latency = latency;
    end
    jitter = abs(latency - prev_latency);
    prev_latency = latency;
    jitter = min(max(jitter, 0.001), 0.05); % Bound jitter between 1ms and 50ms
end
%% Water-Filling Power Allocation Function
function power_allocation = water_filling(P, H)
     % P: Total power available
    % H: Channel gains 
    num_elements = length(H);
    power_allocation = zeros(num_elements, 1);
    [sorted_H, sorted_indices] = sort(H, 'descend');
    total_power = P;
    epsilon = 1e-6;  % Avoid divide-by-zero
    inverse_H = 1 ./ max(sorted_H, epsilon);
    
    % Water-Filling Iteration
    for k = 1:num_elements
        water_level = (total_power + sum(inverse_H(1:k))) / k;
        allocated_power = water_level - inverse_H(1:k);
        if all(allocated_power >= 0)
            break;
        end
    end
    
    power_allocation(sorted_indices(1:k)) = allocated_power;
    remaining_power = total_power - sum(allocated_power);
    
    % Distribute leftover power
    if remaining_power > 0
        power_allocation = power_allocation + (remaining_power / num_elements);
    end
    
    % Normalize if needed
    if sum(power_allocation) > total_power
        power_allocation = power_allocation * (total_power / sum(power_allocation));
    end
end

%% Create a Fuzzy logic
function selected_beam = fuzzy_beam_selection(SINR, path_loss, vehicle_speed, prev_beam, environment)
     % Bound input values within their respective ranges
    SINR = min(max(SINR, 0), 40);
    path_loss = min(max(path_loss, 0), 180);
    vehicle_speed = min(max(vehicle_speed, 0), 120);
    prev_beam = min(max(prev_beam, 0), 1);
    
    % Default environment
    if nargin < 5 || isempty(environment)
        environment = 'urban';
    end
    
    % Initialize FIS
    fis_name = strcat(environment, '_beam_selection');
    fis = mamfis('Name', fis_name);
    
    % Define Fuzzy Inputs
     fis = addInput(fis, [0 40], 'Name', 'SINR');
    fis = addMF(fis, 'SINR', 'trapmf', [0 0 8 15], 'Name', 'Low');
    fis = addMF(fis, 'SINR', 'trimf', [12 22 30], 'Name', 'Medium');
    fis = addMF(fis, 'SINR', 'trapmf', [25 30 40 40], 'Name', 'High');

    fis = addInput(fis, [0 180], 'Name', 'path_loss');
    if strcmp(environment, 'highway')
        fis = addMF(fis, 'path_loss', 'trapmf', [0 0 30 60], 'Name', 'Low');
        fis = addMF(fis, 'path_loss', 'trimf', [50 90 120], 'Name', 'Medium');
        fis = addMF(fis, 'path_loss', 'trapmf', [110 150 180 180], 'Name', 'High');
    else
        fis = addMF(fis, 'path_loss', 'trapmf', [0 0 40 80], 'Name', 'Low');
        fis = addMF(fis, 'path_loss', 'trimf', [70 100 140], 'Name', 'Medium');
        fis = addMF(fis, 'path_loss', 'trapmf', [130 170 180 180], 'Name', 'High');
    end

    fis = addInput(fis, [0 1], 'Name', 'prev_beam');
    fis = addMF(fis, 'prev_beam', 'trimf', [0 0.3 0.5], 'Name', 'Low Stability');
    fis = addMF(fis, 'prev_beam', 'trimf', [0.4 0.7 0.9], 'Name', 'Medium Stability');
    fis = addMF(fis, 'prev_beam', 'trimf', [0.8 1 1], 'Name', 'High Stability');

    % Output: Beam Selection
    fis = addOutput(fis, [0 1], 'Name', 'selected_beam');
    fis = addMF(fis, 'selected_beam', 'trimf', [0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'selected_beam', 'trimf', [0.4 0.7 0.9], 'Name', 'Medium');
    fis = addMF(fis, 'selected_beam', 'trimf', [0.8 1 1], 'Name', 'High');

    % Enhanced Rule Base
    ruleList = [
        1 3 1 1 0.95 1;
        1 2 2 1 0.75 2;
        2 1 3 2 0.85 1;
        3 1 3 3 1.0 1;
        2 2 2 2 0.9 2;
    ];
    fis = addRule(fis, ruleList);

    % Evaluate FIS
    selected_beam = evalfis(fis, [SINR, path_loss, prev_beam]);
    fprintf('Selected Beam Value: %.2f\n', selected_beam);
end

%% Running the Main Function
previous_angle = [0, 0]; optimize_beam_weights = true;
 [optimized_beam_angles, power_allocation, beam_weights, SINR_values, throughput_values, outage_values, energy_values] = ...
        beamforming_with_GWO(previous_angle, true, vehicle_pos, target_pos, speed_factor, vehicle_density, num_elements_az, num_elements_el);
fprintf('Optimized Beam Angle: %.2f degrees\n', optimized_beam_angles); fprintf('Power Allocation: %.2f watts\n', power_allocation); disp('beam_weights');disp(beam_weights);
current_angle = optimized_beam_angles;
hysteresis_threshold = 5;  % Angle threshold for switching
beam_switch = hysteresis_beam_switch(previous_angle, optimized_beam_angles, hysteresis_threshold);
if beam_switch
    disp('Switching beam...');
else
    disp('No beam switch required.');
end
%% Beamforming without GWO (Baseline)
function [beam_angles, power_allocation, SINR] = beamforming_no_GWO(vehicle_pos, target_pos, max_power, num_elements_az, num_elements_el)
    if nargin < 3, max_power = 1; end
    if nargin < 4, num_elements_az = 8; end
    if nargin < 5, num_elements_el = 8; end
    
   noise_power = 1e-12;  % Watts
    threshold_distance = 1000;  % LoS threshold in meters
    distance = norm(vehicle_pos - target_pos);
    is_LoS = distance < threshold_distance;
    % Calculate Beam Angle
    azimuth_angle = atan2d(target_pos(2) - vehicle_pos(2), target_pos(1) - vehicle_pos(1));
    azimuth_angle = max(-60, min(azimuth_angle, 60));  % Bound angle
    elevation_angle = atan2d(target_pos(3) - vehicle_pos(3), distance);
    elevation_angle = max(-30, min(elevation_angle, 30));
    beam_angles = [azimuth_angle, elevation_angle];
    % Beam Gain Calculation
    beam_gain = calculate_beam_gain(azimuth_angle, elevation_angle, num_elements_az, num_elements_el);
    % Path Loss Calculation
    path_loss = calculate_path_loss(distance, 28e9, is_LoS, beam_gain);
    % Signal Power Calculation
    signal_power = max(1e-9, (max_power * beam_gain) / (10^(path_loss / 10)));
    % SINR Calculation
    SINR = 10 * log10(signal_power / noise_power);
    % Dynamic Power Allocation
    target_SINR = 12; % Target SINR in dB
    gamma = 0.5; % Power adaptation factor
    power_allocation = max_power * (target_SINR / (SINR + 1e-6))^gamma;
    power_allocation = max(0.1, min(power_allocation, max_power));
end
%% 
function [optimized_beam_angles, power_allocation, beam_weights, SINR_values, throughput_values, outage_values, energy_values, packet_loss_values, latency_values, jitter_values] = ...
    beamforming_with_GWO(previous_angle, true, vehicle_pos, target_pos, speed_factor, vehicle_density, num_elements_az, num_elements_el)
    num_wolves = 10;  num_iterations = 100;
    optimize_beam_weights = true;
    theta_range = [-60, 60];  % Azimuth range
    phi_range = [-30, 30];    % Elevation range
    power_min = 0.1;
    power_max = min(1.0, dynamic_power_allocation(vehicle_density) + 0.5);
    weight_range = [0.1, 2.0];
    % Define antenna array parameters
    antenna_array_az = phased.ULA('NumElements', num_elements_az, 'ElementSpacing', 0.5);
    antenna_array_el = phased.ULA('NumElements', num_elements_el, 'ElementSpacing', 0.5);
    % Initialize Wolves
    wolves_positions = [ ...
        rand(num_wolves, 1) * diff(theta_range) + theta_range(1), ...
        rand(num_wolves, 1) * diff(phi_range) + phi_range(1), ...
        rand(num_wolves, 1) * (power_max - power_min) + power_min, ...
        rand(num_wolves, 1) * diff(weight_range) + weight_range(1) ...
    ];
    % Initialize Performance Metrics
    SINR_values = zeros(num_iterations, 1);
    throughput_values = zeros(num_iterations, 1);
    outage_values = zeros(num_iterations, 1);
    energy_values = zeros(num_iterations, 1);
    packet_loss_values = zeros(num_iterations, 1);
    latency_values = zeros(num_iterations, 1);
    jitter_values = zeros(num_iterations, 1);
    fitness = zeros(num_wolves, 1);
    % GWO Optimization Loop
    for iter = 1:num_iterations
        a = 2 * (1 - iter / num_iterations);
        % Evaluate fitness for all wolves
        for i = 1:num_wolves
            [fitness(i), SINR_values(iter), throughput_values(iter), outage_values(iter), energy_values(iter), ...
             packet_loss_values(iter), latency_values(iter), jitter_values(iter)] = ...
                fitness_function(wolves_positions(i, :), previous_angle, vehicle_pos, target_pos, speed_factor, vehicle_density, antenna_array_az, antenna_array_el);
        end
        % Sort Wolves by Fitness
        [~, idx] = sort(fitness, 'descend');
        alpha_wolf = wolves_positions(idx(1), :);
        beta_wolf = wolves_positions(idx(2), :);
        delta_wolf = wolves_positions(idx(3), :);
        % Update Wolves' Positions
        for i = 1:num_wolves
            A = 2 * a * rand(1, 4) - a;
            C = 2 * rand(1, 4);
            D_alpha = abs(C .* alpha_wolf - wolves_positions(i, :));
            D_beta  = abs(C .* beta_wolf  - wolves_positions(i, :));
            D_delta = abs(C .* delta_wolf - wolves_positions(i, :));
            X1 = alpha_wolf - A .* D_alpha + 0.02 * randn(1, 4);
            X2 = beta_wolf  - A .* D_beta  + 0.02 * randn(1, 4);
            X3 = delta_wolf - A .* D_delta + 0.02 * randn(1, 4);
            wolves_positions(i, :) = (X1 + X2 + X3) / 3;
            % Apply Constraints
            wolves_positions(i, 1) = max(theta_range(1), min(wolves_positions(i, 1), theta_range(2)));
            wolves_positions(i, 2) = max(phi_range(1), min(wolves_positions(i, 2), phi_range(2)));
            wolves_positions(i, 3) = max(power_min, min(wolves_positions(i, 3), power_max));
            if optimize_beam_weights
                wolves_positions(i, 4) = max(weight_range(1), min(wolves_positions(i, 4), weight_range(2)));
            else
                wolves_positions(i, 4) = 1.0;
            end
        end
        % Store Best Solution
        optimized_beam_angles = alpha_wolf(1:2);
        power_allocation = max(0.05, alpha_wolf(3));
        beam_weights = alpha_wolf(4);
    end
end

function [fitness, SINR, throughput, outage_prob, energy, packet_loss, latency, jitter] = ...
    fitness_function(params, previous_angle, vehicle_pos, target_pos, speed_factor, vehicle_density, antenna_array_az, antenna_array_el)

    azimuth = params(1);  elevation = params(2);  power = max(0.05, params(3));     
   num_elements_az = antenna_array_az.NumElements;
    num_elements_el = antenna_array_el.NumElements;

    distance = max(norm(vehicle_pos - target_pos), 1e-3);
    is_LoS = distance < 1000;

    % Beam Gain Calculation
    beam_gain = max(1.0, calculate_beam_gain(azimuth, elevation, num_elements_az, num_elements_el));

    % Path Loss Calculation
    path_loss = max(1e-3, calculate_path_loss(distance, 2.16e9, is_LoS, beam_gain));

    % SINR Calculation
    SINR = calculate_SINR([azimuth, elevation], power, vehicle_pos, target_pos, is_LoS, num_elements_az, num_elements_el);

    % Power Allocation
    power = max(0.1, water_filling(SINR, power));  

    % Adaptive Beam Weight Selection
    adaptive_weight = fuzzy_beam_selection(SINR, path_loss, speed_factor, 0.5, determine_environment(vehicle_density));
    beam_weight = 0.1 + adaptive_weight * (2.0 - 0.1);

    % Throughput, Outage, Energy, Packet Loss, Latency, and Jitter Calculation
    throughput = calculate_throughput(SINR);
    outage_prob = calculate_outage_probability(SINR);
    energy = calculate_energy_consumption(power);
    packet_loss = max(0, calculate_packet_loss(SINR, throughput));
    latency = calculate_latency(SINR, throughput);
    jitter = calculate_jitter(latency);

    % Normalization
    fitness = 2.0 * (SINR / (SINR + 1)) + 1.5 * (throughput / (throughput + 1)) - ...
              0.7 * (path_loss / (path_loss + 1)) - 1.5 * (energy / (energy + 1)) - ...
              0.2 * (norm(previous_angle - [azimuth, elevation])) - 0.05 * (azimuth^2 + elevation^2);
end

function environment = determine_environment(vehicle_density)
    if vehicle_density > 30  % Threshold can be tuned
        environment = 'urban';
    else
        environment = 'highway';
    end
end

function power_max = dynamic_power_allocation(vehicle_density)
    if vehicle_density < 0.8
        power_max = 1.5;
    elseif vehicle_density < 1.5
        power_max = 2.5;
    else
        power_max = 4.0;
    end
end

% ==== Hysteresis for Beam Switching ====
function beam_switch_angle = hysteresis_beam_switch(previous_angle, current_angle, threshold)
    if nargin < 3
        threshold = 3;
    end
    if abs(current_angle - previous_angle) > threshold
        beam_switch_angle = current_angle;
    else
        beam_switch_angle = previous_angle;
    end
end
function beam_gain = calculate_beam_gain(azimuth_angle, elevation_angle, num_elements_az, num_elements_el)
    % Analog Beamforming Calculation for Two ULAs
    max_gain_dBi = 15;  % Adjusted maximum gain in dBi
    HPBW_az = 12.8;  % Azimuth Half-Power Beamwidth
    HPBW_el = 12.8;  % Elevation Half-Power Beamwidth
    
    % Calculate attenuation based on angles
    attenuation_az = -12 * (azimuth_angle / HPBW_az)^2;
    attenuation_el = -12 * (elevation_angle / HPBW_el)^2;
    
    % Wavelength at 28 GHz (mmWave)
    wavelength = 3e8 / 28e9;
    d = wavelength / 2;  % Element spacing
    
    % Steering vector for Azimuth ULA
    steering_vector_az = exp(1i * 2 * pi * d * (0:num_elements_az-1) * sin(deg2rad(azimuth_angle)) / wavelength);
    
    % Steering vector for Elevation ULA
    steering_vector_el = exp(1i * 2 * pi * d * (0:num_elements_el-1) * sin(deg2rad(elevation_angle)) / wavelength);
    
    % Array Factor Calculation (Azimuth * Elevation)
    array_factor_az = abs(sum(steering_vector_az)) / num_elements_az; % Normalize azimuth
    array_factor_el = abs(sum(steering_vector_el)) / num_elements_el; % Normalize elevation
    array_factor = array_factor_az * array_factor_el;
    
    % Total Gain Calculation for Analog Beamforming
    total_gain_dBi = max_gain_dBi + attenuation_az + attenuation_el + 10 * log10(array_factor);
    total_gain_dBi = max(total_gain_dBi, 0);  % Avoid negative gain
    
    % Convert to linear scale
    beam_gain = 10^(total_gain_dBi / 10);
end
%% Simulation Loop
num_iterations = 100; previous_angle = [0, 0];
hysteresis_threshold = 5; optimize_beam_weights = true; default_gain = 0.1;
%% Simulation Parameters
theta = -60:1:60; % Azimuth angles
phi = -60:1:60;   % Elevation angles
[Theta, Phi] = meshgrid(theta, phi);

% Pre-allocation
azimuth_angles_no_GWO = zeros(1, num_iterations);
elevation_angles_no_GWO = zeros(1, num_iterations);
beam_gain_no_GWO = zeros(1, num_iterations);
SINR_no_GWO = zeros(1, num_iterations);
path_loss_no_GWO = zeros(1, num_iterations);
throughput_no_GWO = zeros(1, num_iterations);
outage_no_GWO = zeros(1, num_iterations);
energy_no_GWO = zeros(1, num_iterations);
beam_delay_no_GWO = zeros(1, num_iterations);
spectral_efficiency_no_GWO = zeros(1, num_iterations);
packet_loss_no_GWO = zeros(1, num_iterations);
latency_no_GWO = zeros(1, num_iterations);
jitter_no_GWO = zeros(1, num_iterations);

azimuth_angles_with_GWO = zeros(1, num_iterations);
elevation_angles_with_GWO = zeros(1, num_iterations);
beam_gain_with_GWO = zeros(1, num_iterations);
SINR_with_GWO = zeros(1, num_iterations);
path_loss_with_GWO = zeros(1, num_iterations);
throughput_with_GWO = zeros(1, num_iterations);
outage_with_GWO = zeros(1, num_iterations);
energy_with_GWO = zeros(1, num_iterations);
beam_delay_with_GWO = zeros(1, num_iterations);
spectral_efficiency_with_GWO = zeros(1, num_iterations);
packet_loss_with_GWO = zeros(1, num_iterations);
latency_with_GWO = zeros(1, num_iterations);
jitter_with_GWO = zeros(1, num_iterations);

%% Simulation Loop
previous_angle = [0, 0]; num_elements = 8;
target_pos_initial = [0, 0, 1.5]; % initial position
target_velocity = [1, 0.5, 0]; % velocity (m/s)
time_step = 0.1; % Time step (seconds)
%% Define Antenna Array
num_elements_az = 8;  % Number of elements in azimuth
num_elements_el = 8;  % Number of elements in elevation

antenna_array_az = phased.ULA('NumElements', num_elements_az, 'ElementSpacing', 0.5);
antenna_array_el = phased.ULA('NumElements', num_elements_el, 'ElementSpacing', 0.5);

% Extract number of elements
num_elements_az = antenna_array_az.NumElements;
num_elements_el = antenna_array_el.NumElements;

for i = 1:num_iterations
    target_idx = randi([1, numVehicles]);
    target_pos = [xPos(target_idx), yPos(target_idx), zPos(target_idx)];

    %current_time = i * time_step;
    %target_pos = target_pos_initial + target_velocity * current_time; % Predict target position
    vehicle_pos = [xPos(i), yPos(i), zPos(i)]; % Calculate vehicle position
    distance = norm(vehicle_pos - target_pos) + 1e-6; % Calculate distance
    is_LoS = distance < 1000;

    %% Without GWO
    [beam_angle_no_GWO, power_no_GWO] = beamforming_no_GWO(vehicle_pos, target_pos, 1);
    azimuth_angles_no_GWO(i) = beam_angle_no_GWO(1);
    elevation_angles_no_GWO(i) = beam_angle_no_GWO(2);
    beam_gain_no_GWO(i) = calculate_beam_gain(azimuth_angles_no_GWO(i), elevation_angles_no_GWO(i), num_elements_az, num_elements_el);
    SINR_no_GWO(i) = calculate_SINR([azimuth_angles_no_GWO(i), elevation_angles_no_GWO(i)], power_no_GWO, vehicle_pos, target_pos, is_LoS, num_elements_az, num_elements_el);
    throughput_no_GWO(i) = calculate_throughput(SINR_no_GWO(i));
    path_loss_no_GWO(i) = calculate_path_loss(distance, 28e9, is_LoS, beam_gain_no_GWO(i));
    outage_no_GWO(i) = calculate_outage_probability(SINR_no_GWO(i));
    energy_no_GWO(i) = calculate_energy_consumption(power_no_GWO);
    beam_delay_no_GWO(i) = calculate_beam_switching_delay(previous_angle, [azimuth_angles_no_GWO(i), elevation_angles_no_GWO(i)]);
    spectral_efficiency_no_GWO(i) = log2(1 + SINR_no_GWO(i));
    packet_loss_no_GWO(i) = calculate_packet_loss(SINR_no_GWO(i), throughput_no_GWO(i));
    latency_no_GWO(i) = calculate_latency(SINR_no_GWO(i), speed_factor);
    jitter_no_GWO(i) = calculate_jitter(latency_no_GWO(i));
    previous_angle = [azimuth_angles_no_GWO(i), elevation_angles_no_GWO(i)];

    % With GWO
    [optimized_beam_angles, power_allocation, ~, SINR_values, throughput_values, outage_values, energy_values] = ...
        beamforming_with_GWO(previous_angle, true, vehicle_pos, target_pos, speed_factor, vehicle_density, num_elements_az, num_elements_el);
    azimuth_angles_with_GWO(i) = optimized_beam_angles(1);
    elevation_angles_with_GWO(i) = optimized_beam_angles(2);
    beam_gain_with_GWO(i) = calculate_beam_gain(azimuth_angles_with_GWO(i), elevation_angles_with_GWO(i), num_elements_az, num_elements_el);
    SINR_with_GWO(i) = calculate_SINR(optimized_beam_angles, power_allocation, vehicle_pos, target_pos, is_LoS, num_elements_az, num_elements_el);
    throughput_with_GWO(i) = calculate_throughput(SINR_with_GWO(i));
    path_loss_with_GWO(i) = calculate_path_loss(distance, 28e9, is_LoS, beam_gain_with_GWO(i));
    outage_with_GWO(i) = calculate_outage_probability(SINR_with_GWO(i));
    energy_with_GWO(i) = calculate_energy_consumption(power_allocation);
    beam_delay_with_GWO(i) = calculate_beam_switching_delay(previous_angle, [azimuth_angles_with_GWO(i), elevation_angles_with_GWO(i)]);
    spectral_efficiency_with_GWO(i) = log2(1 + SINR_with_GWO(i));
    packet_loss_with_GWO(i) = calculate_packet_loss(SINR_with_GWO(i), throughput_with_GWO(i));
    latency_with_GWO(i) = calculate_latency(SINR_with_GWO(i), speed_factor);
    jitter_with_GWO(i) = calculate_jitter(latency_with_GWO(i));
    previous_angle = [azimuth_angles_with_GWO(i); elevation_angles_with_GWO(i)];

end

%% Beamforming Pattern Plotting
valid_no_GWO = azimuth_angles_no_GWO ~= 0;
valid_with_GWO = azimuth_angles_with_GWO ~= 0;

avg_azimuth_no_GWO = mean(azimuth_angles_no_GWO(valid_no_GWO));
avg_azimuth_with_GWO = mean(azimuth_angles_with_GWO(valid_with_GWO));

fprintf('\n--- Beamforming Summary ---\n');
fprintf('Average Azimuth (No GWO): %.2f°\n', avg_azimuth_no_GWO);
fprintf('Average Azimuth (With GWO): %.2f°\n', avg_azimuth_with_GWO);

% Index Mapping
theta_idx_no_GWO = round(interp1(theta, 1:length(theta), azimuth_angles_no_GWO(valid_no_GWO), 'nearest', 'extrap'));
phi_idx_no_GWO = round(interp1(phi, 1:length(phi), elevation_angles_no_GWO(valid_no_GWO), 'nearest', 'extrap'));
theta_idx_with_GWO = round(interp1(theta, 1:length(theta), azimuth_angles_with_GWO(valid_with_GWO), 'nearest', 'extrap'));
phi_idx_with_GWO = round(interp1(phi, 1:length(phi), elevation_angles_with_GWO(valid_with_GWO), 'nearest', 'extrap'));

% Clamping
theta_idx_no_GWO = max(1, min(length(theta), theta_idx_no_GWO));
phi_idx_no_GWO = max(1, min(length(phi), phi_idx_no_GWO));
theta_idx_with_GWO = max(1, min(length(theta), theta_idx_with_GWO));
phi_idx_with_GWO = max(1, min(length(phi), phi_idx_with_GWO));

% Gain Matrices
AF_without_GWO = zeros(size(Theta));
AF_with_GWO = zeros(size(Theta));
for i = 1:length(theta_idx_no_GWO)
    AF_without_GWO(phi_idx_no_GWO(i), theta_idx_no_GWO(i)) = 1.5;
end
for i = 1:length(theta_idx_with_GWO)
    AF_with_GWO(phi_idx_with_GWO(i), theta_idx_with_GWO(i)) = 1;
end

% Smooth and Normalize
AF_without_GWO = imgaussfilt(AF_without_GWO, 2.5);
AF_with_GWO = imgaussfilt(AF_with_GWO, 2.5);

% Equal normalization
max_no_GWO = max(AF_without_GWO(:));
max_with_GWO = max(AF_with_GWO(:));
AF_with_GWO = AF_with_GWO * (max_no_GWO / max_with_GWO);

% dB Conversion
gain_no_GWO_dB = 10 * log10(AF_without_GWO + 1e-6);
gain_with_GWO_dB = 10 * log10(AF_with_GWO + 1e-6);

%% Polar Pattern
beam_pattern_no_GWO = abs(sin(deg2rad(theta - avg_azimuth_no_GWO)));
beam_pattern_with_GWO = abs(sin(deg2rad(theta - avg_azimuth_with_GWO)));
beam_pattern_no_GWO = beam_pattern_no_GWO / max(beam_pattern_no_GWO);
beam_pattern_with_GWO = beam_pattern_with_GWO / max(beam_pattern_with_GWO);

figure;
polarplot(deg2rad(theta), beam_pattern_no_GWO, 'r', 'LineWidth', 1.5); hold on;
polarplot(deg2rad(theta), beam_pattern_with_GWO, 'b--', 'LineWidth', 1.5);
legend('Without GWO', 'With GWO');
title('Beam Radiation Pattern (Azimuth Cut)');
grid on;
%% Plot Beam Pattern Without GWO
figure;
surf(Theta, Phi, gain_no_GWO_dB, 'EdgeColor', 'none');
xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Gain (dB)');
title('Beam Gain Pattern WITHOUT GWO Optimization');
colorbar; axis tight; view(2);
clim([-20 0]); % Gain range for better visibility

%% Plot Beam Pattern With GWO
figure;
surf(Theta, Phi, gain_with_GWO_dB, 'EdgeColor', 'none');
xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Gain (dB)');
title('Beam Gain Pattern WITH GWO Optimization');
colorbar; axis tight; view(2);
clim([-20 0]);  % Same range for comparison

%% Check Results
fprintf('Average Azimuth without GWO: %.2f°\n', avg_azimuth_no_GWO);
fprintf('Average Azimuth with GWO: %.2f°\n', avg_azimuth_with_GWO);
fprintf('Mean Path Loss (no GWO): %.2f dB\n', mean(path_loss_no_GWO));
fprintf('Mean Path Loss (with GWO): %.2f dB\n', mean(path_loss_with_GWO));
fprintf('Energy Consumption (no GWO): %.4f\n', mean(energy_no_GWO));
fprintf('Energy Consumption (with GWO): %.4f\n', mean(energy_with_GWO));
fprintf('Outage Probability (no GWO): %.4f\n', mean(outage_no_GWO));
fprintf('Outage Probability (with GWO): %.4f\n', mean(outage_with_GWO));
%% Debug Angle Visualization
figure; 
plot(azimuth_angles_no_GWO, 'r--'); hold on;
plot(azimuth_angles_with_GWO, 'b-');

% Manually set the legend entries
legend({'Azimuth No GWO', 'Azimuth With GWO'}, 'Location', 'best');
title('Azimuth Angle Comparison');

figure;
plot(elevation_angles_no_GWO, 'g--'); hold on;
plot(elevation_angles_with_GWO, 'm-');
legend({'Elevation No GWO', 'Elevation With GWO'}, 'Location', 'best');
title('Elevation Angle Comparison');
% 3D Beamforming Plots
figure;
surf(Theta, Phi, AF_without_GWO);
xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Beam Gain');
title('3D Beamforming Pattern - Without GWO');
shading interp; colormap jet; colorbar; view(3);

figure;
surf(Theta, Phi, AF_with_GWO);
xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Beam Gain');
title('3D Beamforming Pattern - With GWO');
shading interp; colormap jet; colorbar; view(3);

% Beam Gain Comparison Plot
figure;
plot(theta, gain_no_GWO_dB, 'b--', 'LineWidth', 2); hold on;
plot(theta, gain_with_GWO_dB, 'r', 'LineWidth', 2);
xlabel('Azimuth Angle (°)');
ylabel('Beam Gain (dB)');
legend('Without GWO', 'With GWO');
title('Beam Gain Comparison (Azimuth Cut)');
grid on;

% Polar Plot for Beam Patterns
figure;
polarplot(deg2rad(theta), beam_pattern_no_GWO, 'b--', 'LineWidth', 2); hold on;
polarplot(deg2rad(theta), beam_pattern_with_GWO, 'r', 'LineWidth', 2);
legend('Without GWO', 'With GWO');
title('Beam Radiation Pattern (Azimuth Cut)');

% Plot without GWO
subplot(1, 2, 1);
surf(Theta, Phi, AF_without_GWO); % Plot Amplitude without GWO
title('Beamforming Without GWO'); xlabel('Azimuth Angle (°)'); ylabel('Elevation Angle (°)'); zlabel('Amplitude');
shading interp; colorbar;

% Plot with GWO
subplot(1, 2, 2);
surf(Theta, Phi, AF_with_GWO); % Plot Amplitude with GWO
title('Beamforming With GWO'); xlabel('Azimuth Angle (°)'); ylabel('Elevation Angle (°)'); zlabel('Amplitude');
shading interp; colorbar;

figure;
scatter3(SINR_no_GWO, path_loss_no_GWO, throughput_no_GWO, 'ro', 'filled'); % Red for Without GWO
hold on;
scatter3(SINR_with_GWO, path_loss_with_GWO, throughput_with_GWO, 'bo', 'filled'); % Blue for With GWO
title('3D Beamforming: With and Without GWO');  
xlabel('SINR (dB)'); ylabel('Path Loss (dB)'); zlabel('Throughput (bps)');
legend('Without GWO', 'With GWO');
grid on;
view(3);

%% Plotting Results
figure;

% SINR Comparison
subplot(3, 3, 1);
plot(1:num_iterations, SINR_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, SINR_with_GWO, 'b-', 'LineWidth', 2);
title('SINR Comparison'); xlabel('Simulation'); ylabel('SINR (dB)');
legend('Without GWO', 'With GWO');
grid on;

% Path Loss Comparison
subplot(3, 3, 2);
plot(1:num_iterations, path_loss_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, path_loss_with_GWO, 'b-', 'LineWidth', 2);
title('Path Loss Comparison'); xlabel('Simulation'); ylabel('Path Loss (dB)');
legend('Without GWO', 'With GWO');
grid on;

% Outage Probability
subplot(3, 3, 3);
plot(1:num_iterations, outage_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, outage_with_GWO, 'b-', 'LineWidth', 2);
title('Outage Probability'); xlabel('Simulation'); ylabel('Outage Probability');
legend('Without GWO', 'With GWO');
grid on;

% Throughput Comparison
subplot(3, 3, 4);
plot(1:num_iterations, throughput_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, throughput_with_GWO, 'b-', 'LineWidth', 2);
title('Throughput Comparison'); xlabel('Simulation'); ylabel('Throughput (bps)');
legend('Without GWO', 'With GWO');
grid on;

% Energy Consumption
subplot(3, 3, 5);
plot(1:num_iterations, energy_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, energy_with_GWO, 'b-', 'LineWidth', 2);
title('Energy Consumption'); xlabel('Simulation'); ylabel('Energy (J)');
legend('Without GWO', 'With GWO');
grid on;

% Spectral Efficiency
subplot(3, 3, 6);
plot(1:num_iterations, spectral_efficiency_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, spectral_efficiency_with_GWO, 'b-', 'LineWidth', 2);
title('Spectral Efficiency'); xlabel('Simulation'); ylabel('Spectral Efficiency (bps/Hz)');
legend('Without GWO', 'With GWO');
grid on;

% Beam Switching Delay
subplot(3, 3, 7);
plot(1:num_iterations, beam_delay_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, beam_delay_with_GWO, 'b-', 'LineWidth', 2);
title('Beam Switching Delay'); xlabel('Simulation'); ylabel('Delay (ms)');
legend('Without GWO', 'With GWO');
grid on;

% Packet Loss
subplot(3, 3, 8);
plot(1:num_iterations, packet_loss_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, packet_loss_with_GWO, 'b-', 'LineWidth', 2);
title('Packet Loss'); xlabel('Simulation'); ylabel('Packet Loss (%)');
legend('Without GWO', 'With GWO');
grid on;

% Latency & Jitter Comparison
subplot(3, 3, 9);
plot(1:num_iterations, latency_no_GWO, 'r-', 'LineWidth', 2);
hold on;
plot(1:num_iterations, latency_with_GWO, 'b-', 'LineWidth', 2);
plot(1:num_iterations, jitter_no_GWO, 'g--', 'LineWidth', 2);
plot(1:num_iterations, jitter_with_GWO, 'k--', 'LineWidth', 2);
title('Latency & Jitter Comparison'); xlabel('Simulation'); ylabel('Time (ms)');
legend('Latency (No GWO)', 'Latency (GWO)', 'Jitter (No GWO)', 'Jitter (GWO)');
grid on;


% *3. SINR Distribution*
figure;
histogram(SINR_no_GWO, 'FaceColor', 'r', 'FaceAlpha', 0.5);
hold on;
histogram(SINR_with_GWO, 'FaceColor', 'b', 'FaceAlpha', 0.5);
xlabel('SINR (dB)');
ylabel('Frequency');
title('SINR Distribution');
legend('Without GWO', 'With GWO');
grid on;

% *4. Cumulative Distribution Function (CDF) of SINR*
figure;
cdfplot(SINR_no_GWO);
hold on;
cdfplot(SINR_with_GWO);
xlabel('SINR (dB)');
ylabel('CDF');
title('CDF of SINR');
legend('Without GWO', 'With GWO');
grid on;

figure;
plot(path_loss_no_GWO, SINR_no_GWO, 'ro-', 'LineWidth', 2);
hold on;
plot(path_loss_with_GWO, SINR_with_GWO, 'bo-', 'LineWidth', 2);
xlabel('Path Loss (dB)');
ylabel('SINR (dB)');
title('Path Loss vs SINR');
legend('Without GWO', 'With GWO');
grid on;

%Beam Radiation Patterns
figure;
subplot(1, 2, 1);
polarplot(deg2rad(theta), beam_pattern_no_GWO, 'b--', 'LineWidth', 2);
title('Without GWO');
legend('Without GWO');

subplot(1, 2, 2);
polarplot(deg2rad(theta), beam_pattern_with_GWO, 'r--', 'LineWidth', 2);
title('With GWO');
legend('With GWO');

sgtitle('Comparison of Beam Radiation Patterns');

figure;
hold on;
for i = 1:length(azimuth_angles_with_GWO)
    plot3(azimuth_angles_with_GWO(i), elevation_angles_with_GWO(i), beam_gain_with_GWO(i), 'ro', 'MarkerSize', 8, 'DisplayName', 'Optimized Beam');
end
grid on;
xlabel('Azimuth Angle (°)');
ylabel('Elevation Angle (°)');
zlabel('Beam Gain');
title('Dynamic Beam Steering Optimization');
legend;
hold off;

% *8. Average SINR per Beam Angle*
figure;
plot(azimuth_angles_no_GWO, SINR_no_GWO, 'o-', 'LineWidth', 2);
hold on;
plot(azimuth_angles_with_GWO, SINR_with_GWO, 'bo-', 'LineWidth', 2);
xlabel('Azimuth Angle (°)');
ylabel('Average SINR (dB)');
title('Average SINR per Beam Angle');
legend('Without GWO', 'With GWO');
grid on;

% *10. SINR vs Throughput*
if ~any(isnan(SINR_no_GWO)) && ~any(isnan(SINR_with_GWO))
    figure;
    plot(1:num_iterations, SINR_no_GWO, 'r', 'LineWidth', 2); 
    hold on;
    plot(1:num_iterations, SINR_with_GWO, 'b', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('SINR (dB)');
    title('SINR vs Iteration (Without and With GWO)');
    legend('Without GWO', 'With GWO');
end
disp('Simulation completed successfully!');
% Plot the azimuth and elevation angles over time
figure;
subplot(2,1,1);
plot(1:num_iterations, azimuth_angles_with_GWO, '-b', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Azimuth Angle (radians)');
title('Real-Time Azimuth Angle Feedback');

subplot(2,1,2);
plot(1:num_iterations, elevation_angles_with_GWO, '-r', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Elevation Angle (radians)');
title('Real-Time Elevation Angle Feedback');

% Compute averages for each metric
avg_SINR_no_GWO = mean(SINR_no_GWO); avg_SINR_with_GWO = mean(SINR_with_GWO);
avg_path_loss_no_GWO = mean(path_loss_no_GWO); avg_path_loss_with_GWO = mean(path_loss_with_GWO);
avg_throughput_no_GWO = mean(throughput_no_GWO); avg_throughput_with_GWO = mean(throughput_with_GWO);
avg_outage_no_GWO = mean(outage_no_GWO); avg_outage_with_GWO = mean(outage_with_GWO);
avg_energy_no_GWO = mean(energy_no_GWO); avg_energy_with_GWO = mean(energy_with_GWO);
avg_beam_delay_no_GWO = mean(beam_delay_no_GWO); avg_beam_delay_with_GWO = mean(beam_delay_with_GWO);
avg_packet_loss_no_GWO = mean(packet_loss_no_GWO); avg_packet_loss_with_GWO = mean(packet_loss_with_GWO);
avg_latency_no_GWO = mean(latency_no_GWO); avg_latency_with_GWO = mean(latency_with_GWO);
avg_jitter_no_GWO = mean(jitter_no_GWO); avg_jitter_with_GWO = mean(jitter_with_GWO);

% Create the table
summary_table = table(["Without GWO"; "With GWO"], ...
                      [avg_SINR_no_GWO; avg_SINR_with_GWO], ...
                      [avg_path_loss_no_GWO; avg_path_loss_with_GWO], ...
                      [avg_throughput_no_GWO; avg_throughput_with_GWO], ...
                      [avg_outage_no_GWO; avg_outage_with_GWO], ...
                      [avg_energy_no_GWO; avg_energy_with_GWO], ...
                      [avg_beam_delay_no_GWO; avg_beam_delay_with_GWO], ...
                      [avg_packet_loss_no_GWO; avg_packet_loss_with_GWO], ...
                      [avg_latency_no_GWO; avg_latency_with_GWO], ...
                      [avg_jitter_no_GWO; avg_jitter_with_GWO], ...
                      'VariableNames', {'Scenario', 'Avg_SINR', 'Avg_Path_Loss', 'Avg_Throughput', 'Avg_Outage_Probability', 'Avg_Energy_Consumption', 'Avg_Beam_Switching_Delay', 'Avg_Packet_Loss', 'Avg_Latency', 'Avg_Jitter'});

% Display the table
disp(summary_table);
fprintf('Iteration %d: Selected Beam = %.2f°, Power = %.2fW\n', num_iterations, optimized_beam_angles, power_allocation);
% Ensure optimized_beam_angles has the same number of rows as vehicle_pos
num_vehicles = size(vehicle_pos, 1);

if size(optimized_beam_angles, 1) ~= num_vehicles
    optimized_beam_angles = repmat(optimized_beam_angles(1), num_vehicles, 1);  % Repeat for all vehicles
end

% Quiver plot for beam directions
quiver(vehicle_pos(:,1), vehicle_pos(:,2), cosd(optimized_beam_angles(:,1)), sind(optimized_beam_angles(:,1)), 0.5, 'LineWidth', 2, 'Color', 'b');

% Hold plot for additional elements
hold on;

% Plot the target position
plot(target_pos(1), target_pos(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

% Add title and legend
title('Optimized Beam Direction');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend({'Beam Directions', 'Target'}, 'Location', 'best');

% Enhance visualization
axis equal;
grid on;
hold off;
clear; clc; close all;

% Injection Rates and Corresponding XML Files
injection_rates = [0.5, 1.0, 1.5, 2.0];
xml_files = {'./test_urban_rate_0.5_trace.xml', ...
             './test_urban_rate_1_trace.xml', ...
             './test_urban_rate_1.5_trace.xml', ...
             './test_urban_rate_2_trace.xml'};
% Arrays to store results
silhouette_scores_kmeans = zeros(1, length(injection_rates));
silhouette_scores_fcm = zeros(1, length(injection_rates));
ch_lifetime_kmeans = cell(1, length(injection_rates));
ch_lifetime_fcm = cell(1, length(injection_rates));
ch_stability_kmeans = cell(1, length(injection_rates));
ch_stability_fcm = cell(1, length(injection_rates));
optimal_clusters = zeros(1, length(injection_rates));
cluster_member_dist_kmeans = cell(1, length(injection_rates));
cluster_member_dist_fcm = cell(1, length(injection_rates));

% Loop through each injection rate and XML file
for i = 1:length(injection_rates)
    filename = xml_files{i};
    data = parse_fcd_xml_to_df(filename);  % Ensure this function is defined or available
    x = data(:, 1);
    y = data(:, 2);
    speeds = data(:, 3);  % Assuming speed is the third column in the XML data
    data_points = [x, y];  % Assuming data_points is a combination of x and y

    % Dynamic Cluster Count Calculation
    vehicle_density = length(x) / (max(x) * max(y));
    num_clusters = max(2, round(vehicle_density * 100));
    
    % --- K-Means Clustering ---
    [idx_kmeans, centroids_kmeans] = kmeans(data_points, num_clusters);
    cluster_member_dist_kmeans{i} = histcounts(idx_kmeans, num_clusters);
    [optimized_centroids, optimized_clusters] = gwo_centroid_optimization(data_points, num_clusters);
    % --- FCM Clustering ---
    [centers_fcm, U, ~] = fcm([x, y], num_clusters);
    [~, idx_fcm] = max(U, [], 1);
    cluster_member_dist_fcm{i} = histcounts(idx_fcm, num_clusters);
    centroids_fcm = zeros(num_clusters, 2); % Initialize centroids array
    for k = 1:num_clusters
        % Calculate weighted sum for each centroid (using U^2 as the weight)
        weighted_points = bsxfun(@times, [x, y], U(k, :)'.^2);  % Element-wise multiplication
        numerator = sum(weighted_points, 1);  % Sum along the rows
        denominator = sum(U(k, :).^2);  % Sum of membership weights (U^2)
        centroids_fcm(k, :) = numerator / denominator;  % Compute centroid
    end
    
    % Silhouette Scores (with error checking)
    if num_clusters > 1
        silhouette_scores_kmeans(i) = mean(silhouette([x, y], idx_kmeans));
        silhouette_scores_fcm(i) = mean(silhouette([x, y], idx_fcm'));
    else
        warning('Insufficient clusters for silhouette score calculation at rate %.1f', injection_rates(i));
    end
    
    % --- Cluster Head Selection ---
    candidate_CHs_kmeans = find_cluster_heads(data_points, idx_kmeans, centroids_kmeans);
    disp('K-means candidate CHs:');
    disp(candidate_CHs_kmeans);
    candidate_CHs_fcm = find_cluster_heads(data_points, idx_fcm, centroids_fcm);

    % --- Calculate CH Lifetime and Stability for K-Means ---
    lifetime_kmeans = calculate_ch_lifetime(candidate_CHs_kmeans, data_points, idx_kmeans, centroids_kmeans, speeds);
    stability_kmeans = calculate_ch_stability(candidate_CHs_kmeans, data_points, idx_kmeans, speeds);

    
    % Store Results for K-Means
    ch_lifetime_kmeans{i} = lifetime_kmeans;
    ch_stability_kmeans{i} = stability_kmeans;
    
    % --- Calculate CH Lifetime and Stability for FCM ---
    lifetime_fcm = calculate_ch_lifetime(candidate_CHs_fcm, data_points, idx_fcm, centroids_fcm, speeds);
    stability_fcm = calculate_ch_stability(candidate_CHs_fcm, data_points, idx_fcm, speeds);

    
    % Store Results for FCM
    ch_lifetime_fcm{i} = lifetime_fcm;
    ch_stability_fcm{i} = stability_fcm;
    
    % Store Optimal Cluster Count
    optimal_clusters(i) = num_clusters;
    
   % Visualization for K-Means
    figure;
    gscatter(x, y, idx_kmeans, 'rgbcmk', 'o', 8);  % Scatter plot of data points colored by cluster
    hold on;
    
    % Plot centroids of K-Means
    plot(centroids_kmeans(:,1), centroids_kmeans(:,2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
    
    % Ensure candidate_CHs is defined before plotting
    if exist('candidate_CHs', 'var') && ~isempty(candidate_CHs)
        plot(x(candidate_CHs), y(candidate_CHs), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    else
        disp('Warning: candidate_CHs is empty or undefined for K-Means.');
    end
    
    title(sprintf('K-Means Clustering (Rate: %.1f)', injection_rates(i)));
    legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids', 'Cluster Heads');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    grid on;
    hold off;
    
    % Visualization for FCM
    figure;
    gscatter(x, y, idx_fcm, 'rgbcmk', 'o', 8);  % Scatter plot of data points colored by cluster
    hold on;
    
    % Plot centroids of FCM
    plot(centers_fcm(:,1), centers_fcm(:,2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
    
    title(sprintf('FCM Clustering (Rate: %.1f)', injection_rates(i)));
    legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    grid on;
    hold off;

end
% Plot cluster heads for K-means
figure;
hold on;
plot(x, y, 'bo');  % Plot all data points
plot(x(candidate_CHs_kmeans), y(candidate_CHs_kmeans), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); % Plot K-means CHs
title('K-means Cluster Heads');
xlabel('X-coordinate');
ylabel('Y-coordinate');
hold off;

% Plot cluster heads for FCM
figure;
hold on;
plot(x, y, 'bo');  % Plot all data points
plot(x(candidate_CHs_fcm), y(candidate_CHs_fcm), 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g'); % Plot FCM CHs
title('FCM Cluster Heads');
xlabel('X-coordinate');
ylabel('Y-coordinate');
hold off;
% --- Plot Silhouette Scores ---
figure;
plot(injection_rates, silhouette_scores_kmeans, '-o', 'LineWidth', 2);
hold on;
plot(injection_rates, silhouette_scores_fcm, '-s', 'LineWidth', 2);
title('Silhouette Scores for K-Means and FCM');
xlabel('Injection Rates');
ylabel('Silhouette Score');
legend('K-Means', 'FCM');
grid on;
hold off;

% --- Plot Optimal Clusters ---
figure;
plot(injection_rates, optimal_clusters, '-^', 'LineWidth', 2);
title('Optimal Cluster Count vs Injection Rates');
xlabel('Injection Rates');
ylabel('Optimal Cluster Count');
grid on;

% --- Cluster Member Distribution ---
figure;
for i = 1:length(injection_rates)
    subplot(2, 2, i);
    bar([cluster_member_dist_kmeans{i}; cluster_member_dist_fcm{i}]');
    title(sprintf('Cluster Distribution (Rate: %.1f)', injection_rates(i)));
    xlabel('Clusters');
    ylabel('Number of Members');
    legend('K-Means', 'FCM');
end

% --- Plot CH Lifetime and Stability ---
figure;
for i = 1:length(injection_rates)
    subplot(2, 2, i);
    bar([mean(ch_lifetime_kmeans{i}), mean(ch_lifetime_fcm{i})]);
    title(sprintf('Avg CH Lifetime (Rate: %.1f)', injection_rates(i)));
    xlabel('Clustering Method');
    ylabel('Avg CH Lifetime');
    set(gca, 'XTickLabel', {'K-Means', 'FCM'});
    grid on;
end
sgtitle('Average CH Lifetime across Injection Rates');
figure;
for i = 1:length(injection_rates)
    subplot(2, 2, i);
    bar([mean(ch_stability_kmeans{i}), mean(ch_stability_fcm{i})]);
    title(sprintf('Avg CH Stability (Rate: %.1f)', injection_rates(i)));
    xlabel('Clustering Method');
    ylabel('Avg CH Stability');
    set(gca, 'XTickLabel', {'K-Means', 'FCM'});
    grid on;
end
sgtitle('Average CH Stability across Injection Rates');

%% Function Definitions
function candidate_CHs = find_cluster_heads(data_points, cluster_indices, centroids)
    % Initialize an empty column vector to store candidate cluster heads
    candidate_CHs = [];
    
    % Number of clusters from the centroids
    num_clusters = size(centroids, 1);
    
    for k = 1:num_clusters
        % Find the points in the current cluster
        cluster_members = find(cluster_indices == k);
        
        % Debug: Check if this cluster has members
        if isempty(cluster_members)
            disp(['Cluster ', num2str(k), ' has no members. Skipping...']);
            continue; % Skip if no points are assigned to this cluster
        end
        
        cluster_points = data_points(cluster_members, :);
        
        % Debug: Output the number of points in the current cluster
        disp(['Cluster ', num2str(k), ' has ', num2str(length(cluster_members)), ' members.']);
        
        % Calculate the distances from each point to the centroid of the cluster
        distances = sqrt(sum((cluster_points - centroids(k, :)).^2, 2));
        
        % Debug: Output the distances
        disp(['Distances from cluster points to centroid:']);
        disp(distances);
        
        % Compute the LOS scores (e.g., inverse distance) and connectivity (dummy example)
        los_scores = 1 ./ (1 + distances);
        
        % Debug: Output the LOS scores
        disp(['LOS scores for cluster ', num2str(k), ':']);
        disp(los_scores);
        
        % Assuming connectivity is simply the density of nearby points in the cluster (dummy calculation)
        connectivity_scores = arrayfun(@(i) sum(sqrt(sum((cluster_points - cluster_points(i, :)).^2, 2)) < 100), 1:length(cluster_members));
        
        % Debug: Output the connectivity scores
        disp(['Connectivity scores for cluster ', num2str(k), ':']);
        disp(connectivity_scores);
        
        % Combine the LOS and connectivity scores to determine candidate CH
        combined_scores = 0.6 * los_scores + 0.4 * connectivity_scores;
        
        % Debug: Output the combined scores
        disp(['Combined scores for cluster ', num2str(k), ':']);
        disp(combined_scores);
        
        % Select the candidate CH with the highest combined score
        [~, best_idx] = max(combined_scores);
        
        % Get the candidate CH from the cluster_members list
        final_CH = cluster_members(best_idx);
        
        % Ensure final_CH is a scalar index (if it's a vector, select the first element)
        if numel(final_CH) > 1
            disp(['Selecting the first index from final_CH: ', num2str(final_CH(1))]);
            final_CH = final_CH(1);  % Select the first index if multiple values
        end
        
        % Debugging step: Check the type and size of final_CH before concatenation
        disp(['Selected candidate CH for cluster ', num2str(k), ': ', num2str(final_CH)]);
        
        % Concatenate the candidate CH (final_CH) as a scalar index to candidate_CHs
        % Check if candidate_CHs is empty and initialize it accordingly
        if isempty(candidate_CHs)
            candidate_CHs = final_CH;  % First element
        else
            candidate_CHs = [candidate_CHs; final_CH];  % Add to list
        end
        
        % Debug: Confirm candidate CHs are being added
        disp(['Candidate CHs after cluster ', num2str(k), ':']);
        disp(candidate_CHs);
    end
    
    % Final debug: Confirm the list of candidate CHs at the end
    disp('Final candidate CHs:');
    disp(candidate_CHs);
end

% CH Lifetime Calculation
function lifetime = calculate_ch_lifetime(candidate_CHs, data_points, cluster_indices, centroids, speeds)
    lifetime = zeros(length(candidate_CHs), 1);
    for i = 1:length(candidate_CHs)
        ch_index = candidate_CHs(i);
        
        % Validate index
        if ch_index <= 0 || ch_index > size(data_points, 1)
            warning('Invalid index %d. Skipping this candidate CH.', ch_index);
            lifetime(i) = NaN;  % Assign NaN for invalid CHs
            continue;
        end
        
        ch_coordinates = data_points(ch_index, :);
        cluster_id = cluster_indices(ch_index);
        centroid_coordinates = centroids(cluster_id, :);
        distance_to_centroid = sqrt(sum((ch_coordinates - centroid_coordinates).^2));
        
        % Influence of speed on CH lifetime (optional)
        speed_influence = max(0, 1 - 0.01 * speeds(ch_index));  % Modify the influence as needed
        connectivity_score = sum(sqrt(sum((data_points - ch_coordinates).^2, 2)) < 100);
        
        lifetime(i) = max(0, 100 - distance_to_centroid + 0.5 * connectivity_score + speed_influence);
    end
end


% CH Stability Calculation
function stability = calculate_ch_stability(candidate_CHs, data_points, cluster_indices, speeds)
    stability = zeros(length(candidate_CHs), 1);
    for i = 1:length(candidate_CHs)
        ch_index = candidate_CHs(i);
        ch_coordinates = data_points(ch_index, :);
        cluster_id = cluster_indices(ch_index);
        
        % Get members of the same cluster
        cluster_members = find(cluster_indices == cluster_id);
        cluster_points = data_points(cluster_members, :);
        
        % Calculate distance to all other cluster members
        distances = sqrt(sum((cluster_points - ch_coordinates).^2, 2));
        
        % Calculate stability as the inverse of the average distance
        avg_distance = mean(distances);
        speed_influence = max(0, 1 - 0.01 * speeds(ch_index)); % Include speed influence
        
        % Inverse to ensure higher is better
        stability(i) = (1 / (1 + avg_distance)) * speed_influence;
    end
end
function [optimized_centroids, optimized_clusters] = gwo_centroid_optimization(data_points, num_clusters)
    % Initialize centroids randomly by selecting num_clusters data points from the data
    centroids = data_points(randperm(size(data_points, 1), num_clusters), :);

    % Number of GWO iterations
    num_iterations = 100;

    % Store the best centroids and the corresponding clusters
    best_centroids = centroids;  % Store the best centroids from all iterations
    best_fitness = inf;  % Initialize the best fitness value (lowest error)

    % Iterate through GWO iterations
    for iter = 1:num_iterations
        % Perform k-means clustering with the current centroids as the starting points
        [idx_kmeans, ~] = kmeans(data_points, num_clusters, 'Start', centroids);

        % Calculate the sum of squared errors (fitness function)
        fitness = sum(sum((data_points - centroids(idx_kmeans, :)).^2));  % Sum of squared errors

        % If the current fitness is better (lower error), update the best centroids
        if fitness < best_fitness
            best_fitness = fitness;
            best_centroids = centroids;  % Store the best centroids
        end
        
        % For the next iteration, update the centroids
        % Calculate new centroids based on the current cluster assignments
        new_centroids = zeros(num_clusters, size(data_points, 2));
        for c = 1:num_clusters
            % Calculate the mean of points assigned to each cluster
            new_centroids(c, :) = mean(data_points(idx_kmeans == c, :), 1);
        end

        % Update centroids for the next iteration
        centroids = new_centroids;
    end
    
    % Return the best centroids after optimization
    optimized_centroids = best_centroids;

    % Perform K-means clustering with the optimized centroids to get the final clusters
    [optimized_clusters, ~] = kmeans(data_points, num_clusters, 'Start', optimized_centroids);
end
function [data, data_points] = parse_fcd_xml_to_df(filename)
    % Initialize storage for positions
    x = [];
    y = [];
    speed = [];
    
    % Read XML file
    xmlData = xmlread(filename);
    timesteps = xmlData.getElementsByTagName('timestep');
    
    % Loop through each timestep
    for t = 0:timesteps.getLength()-1
        timestep = timesteps.item(t);
        vehicles = timestep.getElementsByTagName('vehicle');
        
        % Loop through each vehicle in the timestep
        for v = 0:vehicles.getLength()-1
            vehicle = vehicles.item(v);
            
            % Get x, y, and speed attributes
            x_attr = str2double(vehicle.getAttribute('x'));
            y_attr = str2double(vehicle.getAttribute('y'));
            speed_attr = str2double(vehicle.getAttribute('speed'));
            
            % Check if all attributes are present
            if ~isnan(x_attr) && ~isnan(y_attr) && ~isnan(speed_attr)
                x = [x; x_attr];  % Append x-coordinate
                y = [y; y_attr];  % Append y-coordinate
                speed = [speed; speed_attr];  % Append speed
            end
        end
    end
    
    % Create the data matrix (x, y, speed)
    data = [x, y, speed];
    
    % Create the data_points matrix (x, y)
    data_points = [x, y];
end
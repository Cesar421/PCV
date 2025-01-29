% Load stereo images
left = im2double(imread('left.png'));
right = im2double(imread('right.png'));

% Display input images
figure; imshow(left, []); title('Left Image');
figure; imshow(right, []); title('Right Image');

% Optimal parameters (set manually or from optimization in Task (c))
window_radius = 5; % Example optimal window size
search_range = 50; % Example optimal search range
[rows, cols] = size(left);

% Initialize disparity map
disparity_map = NaN(rows, cols);

% Precompute mean values for faster NCC (Task a)
left_means = conv2(left, ones(2*window_radius+1) / (2*window_radius+1)^2, 'same');
right_means = conv2(right, ones(2*window_radius+1) / (2*window_radius+1)^2, 'same');

% Custom NCC function
ncc = @(a, b, ma, mb) sum((a(:) - ma) .* (b(:) - mb)) / ...
                      sqrt(sum((a(:) - ma).^2) * sum((b(:) - mb).^2));

% Compute disparity map (Task a)
for i = 1+window_radius:rows-window_radius
    for j = 1+window_radius:cols-window_radius
        % Extract reference window
        ref_window = left(i-window_radius:i+window_radius, j-window_radius:j+window_radius);
        ref_mean = left_means(i, j);

        % Skip homogeneous regions
        if var(ref_window(:)) == 0
            continue;
        end

        best_ncc = -inf;
        best_disp = 0;

        % Search within the defined range
        for d = -search_range:search_range
            if j+d-window_radius > 0 && j+d+window_radius <= cols
                search_window = right(i-window_radius:i+window_radius, j+d-window_radius:j+d+window_radius);
                search_mean = right_means(i, j+d);
                current_ncc = ncc(ref_window, search_window, ref_mean, search_mean);

                if current_ncc > best_ncc
                    best_ncc = current_ncc;
                    best_disp = d;
                end
            end
        end

        % Record disparity value
        disparity_map(i, j) = best_disp;
    end
end

% Normalize disparity map for visualization
valid_disparities = ~isnan(disparity_map);
min_disp = min(disparity_map(valid_disparities));
max_disp = max(disparity_map(valid_disparities));
disparity_map = (disparity_map - min_disp) / (max_disp - min_disp) * 255;

% Display disparity map
figure; imshow(uint8(disparity_map)); title('Disparity Map');

% Apply median filtering to smooth disparity map
smoothed_disparity_map = medfilt2(disparity_map, [5 5]);
figure; imshow(uint8(smoothed_disparity_map)); title('Smoothed Disparity Map');

% ---- Depth Map Calculation ----
% Parameters for depth calculation
f = 1; % Example focal length (in pixels)
B = 0.1; % Baseline (in meters)

% Ensure disparity values are valid (avoid zero or NaN values)
disparity_map(disparity_map == 0) = NaN; % Handle zero disparities
disparity_map(disparity_map < 1) = NaN; % Avoid very low disparities for stability

% Compute depth map
depth_map = (f * B) ./ abs(disparity_map);

% Replace invalid depths for visualization
depth_map(isnan(depth_map)) = 0;

% Normalize depth map to [0, 255] for visualization
valid_depths = depth_map > 0; % Exclude zeros from normalization
min_depth = min(depth_map(valid_depths));
max_depth = max(depth_map(valid_depths));
depth_map = (depth_map - min_depth) / (max_depth - min_depth) * 255;

% Display depth map
figure; imshow(uint8(depth_map)); title('Depth Map');

% Apply median filtering to smooth depth map
smoothed_depth_map = medfilt2(depth_map, [5 5]);
figure; imshow(uint8(smoothed_depth_map)); title('Smoothed Depth Map');
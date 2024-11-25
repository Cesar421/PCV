% CAMERA CALIBRATION USING DIRECT LINEAR TRANSFORMATION (DLT)

% Step 1: Load and Display the Image
img = imread('a1.jpg'); % Replace with your image file name 
imshow(img);
title('Calibration Object');
pause; % Pause to allow viewing

% Step 2: Manual Selection of Image Points
% Prompt the user to select corresponding points manually
disp('Click on the image to select calibration points (at least 6 points).');
[x, y] = ginput(6); % Adjust the number of points based on your setup
image_points = [x, y]; % Save the manually selected 2D points
hold on;
plot(x, y, 'r+', 'MarkerSize', 10); % Visualize selected points
disp('Selected Image Points:');
disp(image_points);

% Step 3: Define Object Points
% Define 3D object points in the calibration object's coordinate system
object_points = [
    100, 200, 0;  % Bottom-left corner
    300, 200, 0;  % Bottom-right corner
    300, 0, 200;  % Top-left corner
    100, 0, 200;  % Top-right corner
    0, 200, 200;  % Top-left corner of the back face
    0, 200, 100;  % Top-right corner of the back face
];

disp('Object Points:');
disp(object_points);

% Step 4: Compute the Projection Matrix using DLT
function P = computeDLT(object_points, image_points)
    num_points = size(object_points, 1);
    A = zeros(2 * num_points, 12);
    
    for i = 1:num_points
        X = object_points(i, :);
        x = image_points(i, 1);
        y = image_points(i, 2);
        
        A(2 * i - 1, :) = [-X, -1, 0, 0, 0, 0, x * X, x];
        A(2 * i, :) = [0, 0, 0, 0, -X, -1, y * X, y];
    end
    
    % Solve using SVD
    [~, ~, V] = svd(A);
    P = reshape(V(:, end), 4, 3)';
end

P = computeDLT(object_points, image_points);
disp('Projection Matrix (P):');
disp(P);

% Step 5: Decompose the Projection Matrix
function [K, R, C] = decomposeProjectionMatrix(P)
    % Decompose the projection matrix into K (intrinsic), R (rotation), and C (camera center)
    M = P(:, 1:3);
    [K, R] = rq(M); % Perform RQ decomposition
    
    % Normalize intrinsic matrix
    K = K ./ K(3, 3);
    
    % Compute camera center
    C = -inv(M) * P(:, 4);
end

function [R, Q] = rq(A)
    % RQ decomposition using QR decomposition
    [Q, R] = qr(flipud(A)');
    R = flipud(R');
    R = fliplr(R);
    Q = Q';
    Q = flipud(Q);
end

[K, R, C] = decomposeProjectionMatrix(P);
disp('Intrinsic Matrix (K):');
disp(K);
disp('Rotation Matrix (R):');
disp(R);
disp('Camera Center (C):');
disp(C);

% Step 6: Compute Reprojection Errors
function residuals = computeResiduals(P, object_points, image_points)
    num_points = size(object_points, 1);
    residuals = zeros(num_points, 1);
    
    for i = 1:num_points
        X = [object_points(i, :), 1]'; % Homogeneous coordinates
        x_proj = P * X;               % Project point
        x_proj = x_proj ./ x_proj(3); % Normalize
        
        residuals(i) = norm(image_points(i, :) - x_proj(1:2)');
    end
end

residuals = computeResiduals(P, object_points, image_points);
disp('Reprojection Residuals:');
disp(residuals);

% Step 7: Visualize Reprojected Points
% Reproject points
num_points = size(object_points, 1);
reprojected_points = zeros(num_points, 2);

for i = 1:num_points
    X = [object_points(i, :), 1]'; % Homogeneous coordinates
    x_proj = P * X;
    x_proj = x_proj ./ x_proj(3); % Normalize
    reprojected_points(i, :) = x_proj(1:2)';
end

% Plot original and reprojected points
imshow(img); hold on;
plot(image_points(:, 1), image_points(:, 2), 'ro', 'MarkerSize', 10, 'DisplayName', 'Original Points');
plot(reprojected_points(:, 1), reprojected_points(:, 2), 'bx', 'MarkerSize', 10, 'DisplayName', 'Reprojected Points');
legend;
title('Reprojection of Calibration Points');
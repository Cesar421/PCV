% Camera Calibration using Direct Linear Transformation (DLT)
clc;
clear all;

%% Step 1: Load and Display Image
img = imread('Object3.JPG'); % Replace with your calibration image file
imgR = imread('point_selection_reference.JPG');
figure;
subplot(1, 2, 1); % 1 row, 2 columns, first position
imshow(img);
title('Image to calibrate');
subplot(1, 2, 2); % 1 row, 2 columns, first position
imshow(imgR);
title('Reference image real points');

%% Step 2: Select Image Points (Manual Selection)
disp('Click on the image to select at least 6 corresponding points.');
[x, y] = ginput(6); % Allow user to manually select 6 points
image_points = [x, y]; % Save 2D image points
hold on;
plot(x, y, 'r+', 'MarkerSize', 10); % Display selected points on the image
disp('Selected Image Points (2D):');
disp(image_points);

%% Step 3: Define Object Points with \( (0, 0, 0) \) Origin
% Updated object points to match the first code, see the image
object_points = [
    4,  70,  0;  %point1
    28, 46,  0;  %point2
    0,  87, 70;  %point3
    28,  0, 70;  %point4
    0,   4, 87;  %point5
    0,  87, 87;  %point6
];

%% Step 4: Compute Projection Matrix using DLT
function P = computeDLT(object_points, image_points)
    % Ensure homogeneous coordinates
    num_points = size(object_points, 1);
    A = zeros(2 * num_points, 12); % Design matrix for linear equations
    for i = 1:num_points
        X = object_points(i, :); % 3D object point
        x = image_points(i, 1); % Corresponding x-coordinate in image
        y = image_points(i, 2); % Corresponding y-coordinate in image
        
        % Construct the matrix rows
        A(2 * i - 1, :) = [-X, -1, 0, 0, 0, 0, x * X, x];
        A(2 * i, :) = [0, 0, 0, 0, -X, -1, y * X, y];
    end
    
    % Solve using SVD
    [~, ~, V] = svd(A);
    P = reshape(V(:, end), 4, 3)'; % Reshape last column of V into 3x4 matrix
end

P = computeDLT(object_points, image_points);
disp('Projection Matrix (P):');
disp(P);

%% Step 5: Decompose Projection Matrix using RQ Decomposition
function [K, R, C] = decomposeProjectionMatrix(P)
    M = P(:, 1:3); % Extract 3x3 submatrix
    [K, R] = rq(M); % Perform RQ decomposition
    
    % Normalize intrinsic matrix
    K = K ./ K(3, 3); % Scale K to make K(3,3) = 1
    
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

%% Step 6: Calculate Rotation Angles (omega, phi, kappa)
% Rotation angles from the rotation matrix R
omega = atan2(R(3,2), R(3,3)); % Roll angle (omega)
phi = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2)); % Pitch angle (phi)
kappa = atan2(R(2,1), R(1,1)); % Yaw angle (kappa)

% Convert to degrees
omega_deg = omega * 180 / pi;
phi_deg = phi * 180 / pi;
kappa_deg = kappa * 180 / pi;

disp('Rotation Angles (degrees):');
disp(['Omega (Roll): ', num2str(omega_deg)]);
disp(['Phi (Pitch): ', num2str(phi_deg)]);
disp(['Kappa (Yaw): ', num2str(kappa_deg)]);

%% Step 7: Extract and Display Intrinsic Parameters from K
% Principal distance (focal length) a_x
a_x = K(1, 1);
% Skew factor s
s = K(1, 2);
% Principal point (x_0, y_0)
x_0 = K(1, 3);
y_0 = K(2, 3);
% Aspect ratio a_y / a_x
aspect_ratio = K(2, 2) / K(1, 1);

% Display results
disp('Intrinsic Parameters:');
disp(['Principal Distance (a_x): ', num2str(a_x)]);
disp(['Skew Factor (s): ', num2str(s)]);
disp(['Principal Point (x_0, y_0): (', num2str(x_0), ', ', num2str(y_0), ')']);
disp(['Aspect Ratio (a_y / a_x): ', num2str(aspect_ratio)]);

%% Step 8: Compute Reprojection Errors
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

%% Step 9: Visualize Original vs Reprojected Points
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


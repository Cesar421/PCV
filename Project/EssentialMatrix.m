% Load the data from calib_points.dat
data = load('calib_points.dat');
x1 = data(:, 1:2); % Image coordinates in camera 1
x2 = data(:, 3:4); % Image coordinates in camera 2
X = data(:, 5:7);  % 3D object points

% Convert image coordinates to homogeneous coordinates
x1_h = [x1, ones(size(x1, 1), 1)]; % Homogeneous coordinates for x1
x2_h = [x2, ones(size(x2, 1), 1)]; % Homogeneous coordinates for x2

% Compute the calibration matrices K1 and K2 using the DLT method
K1 = computeCalibrationMatrix(x1_h, X);
K2 = computeCalibrationMatrix(x2_h, X);

% Debug: Check for invalid values in K1 and K2
if any(isnan(K1(:))) || any(isinf(K1(:)))
    error('Calibration matrix K1 contains invalid values (NaN or Inf).');
end
if any(isnan(K2(:))) || any(isinf(K2(:)))
    error('Calibration matrix K2 contains invalid values (NaN or Inf).');
end

% Display calibration matrices
disp('Calibration matrix K1:');
disp(K1);

disp('Calibration matrix K2:');
disp(K2);

% Normalize the image coordinates using the calibration matrices
x1_norm = (K1 \ x1_h')'; % Normalized homogeneous coordinates for x1
x2_norm = (K2 \ x2_h')'; % Normalized homogeneous coordinates for x2

% Ensure normalized coordinates are in homogeneous form (3D)
x1_norm = x1_norm(:, 1:3); % Keep only x, y, and w (homogeneous coordinate)
x2_norm = x2_norm(:, 1:3); % Keep only x, y, and w (homogeneous coordinate)

% Estimate the essential matrix E using the 8-point algorithm
E = estimateEssentialMatrix(x1_norm, x2_norm);

disp('Essential matrix E:');
disp(E);

% Resolve the fourfold ambiguity of the essential matrix
[E_correct, R, t] = resolveFourfoldAmbiguity(E, x1_norm, x2_norm);

disp('Corrected essential matrix E_correct:');
disp(E_correct);

disp('Rotation matrix R:');
disp(R);

disp('Translation vector t:');
disp(t);

% Verify orthogonality of the rotation matrix
disp('Check orthogonality of R (R^T * R should be identity):');
disp(R' * R);

% Compute the epipolar lines
epipolarLines1 = computeEpipolarLines(E_correct, x1_norm);
epipolarLines2 = computeEpipolarLines(E_correct', x2_norm);

% Normalize epipolar lines for visualization
epipolarLines1 = epipolarLines1 ./ epipolarLines1(:, 3); % Normalize by third component
epipolarLines2 = epipolarLines2 ./ epipolarLines2(:, 3); % Normalize by third component

disp('Epipolar lines for camera 1:');
disp(epipolarLines1);

disp('Epipolar lines for camera 2:');
disp(epipolarLines2);

% Visualize the epipolar lines and corresponding image points
visualizeEpipolarLines(x1, x2, epipolarLines1, epipolarLines2);

% Helper functions
function K = computeCalibrationMatrix(x_h, X)
    % Direct Linear Transform (DLT) to compute the calibration matrix
    A = [];
    for i = 1:size(x_h, 1)
        X_i = [X(i, :), 1];
        x_i = x_h(i, 1);
        y_i = x_h(i, 2);
        A = [A; X_i, zeros(1, 4), -x_i * X_i; zeros(1, 4), X_i, -y_i * X_i];
    end
    [~, ~, V] = svd(A);
    K = reshape(V(:, end), 3, 4);
    K = K / K(3, 3); % Normalize the matrix
end

function E = estimateEssentialMatrix(x1, x2)
    % 8-point algorithm to estimate the essential matrix
    A = [];
    for i = 1:size(x1, 1)
        x1_i = x1(i, :);
        x2_i = x2(i, :);
        A = [A; x2_i(1) * x1_i, x2_i(2) * x1_i, x2_i(3) * x1_i];
    end
    [~, ~, V] = svd(A);
    E = reshape(V(:, end), 3, 3);

    % Enforce the rank-2 constraint
    [U, S, V] = svd(E);
    S(3, 3) = 0;
    E = U * S * V';
end

function [E_correct, R, t] = resolveFourfoldAmbiguity(E, x1, x2)
    % Resolve the fourfold ambiguity of the essential matrix
    [U, ~, V] = svd(E);
    W = [0 -1 0; 1 0 0; 0 0 1];

    % Four possible solutions
    R1 = U * W * V';
    R2 = U * W' * V';
    t1 = U(:, 3);
    t2 = -U(:, 3);

    % Check which solution is geometrically plausible
    % (Add logic to select the correct solution based on 3D point triangulation)
    E_correct = E;
    R = R1;
    t = t1;
end

function epipolarLines = computeEpipolarLines(E, x)
    % Compute epipolar lines for the given essential matrix and points
    epipolarLines = E * x';
    epipolarLines = epipolarLines';
end

function visualizeEpipolarLines(x1, x2, epipolarLines1, epipolarLines2)
    % Visualize the epipolar lines and corresponding image points
    figure;
    subplot(1, 2, 1);
    imshow(zeros(500, 500)); % Placeholder for image 1
    hold on;
    plot(x1(:, 1), x1(:, 2), 'r*');
    for i = 1:size(epipolarLines1, 1)
        line = epipolarLines1(i, :);
        plot([0, 500], [-(line(1) * 0 + line(3)) / line(2), -(line(1) * 500 + line(3)) / line(2)], 'b-');
    end
    title('Epipolar Lines in Camera 1');

    subplot(1, 2, 2);
    imshow(zeros(500, 500)); % Placeholder for image 2
    hold on;
    plot(x2(:, 1), x2(:, 2), 'g*');
    for i = 1:size(epipolarLines2, 1)
        line = epipolarLines2(i, :);
        plot([0, 500], [-(line(1) * 0 + line(3)) / line(2), -(line(1) * 500 + line(3)) / line(2)], 'b-');
    end
    title('Epipolar Lines in Camera 2');
end

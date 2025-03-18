function EssentialMatrix2()
    % Task 1a: Get Image Points
    [X1, X2] = GetImageCoordinates();
    disp("Image Points (X1):");
    disp(X1);
    disp("Image Points (X2):");
    disp(X2);

    % Task 1b: Compute Fundamental Matrix
    F = GetFundamentalMatrix(X1, X2);
    disp("Fundamental Matrix (F):");
    disp(F);

    % Task 1c: Compute Projection Matrices
    [P1, P2] = CreateProjectionMatrix(F);
    disp("Projection Matrix (P1):");
    disp(P1);
    disp("Projection Matrix (P2):");
    disp(P2);

    % Task 1d: Compute Calibration Matrices
    K1 = CalibrationMatrix(P1);
    K2 = CalibrationMatrix(P2);
    disp("Calibration Matrix (K1):");
    disp(K1);
    disp("Calibration Matrix (K2):");
    disp(K2);

    % Compute the Essential Matrix
    E = ComputeEssentialMatrix(F, K1, K2);
    disp("Essential Matrix (E):");
    disp(E);

    % Compute Epipolar Lines
    [l1, l2] = ComputeEpipolarLines(E, X1, X2);

    % Generate white-background images
    img1 = GenerateWhiteImage();
    img2 = GenerateWhiteImage();

    % Plot Epipolar Lines
    PlotEpipolarLines(X1, X2, l1, l2, img1, img2);

    % ==========================
    % Task 2: Non-Linear Optimization
    % ==========================

    % Compute Initial Geometric Error
    geometricError_initial = ComputeGeometricError(X1, X2, l2);
    disp("Initial Geometric Error:");
    disp(geometricError_initial);

    % Perform Non-Linear Optimization
    [E_optimized, error_optimized] = NonLinearOptimization(X1, X2, F, K1, K2);
    disp("Optimized Essential Matrix:");
    disp(E_optimized);

    % Compute Final Geometric Error
    geometricError_final = ComputeGeometricError(X1, X2, ComputeEpipolarLines(E_optimized, X1, X2));
    disp("Final Geometric Error after optimization:");
    disp(geometricError_final);
end

% Function to Read Image Points
function [X1, X2] = GetImageCoordinates()
    fh = fopen("calib_points.dat", "r");
    A = fscanf(fh, '%f', [4, inf]);
    fclose(fh);
    X1 = [A(1:2, :); ones(1, size(A, 2))];
    X2 = [A(3:4, :); ones(1, size(A, 2))];
end

% Compute Calibration Matrix
function K = CalibrationMatrix(P)
    M = P(:, 1:3);  
    if abs(det(M)) < 1e-10
        M = M + eye(3) * 1e-10; 
    end
    [Q, R] = qr(inv(M));
    K = inv(R);  
    K = K / K(3,3);
end

% Compute Fundamental Matrix
function fundamentalMatrix = GetFundamentalMatrix(X1, X2)
    [X1_normalized, T1] = NormalizePoints(X1);
    [X2_normalized, T2] = NormalizePoints(X2);
    A = zeros(size(X1_normalized, 2), 9);
    for i = 1:size(X1_normalized, 2)
        A(i, :) = [X1_normalized(1, i) * X2_normalized(1, i), ...
                   X1_normalized(2, i) * X2_normalized(1, i), ...
                   X2_normalized(1, i), ...
                   X1_normalized(1, i) * X2_normalized(2, i), ...
                   X1_normalized(2, i) * X2_normalized(2, i), ...
                   X2_normalized(2, i), ...
                   X1_normalized(1, i), ...
                   X1_normalized(2, i), 1];
    end
    [~, ~, V] = svd(A);
    F = reshape(V(:, end), 3, 3)';
    [U, D, V] = svd(F);
    D(3,3) = 0;
    fundamentalMatrix = U * D * V';
    fundamentalMatrix = T2' * fundamentalMatrix * T1;
end

% Compute Essential Matrix
function E = ComputeEssentialMatrix(F, K1, K2)
    E = K2' * F * K1;
    [U, D, V] = svd(E);
    D(3,3) = 0;
    E = U * D * V';
end

% Compute Projection Matrices
function [P1, P2] = CreateProjectionMatrix(F)
    P1 = [eye(3), zeros(3,1)];
    [U, ~, ~] = svd(F);
    e2 = U(:, end);
    e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
    P2 = [e2x * F, e2];
end

% Compute Epipolar Lines
function [l1, l2] = ComputeEpipolarLines(E, X1, X2)
    l1 = (E * X1);
    l2 = (E' * X2);
end

% Compute Geometric Error
function error = ComputeGeometricError(X1, X2, epipolarLines)
    n = size(X2, 2);
    error = 0;
    for i = 1:n
        l = epipolarLines(:, i);
        x = X2(1, i);
        y = X2(2, i);
        error = error + abs(l(1) * x + l(2) * y + l(3)) / sqrt(l(1)^2 + l(2)^2);
    end
    error = error / n;
end

% Perform Non-Linear Optimization
function [E_optimized, error_optimized] = NonLinearOptimization(X1, X2, F, K1, K2)
    objectiveFunc = @(x) ComputeGeometricError(X1, X2, ComputeEpipolarLines(reshape(x,3,3), X1, X2));
    x0 = reshape(F, 1, []);
    options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
    [x_optimized, error_optimized] = fminunc(objectiveFunc, x0, options);
    E_optimized = reshape(x_optimized, 3, 3);
end

% Generate White Background Image
function img = GenerateWhiteImage()
    img = ones(500, 500, 3);  
end

% Normalize Image Points
function [normalizedPoints, T] = NormalizePoints(points)
    meanX = mean(points(1,:));
    meanY = mean(points(2,:));
    stdX = std(points(1,:));
    stdY = std(points(2,:));

    T = [1/stdX, 0, -meanX/stdX; 0, 1/stdY, -meanY/stdY; 0, 0, 1];
    normalizedPoints = T * points;
end

% Draw a Single Epipolar Line
function DrawEpipolarLine(l, imgSize)
    x = [1 imgSize(2)];
    y = (-l(1)*x - l(3)) / l(2);
    plot(x, y, 'g', 'LineWidth', 1.5);
end

% Plot Epipolar Lines
function PlotEpipolarLines(X1, X2, l1, l2, img1, img2)
    figure;
    subplot(1,2,1);
    imshow(img1); hold on;
    scatter(X1(1,:), X1(2,:), 'ro');
    for i = 1:size(l1,2)
        DrawEpipolarLine(l1(:,i), size(img1));
    end
    title('Epipolar Lines in Image 1');

    subplot(1,2,2);
    imshow(img2); hold on;
    scatter(X2(1,:), X2(2,:), 'bo');
    for i = 1:size(l2,2)
        DrawEpipolarLine(l2(:,i), size(img2));
    end
    title('Epipolar Lines in Image 2');
end
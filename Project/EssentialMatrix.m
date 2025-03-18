function EssentialMatrix()
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

    % Generate synthetic images
    img1 = GenerateSyntheticImage();
    img2 = GenerateSyntheticImage();

    % Plot Epipolar Lines
    PlotEpipolarLines(X1, X2, l1, l2, img1, img2);
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
    M = P(:, 1:3);  % Extract the first 3 columns of P
    
    % Regularization: Add a small value to prevent singularity
    if abs(det(M)) < 1e-10
        M = M + eye(3) * 1e-10; 
    end

    % Compute QR decomposition for stability
    [Q, R] = qr(inv(M));
    K = inv(R);  

    % Normalize K so that K(3,3) = 1
    K = K / K(3,3);
end

% Compute the Fundamental Matrix
function fundamentalMatrix = GetFundamentalMatrix(X1, X2)
    % Normalize image coordinates
    [X1_normalized, T1] = NormalizePoints(X1);
    [X2_normalized, T2] = NormalizePoints(X2);

    % Construct matrix A for the 8-point algorithm
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

    % Compute the Fundamental Matrix using SVD
    [~, ~, V] = svd(A);
    F = reshape(V(:, end), 3, 3)';

    % Enforce rank-2 constraint on F
    [U, D, V] = svd(F);
    D(3,3) = 0;
    fundamentalMatrix = U * D * V';

    % Denormalize F
    fundamentalMatrix = T2' * fundamentalMatrix * T1;
end

% Compute Essential Matrix
function E = ComputeEssentialMatrix(F, K1, K2)
    % Compute Essential Matrix
    E = K2' * F * K1;
    
    % Enforce Singularity Constraint (Rank 2 Constraint)
    [U, D, V] = svd(E);
    D(3,3) = 0; % Set smallest singular value to zero
    E = U * D * V';
end

% Compute Projection Matrices
function [P1, P2] = CreateProjectionMatrix(F)
    % Projection Matrix for first camera (Identity + small perturbation)
    P1 = [eye(3) + randn(3) * 1e-4, zeros(3,1)];

    % Compute epipole for second camera
    [U, ~, ~] = svd(F);
    e2 = U(:, end);

    % Compute skew-symmetric matrix
    e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];

    % Compute second camera projection matrix
    P2 = [e2x * F, e2];
end

% Compute Epipolar Lines
function [l1, l2] = ComputeEpipolarLines(E, X1, X2)
    l1 = (E * X1); % Epipolar lines in Image 1
    l2 = (E' * X2); % Epipolar lines in Image 2
end

% Plot Epipolar Lines
function PlotEpipolarLines(X1, X2, l1, l2, img1, img2)
    figure;

    % Plot Image 1 with Epipolar Lines
    subplot(1,2,1);
    imshow(img1); hold on;
    scatter(X1(1,:), X1(2,:), 'ro'); % Plot points
    for i = 1:size(l1,2)
        DrawEpipolarLine(l1(:,i), size(img1)); % Plot lines
    end
    title('Epipolar Lines in Image 1');

    % Plot Image 2 with Epipolar Lines
    subplot(1,2,2);
    imshow(img2); hold on;
    scatter(X2(1,:), X2(2,:), 'bo'); % Plot points
    for i = 1:size(l2,2)
        DrawEpipolarLine(l2(:,i), size(img2)); % Plot lines
    end
    title('Epipolar Lines in Image 2');
end

% Draw a Single Epipolar Line
function DrawEpipolarLine(l, imgSize)
    x = [1 imgSize(2)];
    y = (-l(1)*x - l(3)) / l(2);
    plot(x, y, 'g', 'LineWidth', 1.5);
end

% Generate Synthetic Image (Black Background)
function img = GenerateSyntheticImage()
    img = zeros(500, 500, 3);  % Black background of 500x500 pixels
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
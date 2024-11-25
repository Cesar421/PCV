% Corresponding 3D object points (homogeneous coordinates)
X = [44.7  -142.4  258.7  1;
    -103.6 -146.6  154.4  1;
     47.4  -150.1   59.8  1;
    -152.2   59.4  245.2  1;
    -153.3  -96.9  151.3  1;
    -149.4   52.7   46.9  1];

% Corresponding 2D image points
x = [18.5  46.8;
     99.1 146.5;
     13.8 221.8;
    242.1  52.5;
    151.1 147.1;
    243.1 224.5];

% Normalize 3D and 2D points
% Normalize 3D points
mean_X = mean(X(:, 1:3), 1);
std_X = std(X(:, 1:3), 0, 1);
T_X = [diag(1./std_X), -mean_X' ./ std_X'; 0 0 0 1];
X_normalized = (T_X * X')';

% Normalize 2D points
mean_x = mean(x, 1);
std_x = std(x, 0, 1);
T_x = [diag(1./std_x), -mean_x' ./ std_x'; 0 0 1];
x_h = [x, ones(size(x, 1), 1)];
x_normalized = (T_x * x_h')';

% Construct the design matrix A
n = size(X, 1);
A = zeros(2*n, 12);
for i = 1:n
    X_i = X_normalized(i, :);
    x_i = x_normalized(i, 1);
    y_i = x_normalized(i, 2);
    A(2*i-1, :) = [X_i, zeros(1, 4), -x_i * X_i];
    A(2*i, :)   = [zeros(1, 4), X_i, -y_i * X_i];
end

% Solve using SVD
[~, ~, V] = svd(A);
P_normalized = reshape(V(:, end), 4, 3)';

% Denormalize the projection matrix
P = inv(T_x) * P_normalized * T_X;

% Normalize P (scale so P(3,4) = 1)
P = P / P(3, 4);

% Display the projection matrix
disp('Projection Matrix P:');
disp(P)

% Scale P to ensure P(3,4) = 1
P = P / P(3, 4);

% 1. Extract M (3x3) from P
M = P(1:3, 1:3);

% 2. RQ Decomposition
[R, K] = rq(M);

% Ensure positive diagonal for K
for i = 1:3
    if K(i, i) < 0
        K(:, i) = -K(:, i);
        R(i, :) = -R(i, :);
    end
end

% 3. Intrinsic Parameters
principal_distance = K(1, 1); % α_x
skew = K(1, 2);              % Skew
principal_point = [K(1, 3), K(2, 3)]; % (x0, y0)
aspect_ratio = K(2, 2) / K(1, 1);     % α_y / α_x

% Display Intrinsic Parameters
disp('Intrinsic Parameters (K):');
disp(K);
disp('Principal Distance (α_x):');
disp(principal_distance);
disp('Skew:');
disp(skew);
disp('Principal Point (x0, y0):');
disp(principal_point);
disp('Aspect Ratio (α_y / α_x):');
disp(aspect_ratio);

% 4. Rotation Matrix and Angles
disp('Rotation Matrix (R):');
disp(R);

omega = atan2(R(3, 2), R(3, 3)) * (180 / pi); % X-axis
phi = asin(-R(3, 1)) * (180 / pi);            % Y-axis
kappa = atan2(R(2, 1), R(1, 1)) * (180 / pi); % Z-axis

disp('Rotation Angles (omega, phi, kappa) [degrees]:');
disp([omega, phi, kappa]);

% 5. Projection Center
[~, ~, V] = svd(P);
C_h = V(:, end);
C = C_h(1:3) / C_h(4);

disp('Projection Center (C):');
disp(C');

function [R, K] = rq(M)
    % Reverse the rows of M
    [Q, R_tilde] = qr(flipud(M)');
    R = flipud(R_tilde'); % Reverse rows again for R
    K = flipud(Q');       % Reverse rows again for K
end
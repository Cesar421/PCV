function Assignment4

clc;
clear all;

% Load images into MATLAB
f = imread("image1.jpg");
g = imread("image2.jpg");

% Select corresponding points on both images
figure(1), imshow(f);
[Image1X, Image1Y] = ginput(8); % User selects 8 points on image 1
hold on;
plot(Image1X, Image1Y); % Overlay the selected points
hold off;
figure(2), imshow(g);
[Image2X, Image2Y] = ginput(8); % User selects 8 points on image 2
hold on;
plot(Image2X, Image2Y); % Overlay the selected points

% Create homogeneous coordinates for both images
Euc = ones(8, 1); % Homogeneous coordinate (z = 1)
Image1 = [Image1X, Image1Y, Euc];
Image2 = [Image2X, Image2Y, Euc];

% Compute the centroid (midpoint) of points in both images
t1 = [mean(Image1X), mean(Image1Y)]; % Midpoint for image 1
disp("The midpoint for Image 1 is:");
disp(t1);
t2 = [mean(Image2X), mean(Image2Y)]; % Midpoint for image 2
disp("The midpoint for Image 2 is:");
disp(t2);

% Translate the points so the centroid is at the origin
X1_T = Image1X - t1(1); % X-coordinates translated for image 1
Y1_T = Image1Y - t1(2); % Y-coordinates translated for image 1
Image1Translated = [X1_T, Y1_T];
X2_T = Image2X - t2(1); % X-coordinates translated for image 2
Y2_T = Image2Y - t2(2); % Y-coordinates translated for image 2
Image2Translated = [X2_T, Y2_T];

% Compute the scaling factor for normalization
s1 = [mean(abs(X1_T)), mean(abs(Y1_T))]; % Scaling for image 1
disp("The scaling for Image 1 is:");
disp(s1);
s2 = [mean(abs(X2_T)), mean(abs(Y2_T))]; % Scaling for image 2
disp("The scaling for Image 2 is:");
disp(s2);

% Construct the transformation matrices
s1_TransMatrix = [1/s1(1), 0, 0; 0, 1/s1(2), 0; 0, 0, 1]; % Scaling matrix for image 1
t1_TransMatrix = [1, 0, -t1(1); 0, 1, -t1(2); 0, 0, 1]; % Translation matrix for image 1
Image1Transformation = s1_TransMatrix * t1_TransMatrix; % Transformation for image 1
disp("The coordinate transformation for Image 1 is:");
disp(Image1Transformation);
s2_TransMatrix = [1/s2(1), 0, 0; 0, 1/s2(2), 0; 0, 0, 1]; % Scaling matrix for image 2
t2_TransMatrix = [1, 0, -t2(1); 0, 1, -t2(2); 0, 0, 1]; % Translation matrix for image 2
Image2Transformation = s2_TransMatrix * t2_TransMatrix; % Transformation for image 2
disp("The coordinate transformation for Image 2 is:");
disp(Image2Transformation);

% Compute conditioned coordinates
Image1Conditioned = [X1_T * Image1Transformation(1, 1), Y1_T * Image1Transformation(2, 2)];
disp("The conditioned coordinates of Image 1 are:");
disp(Image1Conditioned);
Image2Conditioned = [X2_T * Image2Transformation(1, 1), Y2_T * Image2Transformation(2, 2)];
disp("The conditioned coordinates of Image 2 are:");
disp(Image2Conditioned);

% Transpose the conditioned coordinates for easier manipulation
conditioned1 = transpose(Image1Conditioned);
conditioned2 = transpose(Image2Conditioned);

% Construct the design matrix A
A = zeros(size(conditioned1, 2), 9); % Initialize the matrix
for i = 1:size(conditioned1, 2)
    % Populate each row of the design matrix
    A(i, :) = [conditioned1(1, i)*conditioned2(1, i), ...
               conditioned1(2, i)*conditioned2(1, i), ...
               conditioned2(1, i), ...
               conditioned1(1, i)*conditioned2(2, i), ...
               conditioned1(2, i)*conditioned2(2, i), ...
               conditioned2(2, i), ...
               conditioned1(1, i), ...
               conditioned1(2, i), ...
               1];
end
disp("The Design Matrix A for Image 1 & 2 is:");
disp(A);

% Perform singular value decomposition to compute the fundamental matrix
[U, D, V] = svd(A);
disp("The singular values (D) are:");
disp(D);

% Compute the fundamental matrix
fMatrix = reshape(V(:, end), 3, 3)'; % Reshape last column of V into 3x3 matrix
fMatrix = Image2Transformation' * fMatrix * Image1Transformation; % Reverse conditioning
fMatrix = fMatrix / fMatrix(end, end); % Normalize

% Enforce singularity constraint on the fundamental matrix
fMatrix_P = singularity(fMatrix);

% Compute epipolar lines for points in both images
for i = 1:size(Image1, 1)
    l1(i, :) = fMatrix_P' * Image2(i, :)'; % Epipolar line for image 1
    l2(i, :) = fMatrix_P * Image1(i, :)'; % Epipolar line for image 2
end

% Display epipolar lines
disp("Image 1 Epipolar Lines:");
disp(l1);
disp("Image 2 Epipolar Lines:");
disp(l2);

% Draw epipolar lines on images
f = imread('image1.jpg'); 
figure, imshow(f), hold on
for i = 1:size(l1, 1)
    hline(l1(i, :)); % Draw line for image 1
end
f = imread('image2.jpg'); 
figure, imshow(f), hold on
for i = 1:size(l2, 1)
    hline(l2(i, :)); % Draw line for image 2
end

% Compute geometric error for the points
for i = 1:size(Image2, 1)
    a = Image2(i, 1); b = Image2(i, 2);
    x = l2(i, 1); y = l2(i, 2); z = l2(i, 3);
    u = Image1(i, 1); v = Image1(i, 2);
    numerator = (a*x + b*y + z)^2;
    denominator = a^2 + b^2 + u^2 + v^2;
    geometricError = numerator / denominator;
    disp("The geometric error is:");
    disp(geometricError);
end
end

% Function to enforce singularity constraint on the fundamental matrix
function fMatrix_P = singularity(fMatrix)
if det(fMatrix) == 0
    fMatrix_P = fMatrix; % Already satisfies constraint
else
    [U, D, V] = svd(fMatrix);
    D(end, end) = 0; % Enforce singularity
    fMatrix_P = U * D * V';
end
end

% Helper function to draw epipolar lines
function hline(l, varargin)
if abs(l(1)) < abs(l(2)) % Line is more horizontal
    xlim = get(gca, 'XLim');
    x1 = cross(l, [1; 0; -xlim(1)]);
    x2 = cross(l, [1; 0; -xlim(2)]);
else % Line is more vertical
    ylim = get(gca, 'YLim');
    x1 = cross(l, [0; 1; -ylim(1)]);
    x2 = cross(l, [0; 1; -ylim(2)]);
end
x1 = x1 / x1(3);
x2 = x2 / x2(3);
line([x1(1) x2(1)], [x1(2) x2(2)], varargin{:});
end

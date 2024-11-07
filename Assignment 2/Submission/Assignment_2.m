function Assignment_2
    % Load images
    img1 = imread('1_left.jpg');
    img2 = imread('2_central.jpg');
    img3 = imread('3_right.jpg');

    % Display images
    figure('Name', 'Three Images', 'NumberTitle', 'off');
    subplot(1, 3, 1);
    imshow(img1);
    title('Left Image');
    subplot(1, 3, 2);
    imshow(img2);
    title('Central Image');
    subplot(1, 3, 3);
    imshow(img3);
    title('Right Image');

    % Select points for Left <-> Central correspondence
    disp('Select 4 points between the Left and Central images');
    subplot(1, 3, 1);
    [x1, y1] = ginput(4);
    subplot(1, 3, 2);
    [x2, y2] = ginput(4);
    points_left = [x1, y1];
    points_central = [x2, y2];

    % Compute homography H_left_to_central
    H_left_to_central = computeHomography(points_left, points_central);
    disp('Homography matrix H_left_to_central:');
    disp(H_left_to_central);

    % Select points for Central <-> Right correspondence
    disp('Select 4 points between the Central and Right images');
    subplot(1, 3, 2);
    [x3, y3] = ginput(4);
    subplot(1, 3, 3);
    [x4, y4] = ginput(4);
    points_central_2 = [x3, y3];
    points_right = [x4, y4];

    % Compute homography H_central_to_right
    H_central_to_right = computeHomography(points_central_2, points_right);
    disp('Homography matrix H_central_to_right:');
    disp(H_central_to_right);
end

function H = computeHomography(points1, points2)
    % Compute the homography matrix from points1 to points2 using SVD
    % Input:
    %   - points1: Nx2 matrix of (x, y) coordinates in the first image
    %   - points2: Nx2 matrix of (x, y) coordinates in the second image
    % Output:
    %   - H: 3x3 homography matrix

    n = size(points1, 1);  % Number of points
    A = [];

    for i = 1:n
        x = points1(i, 1);
        y = points1(i, 2);
        x_prime = points2(i, 1);
        y_prime = points2(i, 2);

        % Construct matrix A for the equation A * h = 0
        A = [A;
             -x, -y, -1,  0,  0,  0, x * x_prime, y * x_prime, x_prime;
              0,  0,  0, -x, -y, -1, x * y_prime, y * y_prime, y_prime];
    end

    % Perform SVD on A
    [~, ~, V] = svd(A);

    % Homography is the last column of V reshaped as a 3x3 matrix
    H = reshape(V(:, end), 3, 3)';
end

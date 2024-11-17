function Assg_2 %ting
    %% Task 1: Image Acquisition
    % Step 1: Load the three images
    img1 = imread('IMG1.jpg');
    img2 = imread('IMG2.jpg');
    img3 = imread('IMG3.jpg');
    
    % Step 2: Display each image in a separate figure to verify they are loaded correctly
    figure;
    subplot(1, 3, 1); imshow(img1); title('Image 1');
    subplot(1, 3, 2); imshow(img2); title('Image 2');
    subplot(1, 3, 3); imshow(img3); title('Image 3');
    
    %% Task 2: Correspondence Analysis with Error Handling
    try
        % Select corresponding points between Image 1 and Image 2
        fprintf('Select at least 4 corresponding points between Image 1 and Image 2.\n');
        figure; imshow(img1); title('Select points in Image 1');
        [x1, y1] = ginput(4); % Select 4 points in Image 1
        points1 = [x1, y1];
        
        figure; imshow(img2); title('Select corresponding points in Image 2');
        [x2, y2] = ginput(4); % Select 4 corresponding points in Image 2
        points2 = [x2, y2];
        
        % Select corresponding points between Image 2 and Image 3
        fprintf('Select at least 4 corresponding points between Image 2 and Image 3.\n');
        figure; imshow(img2); title('Select points in Image 2');
        [x3, y3] = ginput(4); % Select 4 points in Image 2
        points3 = [x3, y3];
        
        figure; imshow(img3); title('Select corresponding points in Image 3');
        [x4, y4] = ginput(4); % Select 4 corresponding points in Image 3
        points4 = [x4, y4];
        
        % Save points for verification, if needed
        save('Points_Image1_Image2.mat', 'points1', 'points2');
        save('Points_Image2_Image3.mat', 'points3', 'points4');
        
    catch ME
        fprintf('An error occurred during point selection: %s\n', ME.message);
        return;
    end
    
    %% Task 3: Homography Computation
    % Compute homography matrices between the image pairs
    H12 = computeHomography(points1, points2);  % Homography from Image 1 to Image 2
    H23 = computeHomography(points3, points4);  % Homography from Image 2 to Image 3
    
    %% Task 4: Projective Rectification and Stitching
    % Stitch Image 1 and Image 2
    img12 = warpImage(img1, img2, H12);
    
    % Stitch the mosaic (Image 1 + Image 2) with Image 3
    panorama = warpImage(img12, img3, H23);
    
    %% Task 5: Visualization
    % Display the final panoramic image
    figure; imshow(panorama); title('Panoramic Image of IMG1, IMG2, and IMG3');
end

%% Helper Function for Task 3: Compute Homography Matrix
function H = computeHomography(pointsA, pointsB)
    % Normalize points for better accuracy
    [pointsA, T] = normalizePoints(pointsA);
    [pointsB, T1] = normalizePoints(pointsB);
    
    % Construct the design matrix A
    A = [];
    for i = 1:size(pointsA, 1)
        X = pointsA(i, :);
        X_prime = pointsB(i, :);
        A = [A;
            -X(1), -X(2), -1, 0, 0, 0, X_prime(1) * X(1), X_prime(1) * X(2), X_prime(1);
             0, 0, 0, -X(1), -X(2), -1, X_prime(2) * X(1), X_prime(2) * X(2), X_prime(2)];
    end
    
    % Solve for H using SVD
    [~, ~, V] = svd(A);
    H = reshape(V(:, end), 3, 3)';  % Last column of V gives the homography matrix
    
    % De-normalize the homography matrix
    H = inv(T1) * H * T;
end

%% Helper Function: Normalize Points
function [normalizedPoints, T] = normalizePoints(points)
    centroid = mean(points, 1);
    shiftedPoints = points - centroid;
    avgDist = mean(sqrt(sum(shiftedPoints.^2, 2)));
    scale = sqrt(2) / avgDist;
    T = [scale, 0, -scale * centroid(1);
         0, scale, -scale * centroid(2);
         0, 0, 1];
    normalizedPoints = (T * [points, ones(size(points, 1), 1)]')';
    normalizedPoints = normalizedPoints(:, 1:2);
end

%% Helper Function for Task 4: Warp and Combine Two Images
function output = warpImage(img1, img2, H)
    % Warp img1 onto the plane of img2 using the homography H
    tform = projective2d(H');
    outputView = imref2d(size(img2));
    warpedImg1 = imwarp(img1, tform, 'OutputView', outputView);
    
    % Overlay img2 onto warpedImg1 to create the mosaic
    output = max(warpedImg1, img2);
end
function assignment2p4
    % Load images
    img1 = imread('1_left.jpg');
    img2 = imread('2_central.jpg');
    img3 = imread('3_right.jpg');
    
    % Select points for Left <-> Central correspondence
    disp('Select 4 points between the Left and Central images');
    figure;
    subplot(1, 3, 1);
    imshow(img1);
    title('Left Image');
    subplot(1, 3, 2);
    imshow(img2);
    title('Central Image');
    subplot(1, 3, 3);
    imshow(img3);
    title('Right Image');
    
    subplot(1, 3, 1);
    [x1, y1] = ginput(4);
    subplot(1, 3, 2);
    [x2, y2] = ginput(4);
    points_left = [x1, y1];
    points_central = [x2, y2];
    
    % Compute homography H_left_to_central using SVD
    H_left_to_central = computeHomography(points_left, points_central);
    disp('Homography matrix H_left_to_central:');
    disp(H_left_to_central);
    
    % Apply transformation to adapt img1 to img2
    img1_transformed = applyHomography(H_left_to_central, img1, size(img2));
    
    % Select points for Central <-> Right correspondence
    disp('Select 4 points between the Central and Right images');
    subplot(1, 3, 2);
    [x3, y3] = ginput(4);  % Select points for the Central image
    subplot(1, 3, 3);
    [x4, y4] = ginput(4);  % Select points for the Right image
    points_central_2 = [x3, y3];
    points_right = [x4, y4];
    
    % Compute homography H_central_to_right using SVD
    H_central_to_right = computeHomography(points_central_2, points_right);
    disp('Homography matrix H_central_to_right:');
    disp(H_central_to_right);
    
    % Apply transformation to img3 to adapt to the intermediate mosaic
    img3_transformed = applyHomography(H_central_to_right, img3, size(img2));
    
    % Stitch the images together (left, central, and right images)
    stitched_image = stitchImages(img1_transformed, img2, img3_transformed);
    
    % Display the resulting panoramic image
    figure;
    imshow(stitched_image);
    title('Panoramic Image');
end

% Function to compute the homography matrix using SVD
function H = computeHomography(points1, points2)
    % Number of points
    n = size(points1, 1);  
    A = zeros(2 * n, 9);
    
    % Construct matrix A for the equation A * h = 0
    for i = 1:n
        x = points1(i, 1);
        y = points1(i, 2);
        x_prime = points2(i, 1);
        y_prime = points2(i, 2);
        
        A(2*i-1, :) = [-x, -y, -1, 0, 0, 0, x*x_prime, y*x_prime, x_prime];
        A(2*i, :) = [0, 0, 0, -x, -y, -1, x*y_prime, y*y_prime, y_prime];
    end
    
    % Perform SVD on A
    [~, ~, V] = svd(A);
    
    % Homography is the last column of V reshaped as a 3x3 matrix
    H = reshape(V(:, end), 3, 3)';
end

% Function to apply homography to an image
function img_transformed = applyHomography(H, img, output_size)
    [rows, cols, ~] = size(img);
    
    % Create an empty image for the transformed result
    img_transformed = zeros(output_size(1), output_size(2), 3, 'uint8');
    
    % Inverse homography
    H_inv = inv(H);
    
    % Loop through every pixel in the output image
    for i = 1:output_size(1)
        for j = 1:output_size(2)
            % Find the corresponding point in the original image
            coords = H_inv * [j; i; 1];
            x = round(coords(1) / coords(3));
            y = round(coords(2) / coords(3));
            
            % If the coordinates are within the bounds of the original image, copy the pixel
            if x > 0 && x <= cols && y > 0 && y <= rows
                img_transformed(i, j, :) = img(y, x, :);
            end
        end
    end
end

% Function to stitch the images together
function stitched_img = stitchImages(img1, img2, img3)
    % Simple method for stitching images along the horizontal axis
    % You may want to apply advanced blending techniques for better results
    stitched_img = max(img1, img2);  % Combine left and central image
    stitched_img = max(stitched_img, img3);  % Combine with the right image
end

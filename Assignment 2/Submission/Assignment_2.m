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
 
    % Compute homography H_left_to_central
    H_1_2 = computeHomography(img1, img2);
    
    % Join the two images using geokor
    %img_left_central = geokor(H_1_2, img1, img2);
    % Display the joined image
    %figure;
    %imshow(img_left_central);
    %title('Joined image of Left and Central images');

    % Compute homography H_central_to_right
    H_2_3 = computeHomography(img2, img3);
    
    % Join the images again
    %img_central_right = geokor(H_2_3, img2, img3);
    % Display the final joined image (Central and Right images)
    %figure;
    %imshow(img_central_right);
    %title('Joined image of Central and Right images');
end

function H = computeHomography(imag_a, imag_b)
    figure('Name', 'Image Correspondence', 'NumberTitle', 'off');
    subplot(1, 2, 1);
    imshow(imag_a);
    title('First Image');
    subplot(1, 2, 2);
    imshow(imag_b);
    title('Second Image');

    % Select 4 points in the first image (imag_a)
    disp('Select 4 points between the two images');
    subplot(1, 2, 1);
    [x1, y1] = ginput(4);  % points in the first image
    subplot(1, 2, 2);
    [x2, y2] = ginput(4);  % points in the second image

    % Storing the points
    points_a = [x1, y1]; 
    points_b = [x2, y2];

    % Calculate centroids and translate points
    t_a = mean(points_a);
    t_b = mean(points_b);
    Xa = points_a - t_a;  % Shift points to centroid
    Xb = points_b - t_b;

    % Scaling factor
    Sa = mean(abs(Xa));
    Sb = mean(abs(Xb));

    % Create transformation matrices for normalization
    T_a = [1/Sa(1), 0, 0; 0, 1/Sa(2), 0; 0, 0, 1] * [1, 0, -t_a(1); 0, 1, -t_a(2); 0, 0, 1];
    T_b = [1/Sb(1), 0, 0; 0, 1/Sb(2), 0; 0, 0, 1] * [1, 0, -t_b(1); 0, 1, -t_b(2); 0, 0, 1];

    % Homogeneous coordinates
    points_a_hom = [Xa, ones(4,1)];  % Convert to homogeneous coordinates
    points_b_hom = [Xb, ones(4,1)];  % Convert to homogeneous coordinates

    % Apply the transformations
    a = T_a * points_a_hom(1, :)'; 
    b = T_a * points_a_hom(2, :)'; 
    c = T_a * points_a_hom(3, :)'; 
    d = T_a * points_a_hom(4, :)';
    
    a1 = T_b * points_b_hom(1, :)'; 
    b1 = T_b * points_b_hom(2, :)'; 
    c1 = T_b * points_b_hom(3, :)'; 
    d1 = T_b * points_b_hom(4, :)';

    A = zeros(4,1);
    Z = zeros(1,3);

    % Create the design matrix components
    A1 = [-a1(3)*a', Z, a1(1)*a'; Z, -a1(3)*a', a1(2)*a' ];
    A2 = [-b1(3)*b', Z, b1(1)*b'; Z, -b1(3)*b', b1(2)*b' ];
    A3 = [-c1(3)*c', Z, c1(1)*c'; Z, -c1(3)*c', c1(2)*c' ];
    A4 = [-d1(3)*d', Z, d1(1)*d'; Z, -d1(3)*d', d1(2)*d' ];
    % Combine the matrices into a final design matrix
    A = [A1; A2; A3; A4];
    % Singular Value Decomposition
    [U, D, V] = svd(A);

    % Homography for conditioned coordinates (using the last column of V)
    H1 = reshape(V(:, end), 3, 3);

    % Transformation matrix for original coordinates
    H = inv(T_b) * H1 * T_a;
end


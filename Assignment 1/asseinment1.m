% Given points
x = [2; 3; 1]; % Adding 1 for homogeneous coordinates
y = [-4; 5; 1]; % Adding 1 for homogeneous coordinates

% Compute line l as the cross product of x and y
l = cross(x, y);
disp('Connecting line l between points x and y:');
disp(l);

% Translation vector
t = [6; -7];

% Rotation angle in radians
phi = 15 * (pi / 180);

% Scaling factor
lambda = 8;

% Translation: move x and y by vector t
x_trans = [x(1) + t(1); x(2) + t(2); 1];
y_trans = [y(1) + t(1); y(2) + t(2); 1];

% Rotation matrix
R = [cos(phi), -sin(phi); sin(phi), cos(phi)];

% Rotate the translated points
x_rot = R * x_trans(1:2);
y_rot = R * y_trans(1:2);

% Scale the rotated points
x_prime = lambda * x_rot;
y_prime = lambda * y_rot;

disp('Transformed point x_prime:');
disp(x_prime);
disp('Transformed point y_prime:');
disp(y_prime);

% Homogeneous transformation matrix for translation
T = [1, 0, t(1); 0, 1, t(2); 0, 0, 1];

% Homogeneous transformation matrix for rotation
R_hom = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1];

% Homogeneous transformation matrix for scaling
S = [lambda, 0, 0; 0, lambda, 0; 0, 0, 1];

% Combined transformation matrix
transformation = S * R_hom * T;

% Apply transformation to line l
l_prime = transformation' * l;
disp('Transformed line l_prime:');
disp(l_prime);


% Convert x_prime and y_prime back to homogeneous coordinates
x_prime_hom = [x_prime; 1];
y_prime_hom = [y_prime; 1];

% Check if x' and y' lie on the line l'
is_x_on_l = abs(l_prime' * x_prime_hom) < 1e-10;
is_y_on_l = abs(l_prime' * y_prime_hom) < 1e-10;

if is_x_on_l
    disp('Transformed point x_prime lies on transformed line l_prime');
else
    disp('Transformed point x_prime does not lie on transformed line l_prime');
end

if is_y_on_l
    disp('Transformed point y_prime lies on transformed line l_prime');
else
    disp('Transformed point y_prime does not lie on transformed line l_prime');
end
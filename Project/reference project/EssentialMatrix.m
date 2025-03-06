
function EssentialMatrix()

%Task 1: Estimation of Essential Matrix

%Get the Image Cordinate Points
[X1, X2] = GetImageCordinates();

%Obtain the Fundamental Matrix
fundamentalMatrix = GetFundamentalMatrix(X1,X2);

%Obtain Projection Matrices for the two Image Points
[PN, P] = CreateProjectionMatrix(fundamentalMatrix);

%Obtaining the Callibration Matrix for Image Point 1
K1 = CallibrationMatrix(PN);

%Obtaining the Callibration Matrix for Image Point 2
K2 = CallibrationMatrix(P);

%Obtaining Essential Matrix
essentialMatrix = K2' * fundamentalMatrix * K1;


%Enforcing Singularity on the Essential Matrix
singularEssentialMatrix = SingularityConstraint(essentialMatrix);

%Obtaining Object Points
ObjectPoint = GetObjectCoordinate();

%Resolving Ambiguity of the Essential Matrix
[R, t] = ResolveEssentialMatrixAmbiguity(singularEssentialMatrix, X1, X2, K1, K2, ObjectPoint);

%Estimating the Epipolar Lines
[l1, l2] = ComputeEpipolarLines(singularEssentialMatrix, X1, X2);

% Plotting the results
%PlotEpipolarLinesAndPoints(X1, X2, l1, l2);


%Task 2: Non-Linear Optimization

% Calculating the initial geometric error
geometricError_initial = ComputeGeometricError(X1, X2, l2);

% Perform non-linear optimization
 [optimizedEssentialMatrix, error_optimized] = NonLinearOptimization(X1, X2, fundamentalMatrix, K1, K2, ObjectPoint);

% Re-calculate the geometric error
geometricError_final = ComputeGeometricError(X1, X2, ComputeEpipolarLines(optimizedEssentialMatrix, X1, X2));

% Report and comment on the results
disp("Initial Geometric Error:");
disp(geometricError_initial);
disp("Final Geometric Error after optimization:");
disp(geometricError_final);

end


function [optimizedEssentialMatrix, error_optimized] = NonLinearOptimization(X1, X2, fundamentalMatrix, K1, K2, ObjectPoint)

% Define the optimization objective function
objectiveFunc = @(x) ComputeGeometricError(X1, X2, ComputeEpipolarLines(reshape(x,3,3), X1, X2));

% Initial guess for the optimization
x0 = reshape(fundamentalMatrix, 1, []);

% Perform non-linear optimization using fminunc
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
[x_optimized, error_optimized] = fminunc(objectiveFunc, x0, options);

% Reshape the optimized parameters to get the essential matrix
optimizedEssentialMatrix = reshape(x_optimized, 3, 3);

end



%Reading the Image point from the provided file 
function [ImageCoordinate1, ImageCoordinate2] = GetImageCordinates()

fh = fopen("calib_points.dat","r");
A = fscanf(fh, '%f%f%f%f', [7 inf]);
fclose(fh);
ImageCoordinate1 = A(1:2, :);
ImageCoordinate2 = A(3:4, :);
ImageCoordinate1(3,:) = 1;
ImageCoordinate2(3,:) = 1;

disp("Image 1 Coordinates:");
disp(ImageCoordinate1);

disp("Image 2 Coordinates:");
disp(ImageCoordinate2);


end

function objectCoordinate = GetObjectCoordinate()

fh = fopen("calib_points.dat","r");
A = fscanf(fh, '%f%f%f%f', [7 inf]);
fclose(fh);

objectCoordinate = A(5:7, :);

disp("Object Coordinates:");
disp(objectCoordinate);

end


%Estimation of the fundamental matrix
function fundamentalMatrix = GetFundamentalMatrix(ImageCoordinate1, ImageCoordinate2)

%Performing translation, scaling and conditioning

%Translation
Image1Midpoint = [mean(ImageCoordinate1(1,:)); mean(ImageCoordinate1(2,:)); mean(ImageCoordinate1(3,:))];
Image2Midpoint = [mean(ImageCoordinate2(1,:)); mean(ImageCoordinate2(2,:)); mean(ImageCoordinate2(3,:))];

Image1Translated = [(ImageCoordinate1(1,:)-Image1Midpoint(1, 1)); (ImageCoordinate1(2,:)-Image1Midpoint(2, 1)); (ImageCoordinate1(3,:)-Image1Midpoint(3, 1))];
Image2Translated = [(ImageCoordinate2(1,:)-Image2Midpoint(1, 1)); (ImageCoordinate2(2,:)-Image2Midpoint(2, 1)); (ImageCoordinate2(3,:)-Image2Midpoint(3, 1))];

%scaling
Image1Scaled = [mean(abs(Image1Translated(1,:))); mean(abs(Image1Translated(2,:))); mean(abs(Image1Translated(3,:)))];
Image2Scaled = [mean(abs(Image2Translated(1,:))); mean(abs(Image2Translated(2,:))); mean(abs(Image2Translated(3,:)))];

%Translation Matrix
Image1TranslatedMatrix = [1, 0, -Image1Midpoint(1,1); 0, 1, -Image1Midpoint(2,1); 0, 0, 1];
Image2TranslatedMatrix = [1, 0, -Image2Midpoint(1,1); 0, 1, -Image2Midpoint(2,1); 0, 0, 1];

%Scaled Matrix
Image1ScaledMatrix = [1/Image1Scaled(1,1), 0, 0; 0, 1/Image1Scaled(2,1), 0; 0, 0, 1];
Image2ScaledMatrix = [1/Image2Scaled(1,1), 0, 0; 0, 1/Image2Scaled(2,1), 0; 0, 0, 1];

%Image 1 and 2 Transformation Matrix
Image1TransformationMatrix = Image1ScaledMatrix * Image1TranslatedMatrix;
Image2TransformationMatrix = Image2ScaledMatrix * Image2TranslatedMatrix;

%Image 1 and 2 Conditioned Coordinates
Image1Conditioned = [Image1Translated(1,:)*Image1TransformationMatrix(1,1); Image1Translated(2,:)*Image1TransformationMatrix(2,2); Image1Translated(3,:)];
Image2Conditioned = [Image2Translated(1,:)*Image2TransformationMatrix(1,1); Image2Translated(2,:)*Image2TransformationMatrix(2,2); Image2Translated(3,:)];

%Formulating the Design Matrix

A = zeros(size(Image1Conditioned, 2), 9);

for i = 1:size(Image1Conditioned,2)
    
    A(i, :) = [Image1Conditioned(1, i)*Image2Conditioned(1, i), Image1Conditioned(2, i)*Image2Conditioned(1, i), Image2Conditioned(1, i), Image1Conditioned(1, i)*Image2Conditioned(2, i), Image1Conditioned(2, i)*Image2Conditioned(2, i), Image2Conditioned(2, i), Image1Conditioned(1, i), Image1Conditioned(2, i), 1];
    
end

%Applying SVD
[U, D, V] = svd(A);

H = reshape(V(:,end), 3, 3)';

fMatrix = Image2TransformationMatrix' * H * Image1TransformationMatrix; %Computing the reverse conditioning

fundamentalMatrix = fMatrix(:,:)/fMatrix(end,end);  %Normalizing
end


function [PN, P] = CreateProjectionMatrix(FundamentalMatrix)
%Projection Matrix Perpendicular to the Object for Image 1
PN = [eye(3) zeros(3,1)]; 

%Calulating the epipoles of Image 2
[U, D, V] = svd(FundamentalMatrix);

Image2Ep = U(:,end);
disp("The Epipole of Image 2 is:");
disp(Image2Ep);

%Obtaining the skew-symmetric Matrix
aX = [0, -Image2Ep(3,1), Image2Ep(2,1);
    Image2Ep(3,1), 0, -Image2Ep(1,1);
    -Image2Ep(2,1), Image2Ep(1,1), 0];

%Obtaining the projection Matrix for Image 2
firstPart = (aX * FundamentalMatrix);
P = firstPart + [Image2Ep, Image2Ep, Image2Ep];
P = [P, Image2Ep];

disp("Image1 Projective Matrix is:");
disp(PN);
disp("Image2 Projection Matrix is:");
disp(P);
end

function ObjectPoints = LinearTriangulation(ImageCoordinate1, ImageCoordinate2, PN, P)

% Creating the Object Point Design Matrix
A = [];

% Initialize ObjectPoints to an empty matrix
ObjectPoints = [];

for i = 1:size(ImageCoordinate1, 2)
    A = [ImageCoordinate1(1,i)*PN(3,:)-PN(1,:)
         ImageCoordinate1(2,i)*PN(3,:)-PN(2,:)
         ImageCoordinate2(1,i)*P(3,:)-P(1,:)
         ImageCoordinate2(2,i)*P(3,:)-P(2,:)];

    % Performing Singular Value Decomposition
    [Uo, Do, Vo] = svd(A);

    % Extracting the 3D coordinates from the last column of Vo
    % Accumulate the results in ObjectPoints
    ObjectPoints = [ObjectPoints, Vo(:, size(Vo,2)) ./ Vo(4, size(Vo,2))];
    end
end


function KMatrix = CallibrationMatrix(projectionMatrix)

%converting the Projection matrix PMatrix into Matrix M

M = [projectionMatrix(1,1:3); projectionMatrix(2,1:3); projectionMatrix(3,1:3)];


%Normalization of the matrix M we use the sign of the Determinant and
%Magnitude of the last row

Determinant = det(M);

LastRow = M(end,:);

Magnitude = norm(LastRow) * sign(Determinant);

normalizedM = M/Magnitude;

%Applying the QR Theorem to obtain the callibration and rotation matrix

[Q, R] = qr(inv(normalizedM));

RMatrix = inv(Q);
KMatrix = inv(R);

for i = 1:min(size(KMatrix))

    if KMatrix(i, i) < 0;

        KMatrix(i, i) = -KMatrix(i,i);

        RMatrix(:, i) = -RMatrix(:, i);
    end
end

disp("The Calibration Matrix is:");
disp(KMatrix);
end

function singularEssentialMatrix = SingularityConstraint(essentialMatrix)
%Enforcing Singularity Constraint and Trace Constraint

[U, D, V] = svd(essentialMatrix);

%Enforcing trace constraint
traceConstraint = (D(1,1) + D(2,2))/2;

%Ensuring the diagonals satify trace constraint and singularity constraint
DModified = [traceConstraint, 0, 0; 0, traceConstraint, 0; 0, 0, 0];


%Reconstructing the essential matrix
singularEssentialMatrix = U * DModified * V';

disp("The Essential Matrix is");
disp(singularEssentialMatrix);

end 

function [rotationMatrix, translationMatrix] = ResolveEssentialMatrixAmbiguity(EssentialMatrix, Image1, Image2, K1, K2, ObjectPoints)

[U, D, V] = svd(EssentialMatrix);

%Constructing four possible choices of rotation and translation
Z = [0 -1 0; 1 0 0; 0 0 1];

R1 = U * Z * V';
R2 = U * Z' * V';
t1 = U(:, 3);
t2 = -U(:, 3);

% Obtaining the Projection matrix 
P1 = K2 * [eye(3) zeros(3, 1)];
P2_1 = K1 * [R1 t1];
P2_2 = K1 * [R1 t2];
P2_3 = K1 * [R2 t1];
P2_4 = K1 * [R2 t2];

%Performing Linear Triangulation
a = LinearTriangulation(Image1, Image2, P1, P2_1);
b = LinearTriangulation(Image1, Image2, P1, P2_2);
c = LinearTriangulation(Image1, Image2, P1, P2_3);
d = LinearTriangulation(Image1, Image2, P1, P2_4);


% Choose the best solution based on the majority of reconstructed points having positive depth
depth1 = sum(a(3) > 0 & ObjectPoints(3, :) > 0);
depth2 = sum(b(3) > 0 & ObjectPoints(3, :) > 0);
depth3 = sum(c(3) > 0 & ObjectPoints(3, :) > 0);
depth4 = sum(d(3) > 0 & ObjectPoints(3, :) > 0);

if depth1 > depth2 && depth1 > depth3 && depth1 > depth4
    R = R1;
    t = t1;
elseif depth2 > depth1 && depth2 > depth3 && depth2 > depth4
    R = R1;
    t = t2;
elseif depth3 > depth1 && depth3 > depth2 && depth3 > depth4
    R = R2;
    t = t1;
else
    R = R2;
    t = t2;
end

% Output the rotation and translation matrices
rotationMatrix = R;
translationMatrix = t;

disp("The Rotational Matrix is:");
disp(rotationMatrix);

disp("The translational matrix is:");
disp(translationMatrix);

end

%Calculating the Epipolar Lines and Drawing 
function [epipolarline1, epipolarline2] = ComputeEpipolarLines(EssentialMatrix, ImageCoordinate1, ImageCoordinate2)

epipolarline1 = [];
epipolarline2 = [];

% Calculating Epipolar Line 1
for i = 1:size(ImageCoordinate1, 2)
    l = EssentialMatrix * ImageCoordinate1(:, i);
    epipolarline1(i, :) = l';
end

% Calculating Epipolar Line 2
for i = 1:size(ImageCoordinate2, 2)
    lp = EssentialMatrix' * ImageCoordinate2(:, i);
    epipolarline2(i, :) = lp';
end



disp("Image 1 Epipolar Lines are:");
disp(epipolarline1);

disp("Image 2 Epipolar Lines are:");
disp(epipolarline2);

% Calling function to display epipolar line 1
for i = 1:size(epipolarline1, 1)
    hline(epipolarline1(i, :)');
end

% Calling function to display epipolar line 2
for i = 1:size(epipolarline2, 1)
    hline(epipolarline2(i, :)');
end

end




%Calculating the Geometric Error
function geometricError = ComputeGeometricError(Image1Coordinate, Image2Coordinate, epipolarline2)
geometricError = 0; % Initialize the error

for i = 1:size(Image2Coordinate, 1)
    a = Image2Coordinate(i, 1);
    b = Image2Coordinate(i, 2);

    x = epipolarline2(i, 1)';
    y = epipolarline2(i, 2)';
    z = epipolarline2(i, 3)';

    u = Image1Coordinate(i, 1);
    v = Image1Coordinate(i, 2);

    numerator = (a*x + b*y + z)^2;
    denominator = a^2 + b^2 + u^2 + v^2;
    geometricError = geometricError + numerator/denominator;

    disp("The Geometric Error is calculated as:")
    disp(geometricError);
end
end



function  h = hline(l, varargin)
%        ==================
    if abs(l(1)) < abs(l(2))                                  % More horizontal
        xlim = get(get(gcf, 'CurrentAxes'), 'Xlim');
        if numel(xlim) < 2
            xlim = [-1, 1];  % Default xlim if not set
        end
        x1 = cross(l, [1; 0; -xlim(1)]);
        x2 = cross(l, [1; 0; -xlim(2)]);
    else                                                        % More vertical
        ylim = get(get(gcf, 'CurrentAxes'), 'Ylim');
        if numel(ylim) < 2
            ylim = [-1, 1];  % Default ylim if not set
        end
        x1 = cross(l, [0; 1; -ylim(1)]);
        x2 = cross(l, [0; 1; -ylim(2)]);
    end
    x1 = x1 / x1(3);
    x2 = x2 / x2(3);
    h = line([x1(1) x2(1)], [x1(2) x2(2)], varargin{:});
end


function PlotEpipolarLinesAndPoints(X1, X2, epipolarline1, epipolarline2)
    % Plot the first set of image coordinates and its corresponding epipolar lines
    subplot(1, 2, 1);
    scatter(X1(1, :), X1(2, :), 'r*');  % Plot image points
    hold on;
    for i = 1:size(epipolarline1, 1)
        hline(epipolarline1(i, :)');
    end
    hold off;
    title('Image 1 Epipolar Lines');

    % Plot the second set of image coordinates and its corresponding epipolar lines
    subplot(1, 2, 2);
    scatter(X2(1, :), X2(2, :), 'g*');  % Plot image points
    hold on;
    for i = 1:size(epipolarline2, 1)
        hline(epipolarline2(i, :)');
    end
    hold off;
    title('Image 2 Epipolar Lines');

    % Adjust figure properties
    sgtitle('Epipolar Lines and Corresponding Image Points');
end






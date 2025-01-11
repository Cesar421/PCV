function ProjectiveReconstruction

%________________PART 1____________________________________

%Get the Image Cordinate Points
[X1, X2] = GetImageCordinates();

%Obtain the Fundamental Matrix
fMatrix = GetFundamentalMatrix(X1,X2);

%Perform Singularity Constraint on the Fundamental Matrix
FundamentalMatrix = SingularityConstraint(fMatrix);

%Obtain Projection Matrices for the two Image Points
[P1, P2] = CreateProjectionMatrix(FundamentalMatrix);

%Obtain the Object Coordinates by Linear Triangulation
ObjectPoints1 = LinearTriangulation(X1, X2, P1, P2);

%Plot Object Points
PlotObjectPoint(ObjectPoints1);




%________________PART 2____________________________



%Reading Coordinate points from Data provided
[C1, C2, C3] = ReadControlPoint();

%Obtaining Object Points by Linear Triangulation
ObjectPoint2 = LinearTriangulation(C1, C2, P1, P2);
disp(ObjectPoint2);

%Performing 3D Homography on the Object Points
H = Get3DHomography(ObjectPoint2, C3);

% Applying the homography to Object Points 1
EuclieanObjectPoint = GetEuclieanObjectPoint(H, ObjectPoints1);

%Visualize ObjectPoint 1
PlotObjectPoint(EuclieanObjectPoint);

end



%______________PART 1 METHODS__________________________________________

function [ImageCoordinate1, ImageCoordinate2] = GetImageCordinates()
    % Open the file "bh.dat" for reading
    fh = fopen("bh.dat", "r");
    
    % Read the data from the file into matrix A (4 rows, N columns)
    A = fscanf(fh, '%f%f%f%f', [4 inf]);
    
    % Close the file after reading
    fclose(fh);
    
    % Extract the first two rows of A as ImageCoordinate1 (coordinates from the first image)
    ImageCoordinate1 = A(1:2, :);
    
    % Extract the last two rows of A as ImageCoordinate2 (coordinates from the second image)
    ImageCoordinate2 = A(3:4, :);
    
    % Set the third row of both ImageCoordinate1 and ImageCoordinate2 to 1 for homogeneous coordinates
    ImageCoordinate1(3, :) = 1;
    ImageCoordinate2(3, :) = 1;
end

function fMatrix = GetFundamentalMatrix(ImageCoordinate1, ImageCoordinate2)
    % Function to compute the fundamental matrix (to be implemented)


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

% Perform Singular Value Decomposition (SVD) on matrix A
[U, D, V] = svd(A);

% Reshape the last column of V into a 3x3 matrix and transpose it
H = reshape(V(:,end), 3, 3)';

% Compute the fundamental matrix by applying reverse conditioning
fMatrix = Image2TransformationMatrix' * H * Image1TransformationMatrix; 

% Normalize the fundamental matrix by dividing by its last element
fMatrix = fMatrix(:,:)/fMatrix(end,end);  
end

function FundamentalMatrix = SingularityConstraint(fMatrix)
%Enforcing Singularity Constraint
if det(fMatrix) == 0

    FundamentalMatrix = fmatrix;
    disp("The Fundamental matrix is:");
    disp(FundamentalMatrix);

else

    % Perform Singular Value Decomposition (SVD) on the fundamental matrix
    [U, D, V] = svd(fMatrix);
    
    % Set the last diagonal element of D to zero to enforce the singularity constraint
    D(end, end) = 0;
    % Reconstruct the fundamental matrix by multiplying U, D, and V'
    FundamentalMatrix = U * D * V';

    % Display the computed fundamental matrix
    disp("The Fundamental matrix is:");
    disp(FundamentalMatrix);
end 
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
P = (aX * FundamentalMatrix) + [Image2Ep, Image2Ep, Image2Ep];
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

 % Create a new figure window for plotting
function PlotObjectPoint(ObjectPoints)
    % Scatter plot for 3D object points with filled markers
figure; scatter3(ObjectPoints(1,:), ObjectPoints(2,:), ObjectPoints(3,:), 10, "filled");
axis square; view(32, 75); % Set axis to square for equal scaling in all dimensions

end



%________________PART 2 METHODS__________


function [X1, X2, X3] = ReadControlPoint()
fh = fopen("pp.dat","r");
A = fscanf(fh, '%f%f%f%f', [7 inf]);
fclose(fh);
X1 = A(1:2, :);
X2 = A(3:4, :);
X3 = A(5:7, :);
X1(3,:) = 1;
X2(3,:) = 1;
X3(4, :) = 1;

end

function Homography = Get3DHomography(ImageCoordinate1, ImageCoordinate2)

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
Image1TranslatedMatrix = [1, 0, 0, -Image1Midpoint(1,1); 0, 1, 0, -Image1Midpoint(2,1); 0, 0, 1, -Image1Midpoint(3,1); 0, 0, 0, 1];
Image2TranslatedMatrix = [1, 0, 0, -Image2Midpoint(1,1); 0, 1, 0, -Image2Midpoint(2,1); 0, 0, 1, -Image2Midpoint(3,1); 0, 0, 0, 1];

%Scaled Matrix
Image1ScaledMatrix = [1/Image1Scaled(1,1), 0, 0, 0; 0, 1/Image1Scaled(2,1), 0, 0; 0, 0, 1/Image1Scaled(3,1), 0; 0, 0, 0, 1];
Image2ScaledMatrix = [1/Image2Scaled(1,1), 0, 0, 0; 0, 1/Image2Scaled(2,1), 0, 0; 0, 0, 1/Image2Scaled(3,1), 0; 0, 0, 0, 1];

%Image 1 and 2 Transformation Matrix
Image1TransformationMatrix = Image1ScaledMatrix * Image1TranslatedMatrix;
Image2TransformationMatrix = Image2ScaledMatrix * Image2TranslatedMatrix;
disp(Image1TransformationMatrix);
disp(Image2TransformationMatrix);

%Image 1 and 2 Conditioned Coordinates
Image1Conditioned = [Image1Translated(1,:)*Image1TransformationMatrix(1,1); Image1Translated(2,:)*Image1TransformationMatrix(2,2); Image1Translated(3,:)*Image1TransformationMatrix(3,3)];
Image2Conditioned = [Image2Translated(1,:)*Image2TransformationMatrix(1,1); Image2Translated(2,:)*Image2TransformationMatrix(2,2); Image2Translated(3,:)*Image2TransformationMatrix(3,3)];

% Formulate the design matrix A for the homography computation

A = [];

% Loop through each point in the conditioned image points
for i = 1:size(Image1Conditioned, 2)
    
    % Append a row to A based on the point correspondences
    % Each row represents an equation derived from the homography constraints
    
    A = [A; 
        -Image1Conditioned(1, i), -Image1Conditioned(2, i), -Image1Conditioned(3, i), -1, zeros(1, 4), zeros(1, 4), Image1Conditioned(1, i) * Image2Conditioned(1, i), Image1Conditioned(2, i) * Image2Conditioned(1, i), Image1Conditioned(3, i) * Image2Conditioned(1, i), Image2Conditioned(1, i);
        zeros(1, 4), -Image1Conditioned(1, i), -Image1Conditioned(2, i), -Image1Conditioned(3, i), -1, zeros(1, 4), Image1Conditioned(1, i) * Image2Conditioned(2, i), Image1Conditioned(2, i) * Image2Conditioned(2, i), Image1Conditioned(3, i) * Image2Conditioned(2, i), Image2Conditioned(2, i);
        zeros(1, 4), zeros(1, 4), -Image1Conditioned(1, i), -Image1Conditioned(2, i), -Image1Conditioned(3, i), -1, Image1Conditioned(1, i) * Image2Conditioned(3, i), Image1Conditioned(2, i) * Image2Conditioned(3, i), Image1Conditioned(3, i) * Image2Conditioned(3, i), Image2Conditioned(3, i)];
end

% Display the final design matrix A
disp(A);



%Applying SVD
[U, D, V] = svd(A); % Perform Singular Value Decomposition (SVD) on matrix A
disp(V);% Display the matrix V from the SVD decomposition

H = reshape(V(:,end), 4, 4)';% Reshape the last column of V into a 4x4 matrix and transpose it

% Compute the homography matrix by applying reverse conditioning
Homography = inv(Image2TransformationMatrix) * H * Image1TransformationMatrix; %Computing the reverse conditioning

Homography = Homography(:,:)/Homography(end,end);  %Normalizing
% Display the computed 3D homography matrix
disp("The 3D Homography Matrix is:");
disp(Homography);
end


function EuclideanObjectPoints = GetEuclieanObjectPoint(Homography, ObjectPoint)
    % Initialize an empty array to store the Euclidean object points
    EuclideanObjectPoints = [];

    % Apply the homography transformation to the object points
    Points = Homography * ObjectPoint;

    % Loop through each column of the ObjectPoint matrix (each point)
    for i = 1:size(ObjectPoint, 2)
        % Normalize the point by dividing by the 4th element (homogeneous coordinate)
        EuclideanObjectPoints = [EuclideanObjectPoints, Points(:, i)./Points(4,i)];
    end
end

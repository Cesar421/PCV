%  Example call of the GEOKOR function
%  -------------------------------------------
% function test
%          ====
% f = imread('image1.ppm');                          % Read example images
% g = imread('image2.ppm');
% H = [1 0.1 10; 0 1 20; 0 0 1];   % Example for projective Transformation
% i = geokor_octave(H, f, g);     % Transform image f using H and combine with g
% figure; imshow(i);                                   % Show image mosaic

function i = geokor_octave(H, f, g)
%        ==========================
% pkg load image
[fy, fx, fc] = size(f);                                   % f original size
[gy, gx, ~] = size(g);                                    % g original size
l = H * [1 1 fx fx; 1 fy 1 fy; 1 1 1 1];
l = round(l./l(3,:));                           % position of transformed f

bx1 = min([l(1, :) 1]); bx2 = max([l(1, :) gx]);  % common boundaries for g
by1 = min([l(2, :) 1]); by2 = max([l(2, :) gy]);        % and transformed f
i = zeros(by2-by1+1, bx2-bx1+1, fc, "uint8");          % empty common image

ix = (0 : max(l(1,:))-min(l(1,:))) + max(min(l(1,:)), 1);   % position of f
iy = (0 : max(l(2,:))-min(l(2,:))) + max(min(l(2,:)), 1);  % in common area
i(iy, ix, :) = imperspectivewarp(f, H, 'cubic'); % transform image f with H

rx = (1 : gx) - bx1+1;                       % position of g in common area
ry = (1 : gy) - by1+1;
i(ry, rx, :) = max(g, i(ry, rx, :));             % combine brightest pixels

addpath('C:\Users\Dave\Desktop\Computer Vision\Final Project\Constellation Detector\');
I1 = imread('nightsky_2.png');
%I1 = im2double(I1);
%I1 = preprocess_image(I1,20);

%I1 = imresize(I1,0.1);
gauss = fspecial('gaussian',[2 2],10);
gauss_conv = imfilter(I1,gauss);
laplacian = I1 - gauss_conv;
laplacian_gray = rgb2gray(laplacian);
star_corners = corner(laplacian_gray,'harris');
size(star_corners(:,1))
for i = 1:size(star_corners(:,1))
   for j=1:size(star_corners(:,1))
   similarity_matrix(i,j) = 1 + cos(pdist2(star_corners(i),star_corners(j),'squaredeuclidean'));
   end
end

se = strel('sphere',5);
laplacian_dilate = imdilate(laplacian,se);
stats = regionprops(laplacian_dilate,'pixelIdxlist');
[center radii] = imfindcircles(laplacian_dilate,[5 30]);

K = sort(radii,'ascend');
%%%%%%%%%%%% 
%%%%%%%%%%%% radius = ( (1/2.5)^(mag) )/(2*pi));
%%%%%%%%%%%% 
radii_bright = find(radii > median(radii));
radii_dark = find(radii <= median(radii));

center_bright = zeros(size([center(length(radii_bright),:)]));
center_dark = zeros(size([center(length(radii_dark),:)]));

for i=1:length(radii_dark)
    for j = 1:2
    center_dark(i,j) = center(i,j);
    end
end

for i=1:length(radii_bright)
    for j = 1:2
        center_bright(i,j) = center((i+length(radii_dark)),j);
    end
end
% Plot bright circles in blue
viscircles(center_bright, radii_bright,'Color','b');
% Plot dark circles in dashed red boundaries
viscircles(center_dark, radii_dark,'Color','r');

for i = 1:size(center_bright(:,1))
   for j=1:size(center_bright(:,1))
   similarity_matrix_bright_magnitude(i,j) = 1 + cos(pdist2(center_bright(i),center_bright(j),'squaredeuclidean'));
   end
end
bright_constellations = slowdbscan(similarity_matrix_bright_magnitude,3,5);

for i = 1:size(center_dark(:,1))
   for j=1:size(center_dark(:,1))
   similarity_matrix_dark_magnitude(i,j) = 1 + cos(pdist2(center_dark(i),center_dark(j),'squaredeuclidean'));
   end
end
dark_constellations = slowdbscan(similarity_matrix_dark_magnitude,3,5);






% DAVID MURRAY
% 113818946
% HOG - PROJECT 4
addpath('C:\Users\Dave\Desktop\Computer Vision\Final Project');
oimg = imread('nightsky_3.png');
oimg = im2double(oimg);
%laplacian_blob = oimg - imfilter(oimg,fspecial('gaussian',[2 2],10));
figure(1)
%imshow(laplacian_blob);
oimg = rgb2gray(oimg);
imgheight=512;
img=imresize(oimg,[imgheight,imgheight/size(oimg,1)*size(oimg,2)]);


[x y]=ndgrid(1:8:776,1:8:512);   %max_height 300px -- need square list = size(:,1)*dims
list=[x(:) y(:)];
tic;
%data = load('gantrycrane_hog.mat');
h = hog1(oimg,list,8);
load('h_star_data');
toc; % compute hog at the locations in the list
visualize_hog_list(h,list,img);






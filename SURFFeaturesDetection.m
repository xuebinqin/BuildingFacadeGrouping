% SURF Features
% Detect interest points and mark their locations
close
clear all

I = imread('image1.jpg');
I1 = I;
I = rgb2gray(I);
points = detectSURFFeatures(I,'MetricThreshold',1000.0,'NumOctaves',3,'NumScaleLevels',4);
imshow(I1); hold on;
X = points.Location;
% plot(points.selectStrongest(100));
plot(X(:,1),X(:,2),'r+','MarkerSize',25,'LineWidth',3);
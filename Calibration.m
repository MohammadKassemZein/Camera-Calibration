
%{
Application of camera calibration using Zhang's algorithm
Non-linear optimization is not considered
The user can enter the number of Images of 2D plane checker Boards to be used
Author: Mohammad Kassem Zein
Institution: American University of Beirut-Vision and Robotics Lab
Email:mhamadk.zein@gmail.com
%}

%fullfile: builds a full file path for each image.
%Note: Adjust Accordingly.

nimages=input('Enter the number of images: ');
for i = 1:nimages;
  imageFileName = sprintf('Image%d.tif', i);
  imageFileNames{i} = fullfile('C:\','Users','MohammadZein','Desktop','CameraCalibration-MATLAB',imageFileName);
end
%imagePoints: points detected on the images
%boardSize: Checker board dimensions in terms of number of square
[imagePoints,boardSize,imagesUsed] = detectCheckerboardPoints(imageFileNames);
imageFileNames = imageFileNames(imagesUsed);
% for i = 1:numel(imageFileNames)   %% numel: number of Array elements.
%   I = imread(imageFileNames{i});
%   subplot(5,4,i);
%   imshow(I);
%   hold on;
%   plot(imagePoints(:,1,i),imagePoints(:,2,i),'ro');
% end
%worldPoints: Return M-by-2 matrix containing M[x,y] corner coordinates for
%the squares on the checker board
squareSize_mm=30;
worldPoints=generateCheckerboardPoints(boardSize,squareSize_mm);
[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints);
% disp(cameraParams);
%%% info=imfinfo(imageFileNames{1});
%%% disp(info);

%Normalization of All image points:
readimage=imread(imageFileNames{1});
[h,w]=size(readimage);
N=[2/w 0 -1; 0 2/h -1; 0 0 1]; %%normalization matrix
h_imagePoints=zeros(size(imagePoints,1),3,numel(imageFileNames));%%image points in Homogeneous coordinates
h_worldPoints=zeros(size(worldPoints,1),3);%%world points in Homogeneous coordinates
for i = 1:numel(imageFileNames)
    for j=1:size(imagePoints,1)
        h_imagePoints(j,1,i)=imagePoints(j,1,i);
        h_imagePoints(j,2,i)=imagePoints(j,2,i);
        h_imagePoints(j,3,i)=1;
        h_worldPoints(j,1)=worldPoints(j,1);
        h_worldPoints(j,2)=worldPoints(j,2);
        h_worldPoints(j,3)=1;
    end
end
% disp(worldPoints);
% disp(imagePoints);
% disp(h_worldPoints);
%%%Calculation of the homography matrix using Direct Linear Transformation
%%%(DLT):
npoints=size(imagePoints,1);   %%number of points (i.e. 156 points)
x_worldPoints=zeros(1,size(imagePoints,1));
y_worldPoints=zeros(1,size(imagePoints,1));
x_imagePoints=zeros(1,size(imagePoints,1),numel(imageFileNames)); %Array of matrix
y_imagePoints=zeros(1,size(imagePoints,1),numel(imageFileNames));

for i=1:npoints
    x_worldPoints(i)=h_worldPoints(i,1);  %% first column corresponding to the x coordinates of the points
    y_worldPoints(i)=h_worldPoints(i,2);  %% second column corresponding to the y coordinates of the points
end
% disp(x_worldPoints);
for i=1:numel(imageFileNames)
    for j=1:npoints
        x_imagePoints(1,j,i)=h_imagePoints(j,1,i);
        y_imagePoints(1,j,i)=h_imagePoints(j,2,i);
    end
end
%Setting up the Matrix system for the points such that Lh=0;
L=zeros(2*npoints,9,numel(imageFileNames));
for i=1:numel(imageFileNames)
    for j=1:npoints
        L(2*j-1,:,i)=[x_worldPoints(j),y_worldPoints(j),1,0,0,0,-x_worldPoints(j)*x_imagePoints(1,j,i),-y_worldPoints(j)*x_imagePoints(1,j,i),-x_imagePoints(1,j,i)];
        L(2*j,:,i)=[0,0,0,x_worldPoints(j),y_worldPoints(j),1,-x_worldPoints(j)*y_imagePoints(1,j,i),-y_worldPoints(j)*y_imagePoints(1,j,i),-y_imagePoints(1,j,i)];
    end
end
% disp(L(1,:,1));
%Solve for the homography using svd which separates L to the product of
%three matrices U,S,V
h=zeros(1,9,numel(imageFileNames));
if npoints==4
    h=null(L);
else
    for i=1:numel(imageFileNames)
        [U,S,V]=svd(L(:,:,i));
        h(:,:,i)=V(:,9);
    end
end
% disp(h);
%Reshape Homography Matrices
H=zeros(3,3,numel(imageFileNames));
for i=1:numel(imageFileNames)
    H(:,:,i)=[h(1,1,i),h(1,2,i),h(1,3,i);h(1,4,i),h(1,5,i),h(1,6,i);h(1,7,i),h(1,8,i),h(1,9,i)];
end
% disp(H);
Hi=zeros(3,3,numel(imageFileNames));

for i=1:numel(imageFileNames)
    Hi(:,:,i)=(H(:,:,i))';
end
% disp(Hi); 
% disp(Hi(1,1,1)*Hi(1,2,1));
%Setting up V to solve Vb=0 after stacking n equations

v=zeros(2,6,numel(imageFileNames));
for i=1:numel(imageFileNames)
    v11=[Hi(1,1,i)*Hi(1,1,i),Hi(1,1,i)*Hi(1,2,i)+Hi(1,2,i)*Hi(1,1,i),Hi(1,2,i)*Hi(1,2,i),Hi(1,3,i)*Hi(1,1,i)+Hi(1,1,i)*Hi(1,3,i),Hi(1,3,i)*Hi(1,2,i)+Hi(1,2,i)*Hi(1,3,i),Hi(1,3,i)*Hi(1,3,i)];
    v22=[Hi(2,1,i)*Hi(2,1,i),Hi(2,1,i)*Hi(2,2,i)+Hi(2,2,i)*Hi(2,1,i),Hi(2,2,i)*Hi(2,2,i),Hi(2,3,i)*Hi(2,1,i)+Hi(2,1,i)*Hi(2,3,i),Hi(2,3,i)*Hi(2,2,i)+Hi(2,2,i)*Hi(2,3,i),Hi(2,3,i)*Hi(2,3,i)];
    v12=[Hi(1,1,i)*Hi(2,1,i),Hi(1,1,i)*Hi(2,2,i)+Hi(1,2,i)*Hi(2,1,i),Hi(1,2,i)*Hi(2,2,i),Hi(1,3,i)*Hi(2,1,i)+Hi(1,1,i)*Hi(2,3,i),Hi(1,3,i)*Hi(2,2,i)+Hi(1,2,i)*Hi(2,3,i),Hi(1,3,i)*Hi(2,3,i)];
    v(1,:,i)=v12;
    v(2,:,i)=v11-v22;
end

% disp(v);
%Stacking of N Equations:
Vi=zeros(2*numel(imageFileNames),6);
for i=1:2*numel(imageFileNames)
    c=mod(i,2); % check if even
    if(c==0)
        Vi(i,:)=v(2,:,ceil(i/2));
    else
        Vi(i,:)=v(1,:,ceil(i/2));
    end
end
% disp(Vi);
[U,S,V]=svd(Vi);
b=V(:,6);
% disp(b);

%Initialize Intrinsic Parameters:

v0=((b(2)*b(4))-(b(1)*b(5))/(b(1)*b(3))-(b(2))^2);
s=b(6)-((b(4)^2)+v0*(b(2)*b(4)-b(1)*b(5)))/b(1);
alpha=sqrt(s/b(1));
beta=sqrt((s*b(1))/(b(1)*b(3)-b(2)^2));
skew=(-b(2)*(alpha^2)*beta)/s;
u0=((skew*v0)/beta)-((b(4)*(alpha^2))/s);

% disp(v0);
% disp(u0);
% disp(alpha);
% disp(beta);
% disp(skew);

%Formation of Intrinsic Matrix:
A=[alpha,skew,u0;0,beta,v0;0,0,1];
disp(A);
%Intrinsic Parameters using MATLAB Function:
disp(cameraParams);
%Calculation of Extrinsic Parametres (Rotation and Translation Matrices)


















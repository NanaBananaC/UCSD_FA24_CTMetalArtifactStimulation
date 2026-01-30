clear
clc

load("simulate_sample_v2.mat")
subimage_switch = 0;
y_edge = [30, 230]; x_edge = [170, 370];  % For plotting subimage
H = simulate_sample_2.inv_40_metal;
figure
subplot(1,2,1)
imagesc(H)
title("Original")

subplot(1,2,2)
H = rescale(H);
imagesc(H)
title('H rescaled 0 to 1')
colormap gray;
%imread('Original-CT-image.png');
%H = imread('HipMetalArtifact.jpeg');

%Taking the histogram
[counts,x] = imhist(H);
figure(2)
stem(x,counts);
grid on;


%%
%{
H = imread('Original-CT-image.png');
%H = rgb2gray(H);
%H = imread('HipMetalArtifact.jpeg');
figure(1)
imagesc(H);
colormap gray;

%Taking the histogram
[counts,x] = imhist(H);
figure(2)
stem(x,counts);
grid on;

%%NEW EDITS
ed = edge(H, 'canny');
figure
imagesc(ed)
colormap gray;

SE = strel("rectangle",[2 6]);
closeim = imclose(ed, SE);
figure
imagesc(closeim);
colormap gray;

boxfilter = (1/81).*ones(9,9);

%
%ha = inpaintExemplar(H,logical(closeim));
H = roifilt2(boxfilter,H,logical(closeim));
figure
imagesc(H)
colormap gray;
%}
%}


%%

%%Otsu's method chooses a threshold that minimizes the intraclass variance of the thresholded black and white pixels.
%T = otsuthresh(counts); 
%{
T = adaptthresh(H);
%figure(3)
%imshow(T)
BW = imbinarize(H,T);
figure
%figure
subplot(1,3,1)
imshow(BW)
title("Adaptive Threshold")
axis('tight')
%}


%Percentile thresholding
%highPercentile = 0.99; % Adjust this based on your image
thresholdValue = 0.5;
%180; original dental image
%prctile(H(:), highPercentile * 100);
BW_metal = H >= thresholdValue;
BW_metal = bwareaopen(BW_metal, 80); % Removes small bone-like structures
BW_metal = imclose(BW_metal, strel('disk', 3)); % Close gaps in the metal mask

figure
subplot(1,2,1)
imagesc(BW_metal), title('Histogram Threshold');
axis('tight')
%Connected Component analysis to deal with extra bone
%Temporarily don't need connected components

%******************************
%{
CC = bwconncomp(BW_metal);
% Filter based on size or shape to retain only large regions
stats = regionprops(CC, 'Area', 'PixelIdxList');
minSize = 650; % Adjust based on expected metal artifact size
for i = 1 : numel(stats)
    if stats(i).Area < minSize
        BW_metal(stats(i).PixelIdxList) = 0;
    end
end
subplot(1,2,2) 
imshow(BW_metal), title('Connected Components');
axis('tight')
%}

metal = zeros(size(H));
metal(BW_metal==1) = H(BW_metal==1);

%{
boxfilter = (1/9).*ones(3,3);
cleanbox = filter2(boxfilter,H);
imshow(cleanbox)

%
filtered = medfilt2(H);
figure(5)
imshow(filtered);
%}

%% Creating the interpolated image

theta = 0:0.1:179;

% Generate sinograms for each image
sinogram_input = radon(H, theta); 
sinogram_metal_only = radon(metal, theta);  % Metal-only image sinogram
metal_sin_mask = zeros(size(sinogram_metal_only));
metal_sin_mask(sinogram_metal_only~=0) = 1;


%Copying over
sin_interp = sinogram_input;

for ang = 1:size(sin_interp, 2)
    % Extract the column for this angle
    projection = sin_interp(:, ang);
    mask = metal_sin_mask(:, ang); % Get the metal mask for this angle
    
    % Find metal indices
    metal_indices = find(mask); % Indices where metal is present
    
    % Interpolate values only at the metal indices
    non_metal_indices = find(~mask); % Indices where metal is not present
    projection(metal_indices) = interp1(non_metal_indices, projection(non_metal_indices), metal_indices, 'linear','extrap');
   
    % Replace the projection with the interpolated values
    sin_interp(:, ang) = projection;
end


%test = inpaintCoherent(sinogram_input, logical(metal_sin_mask));

interp_img = iradon(sin_interp, theta,'linear', 'Ram-Lak');

figure;
subplot(1,2,1)
imagesc(sinogram_input);
title('Original Sinogram');

subplot(1,2,2)
imagesc(sin_interp);
%sin_interp
colormap(gray);
title('Interpolated Sinogram (Metal Replaced)');

figure
subplot(1,2,1)
imagesc(H);
title('Original Image');

subplot(1,2,2)
imagesc(interp_img);
colormap(gray);
title('Interpolated Image (Metal Replaced)');


%% Tissue Segmentation

%{
tissueThreshold = 0;
BW_tissue = H;
BW_tissue(BW_metal==1) = NaN;

figure(5)
subplot(2,2,1)
imagesc(BW_metal)
colormap('gray')
axis('square')
title("Metal Region")
subplot(2,2,2)
imshow(BW_tissue)
title("Image w/o metal")
axis('square')
%}

% Define a lower threshold for tissue (experiment with the value)
scaled_interp_img = rescale(interp_img);
figure
subplot(1,3,1)
imagesc(interp_img)
colormap gray;
title("Interpolated")
%[counts,x] = imhist(scaled_interp_img);
%figure(2)
%stem(x,counts);
%grid on;

sigma = 3;
blurred = imgaussfilt(interp_img, sigma);
subplot(1,3,2)
imagesc(blurred)
colormap gray;
title("Interpolated + Blur")

tissueThreshold = 0.07;
%= adaptthresh(interp_img)*0.8;  % Adjust factor as necessary

BW_tissue = blurred > tissueThreshold;
%Exclude the outer bars

top_border = 20; % Start row for the top border
bottom_border = size(BW_tissue, 1); % End row for the bottom border
left_border = 1; % Start column for the left border
right_border = size(BW_tissue, 2) - 20; % End column for the right border

% Set everything outside this defined ROI to zero
BW_tissue(1:top_border, :) = 0; % Clear rows above the head
BW_tissue(:, right_border:end) = 0; % Clear columns on the right side
BW_tissue(:, 1:left_border) = 0; % Clear columns on the left side
BW_tissue(bottom_border:end, :) = 0; % Clear rows below the head

BW_tissue = BW_tissue(2:end-1,2:end-1);
%Remove metal teeth from mask
%BW_tissue(BW_metal==1) = 0;

%BW_tissue = imbinarize(interp_img, tissueThreshold);
% Remove metal regions from tissue mask
%BW_tissue = BW_tissue(2:end-1, 2:end-1) & ~BW_metal;
%BW_tissue = imbinarize(BW_tissue,0.5);
subplot(1,3,3)
imagesc(BW_tissue)
title('Tissue Mask');
axis('square')
colormap gray;
%title("Interpolated")



%Creating tissue only image
%{
tissues_mean = mean(interp_img(BW_tissue==1));
%tissues = zeros(size(interp_img));
%tissues(BW_tissue==1) = tissues_mean;
tissue_classified_img = interp_img;
tissue_classified_img(BW_tissue == 1) = tissues_mean; % Set tissue regions to average intensity
subplot(1,2,2)
imagesc(tissue_classified_img)
title('Tissue Classified Image');

imagesc(interp_img)
%}


%% Tissue only Image
interp_img = interp_img(2:end-1,2:end-1);
tissues_mean = mean(mean(interp_img.*BW_tissue));
%tissues = zeros(size(interp_img));
%tissues(BW_tissue==1) = tissues_mean;
tissue_classified_img = interp_img.*BW_tissue;
tissue_classified_img(BW_tissue == 1) = tissues_mean; % Set tissue regions to average intensity
%tissue_classified_img(~BW_tissue) = 0;
figure
imagesc(tissue_classified_img)
title('Tissue Classified Image');
colormap gray;



%% Forward Projections
% Define angles for the Radon transform (common choice is 0 to 180 degrees)
theta = 0:0.1:179;

% Generate sinograms for each image
sinogram_input = radon(H, theta);  % Original image sinogram
sinogram_metal_only = radon(metal, theta);  % Metal-only image sinogram
sinogram_tissue_classified = radon(tissue_classified_img, theta);  % Tissue-classified image sinogram

figure;
subplot(1, 3, 1);
imagesc(sinogram_input);
colormap gray; 
title('Input Image Sinogram');
colorbar;
%clim([0,10e4])
subplot(1, 3, 2);
imagesc(sinogram_metal_only);
colormap gray; 
colorbar
title('Metal-Only Sinogram');
%clim([0,10e4])
subplot(1, 3, 3);
imagesc(sinogram_tissue_classified);
colormap gray; 
title('Tissue-Classified Sinogram');
%clim([0,10e4])
c = colorbar;

%Error Sinogram
error_sin = zeros(size(sinogram_input));
error_sin = sinogram_input-sinogram_tissue_classified;
figure
imagesc(error_sin);
title("Error sinogram")
colormap gray;

%Making metal sinogram mask
metal_sin_mask = zeros(size(sinogram_metal_only));
metal_sin_mask(sinogram_metal_only~=0) = 1;
%metal_sin_mask = sinogram_metal_only > 0;
%sum(sum(sinogram_metal_only == 1|0))

%error_mask = error_sin.*metal_sin_mask;
error_mask = error_sin;
error_mask(metal_sin_mask == 0) = 0;  % Zero out non-metal regions

figure
subplot(1,4,1)
imagesc(sinogram_metal_only)
title("Metal Sinogram")

subplot(1,4,2)
imagesc(metal_sin_mask);
title("Mask of metal Sinogram")

subplot(1,4,3)
imagesc(error_sin);
title("Error Sinogram")

atten_thresh = 0.7;
error_mask(error_mask>atten_thresh) = atten_thresh;

subplot(1,4,4)
imagesc(error_mask);
title("Mask Applied")
colormap gray;


%{
error_mask(error_mask>atten_thresh) = atten_thresh;
figure
imagesc(error_mask);
title("Mask Applied w/ Attenuation")
colormap gray;
%}

%{
%Compare the generated metal sinogram and masked
figure
subplot(1,2,1)
imagesc(error_mask);
title("Mask Generated")

subplot(1,2,2)
imagesc(sinogram_metal_only);
title("Direct from metal-only image")
colormap gray;
%}

%% Correction image

artifact_img = iradon(error_mask, theta,'linear', 'Ram-Lak');
size(artifact_img)
figure
subplot(1,3,1)
imagesc(artifact_img)
title("Correction image")
subtitle('Image of artifacts only')
colormap gray;

%artifact_img = artifact_img*0.7;

H = double(H);
%artifact_img = artifact_img(2:end-1,2:end-1);
%rescale(artifact_img, 0, 255);

final = H-artifact_img(2:end-1,2:end-1);

subplot(1,3,2)
imagesc(final)
colormap gray;
title("MAR Image")

subplot(1,3,3)
imagesc(H)
colormap gray;
title("Original Image")

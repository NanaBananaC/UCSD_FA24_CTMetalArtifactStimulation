%% [simulate_v2.m]
% Change log: 
%   - Change angular resolution (delta theta) from 1 deg to 0.1 deg.
%   - Determine attenuation threshold based on the max value of the sinogram w/o metal: atten_thres = max(max(sino_40)) * 1.05;

close all; clear all; clc;

%% Load sample image (in HU unit) & metal mask
load('./sample_data/sample_2.mat');
image_HU = sample.image;
metal_mask = sample.metal;

%% Convert mu for energy level 1
xray_characteristic_data = readtable('./src/xray_characteristic_data.csv');
mu_water = xray_characteristic_data{40, 'Water'};
mu_air = 0;

image_40 = zeros(size(image_HU));
for i = 1:size(image_HU, 1)
    for j = 1:size(image_HU, 2)
        image_40(i, j) = image_HU(i, j) * mu_water / 1000 + mu_water;
    end
end

% Generate image_mu_metal
image_40_metal = image_40;
for i = 1:size(image_HU, 1)
    for j = 1:size(image_HU, 2)
        if metal_mask(i, j) == 1
            image_40_metal(i, j) = xray_characteristic_data{40, 'Iron'};
        end
    end
end

%% Use Radon transform to get monoenergy sinogram
theta = 0:0.1:180;
sino_40 = radon(image_40, theta);
sino_40_metal = radon(image_40_metal, theta);

% Attenuation
atten_thres = max(max(sino_40)) * 1.05;  % Determine attenuation threshold based on the max value of sinogram w/o metal
sino_40_metal_a = sino_40_metal;
for i = 1:size(sino_40_metal, 1)
    for j = 1:size(sino_40_metal, 2)
        if sino_40_metal(i, j) > atten_thres
            sino_40_metal_a(i, j) = atten_thres;
        end
    end
end

% inverse Radon
iimage_40 = iradon(sino_40, theta, "linear", "Hann");
iimage_40_metal = iradon(sino_40_metal_a, theta, "linear", "Hann");

%% Visualization
datacursormode on; colormap(gray(256));
subplot(2, 3, 1); imagesc(image_40, [0, 1.5]); axis('square'); title("μ Image (Energy Level 40)");
subplot(2, 3, 4); imagesc(image_40_metal, [0, 1.5]); axis('square'); title("μ Image w/ Metal (Energy Level 40)");
subplot(2, 3, 2); imagesc(theta, 1:size(sino_40, 1), sino_40); axis('square'); title("Sinogram (Energy Level 40)");
xlabel('projection angle θ [deg]'); ylabel('linear displacement - 𝓁');
subplot(2, 3, 5); imagesc(theta, 1:size(sino_40_metal_a), sino_40_metal_a); axis('square'); title("Attenuated Sinogram w/ Metal (Energy Level 40)");
xlabel('projection angle θ [deg]'); ylabel('linear displacement - 𝓁');
subplot(2, 3, 3); imagesc(iimage_40, [0, 1.5]); axis('square'); title("Back-projected Image");
subplot(2, 3, 6); imagesc(iimage_40_metal, [0, 1.5]); axis('square'); title("Back-projected Image w/ Metal Artifacts");

%% Beam-hardening Simulation
% Ref: Correction for beam hardening in computed tomography
energy_composition = 1:120;
sino_poly = zeros(size(sino_40));

for i = energy_composition
    disp(i)
    mu_water_E = xray_characteristic_data{i, 'Water'};
    intensity = xray_characteristic_data{i, 'Intensity'};
    d = radon(image_40 / mu_water, theta);
    if intensity ~= 0
        T = intensity * exp(-mu_water_E * d);
        sino_poly = sino_poly + T;
    end
end
sino_poly = -log(sino_poly);
iimage_poly = iradon(sino_poly, theta, "linear", "Hann");

%% Visualization
figure;
datacursormode on;
subplot(2, 2, 1); imagesc(sino_40); axis('square');
xlabel('projection angle θ [deg]'); ylabel('linear displacement - 𝓁');
subplot(2, 2, 2); imagesc(iimage_40); axis('square');
subplot(2, 2, 3); imagesc(sino_poly); axis('square');
xlabel('projection angle θ [deg]'); ylabel('linear displacement - 𝓁');
subplot(2, 2, 4); imagesc(iimage_poly); axis('square');
colormap(gray(256))

%% Beam Hardening Simulation with Metal Data
energy_composition = 1:120;
sino_poly_metal = zeros(size(sino_40_metal));

for i = energy_composition
    disp(i)
    mu_water_E = xray_characteristic_data{i, 'Water'};
    intensity = xray_characteristic_data{i, 'Intensity'};
    d = radon(image_40_metal / mu_water, theta);
    if intensity ~= 0
        T = intensity * exp(-mu_water_E * d);
        sino_poly_metal = sino_poly_metal + T;
    end
end

sino_poly_metal = -log(sino_poly_metal);
% Attenuation
atten_thres = max(max(sino_poly)) * 1.05;
sino_poly_metal_a = sino_poly_metal;
for i = 1:size(sino_poly_metal, 1)
    for j = 1:size(sino_poly_metal, 2)
        if sino_poly_metal(i, j) > atten_thres
            sino_poly_metal_a(i, j) = atten_thres;
        end
    end
end
iimage_40_metal = iradon(sino_40_metal_a, theta, 'linear', 'Hann');
iimage_poly_metal = iradon(sino_poly_metal_a, theta, "linear", "Hann");

%% Visualization
figure;
datacursormode on;
subplot(2, 3, 1); imagesc(theta, 1:size(sino_40_metal_a), sino_40_metal_a); axis('square');
xlabel('projection angle θ [deg]'); ylabel('linear displacement - 𝓁'); title("Sinogram w/ Metal (Monochromatic Energy)");
subplot(2, 3, 2); imagesc(iimage_40_metal, [-0.1 0.5]); axis('square'); title("Back-projected Image w/ Metal (Polychromatic Energy)");

subplot(2, 3, 4); imagesc(theta, 1:size(sino_poly_metal_a), sino_poly_metal_a); axis('square');
xlabel('projection angle θ [deg]'); ylabel('linear displacement - 𝓁'); title("Sinogram w/ Metal (Polychromatic Energy)");
subplot(2, 3, 5); imagesc(iimage_poly_metal, [-0.1 0.4]); axis('square'); title("Back-projected Image w/ Metal (Polychromatic Energy)");
colormap(gray(256))

%% Save data
simulate_sample_2.raw_image_HU = image_HU;
simulate_sample_2.mu_40_clean = image_40;
simulate_sample_2.sino_40_clean = sino_40;
simulate_sample_2.sino_40_metal = sino_40_metal_a;
simulate_sample_2.sino_poly_clean = sino_poly;
simulate_sample_2.sino_poly_metal = sino_poly_metal_a;
simulate_sample_2.inv_40_clean = iimage_40;
simulate_sample_2.inv_40_metal = iimage_40_metal;
simulate_sample_2.inv_poly_clean = iimage_poly;
simulate_sample_2.inv_poly_metal = iimage_poly_metal;
save('simulate_sample_v2.mat', "simulate_sample_2");

% save("simulate_result.mat")
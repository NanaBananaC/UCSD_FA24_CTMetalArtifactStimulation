%% Load data
clear all; close all; clc; 
load("../metal_artifact_simulation-master/simulate_sample_v2.mat")

subimage_switch = 1;
y_edge = [30, 230]; x_edge = [170, 370];  % For plotting subimage

% Project parameters
theta = 1:0.5:180;

%% Load data
inv_40_metal = simulate_sample_2.inv_40_metal;
colormap(gray(256))
subplot(3, 3, 1)
imagesc(inv_40_metal, [0, 1]); axis square; title("(a) Simulated image - inv 40 metal");
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end

%% Generate metal masks
metal_mask = inv_40_metal > 2;  % manually tuned
subplot(3, 3, 4);
imagesc(metal_mask); title("(d) metal mask 1 of (a)"); axis('square');
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end

metal_mask_2 = inv_40_metal > 0.75;  % manually tuned
subplot(3, 3, 7); 
imagesc(metal_mask_2); title("(g) metal mask 2 of (a)"); axis('square');
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end

%% Forward projection of the metal regions
sino_metal = radon(inv_40_metal .* metal_mask, theta);
subplot(3, 3, 5);
imagesc(sino_metal); title("(e) radon of (a) * (d)"); axis('square');

thres = 125;  % manually tuned
for i = 1:size(sino_metal, 1)
    for j = 1:size(sino_metal, 2)
        if sino_metal(i, j) > thres
            sino_metal(i, j) = thres;
        end
    end
end

sino_metal_2 = radon(inv_40_metal .* metal_mask_2, theta);
subplot(3, 3, 8);
imagesc(sino_metal_2); title("(h) radon of (a) * (g)"); axis('square');

thres_2 = 50;  % manually tuned
for i = 1:size(sino_metal_2, 1)
    for j = 1:size(sino_metal_2, 2)
        if sino_metal_2(i, j) > thres_2
            sino_metal_2(i, j) = thres_2;
        end
    end
end

%% Back project metal
% Metal mask 1
H_metal = iradon(sino_metal, theta, "linear", "Hann");
H_metal = H_metal(2:end - 1, 2:end - 1);  % Warping
subplot(3, 3, 6)
imagesc(H_metal, [-0.1, 0.4]);  title("(f) iradon of attenuated (e)"); axis("square");
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end


% Metal mask 2
H_metal_2 = iradon(sino_metal_2, theta, "linear", "Hann");
H_metal_2 = H_metal_2(2:end - 1, 2:end - 1);  % Warping

subplot(3, 3, 9);
imagesc(H_metal_2, [-0.5, 1.0]); title("(i) iradon of attenuated (h)"); axis("square"); 
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end

% Subtraction
H_sub = inv_40_metal;
for i = 1:size(H_sub, 1)
    for j = 1:size(H_sub, 2)
        if metal_mask(i, j) == 0
            H_sub(i, j) = H_sub(i, j) - 0.5 * H_metal(i, j);  % 0.5 is manually tuned. 
        end
    end
end

H_sub_2 = H_sub;
for i = 1:size(H_sub, 1)
    for j = 1:size(H_sub, 2)
        if metal_mask_2(i, j) == 0
            H_sub_2(i, j) = H_sub(i, j) - 0.4 * H_metal_2(i, j);  % 0.4 is manually tuned. 
        end
    end
end

subplot(3, 3, 2);
imagesc(H_sub, [-0.1, 1]); title("(b) Result image: (a) - (d)"); axis("square");
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end

subplot(3, 3, 3);
imagesc(H_sub_2, [-0.1, 1]); title("(c) Result image: (a) - (d) - (h)"); axis("square");
if subimage_switch == 1
    ylim(y_edge); xlim(x_edge);
end
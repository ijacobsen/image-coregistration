clear all; close all; clc;

%% load images
source_path = 'C:\cygwin64\home\ijacobsen\image_corr\working\data_files\';

% lists of filenames
hist_list = ls(strcat(source_path, '*.ndpi'));
SAM_list = ls(strcat(source_path, '*.mat'));

% corresponding images [hist, sam]
ind = [09, 11]; % 4469 5 L  11

% can loop over i if the ind is a list
i = 1;

% load histology image (from nanozoomer file)
filename = strtrim(hist_list(ind(i, 1), :));
NDPIInfo = imfinfo([source_path, filename]);
h_ind = 2;  % which level histology image to use
img = NDPILoad(NDPIInfo, h_ind); 
h_img = img.Data;
h_dx_dy = [1.83, 1.83]; % pixel size in um [dx, dy]

% load SAM images (from Daniel's .mat file)
filename = strtrim(SAM_list(ind(i, 2), :));
SAM_obj = load(strcat(source_path, filename));
s_amp = SAM_obj.sammap.get('amp');
s_amp = uint8(255*(s_amp - min(s_amp(:)))/(max(s_amp(:)) - ...
                   min(s_amp(:))));
s_dx_dy = [2.0, 2.0]; % pixel size in um [dx, dy]

%% register images
% instantiate a registration object
reg = RegClass(s_amp, h_img, s_dx_dy, h_dx_dy);

% do registration
hist = reg.Register();

%% visualize outputs
% load and preprocess parameter images
trim = floor(size(s_amp)/16);
s_amp = s_amp(1:trim(1)*16, 1:trim(2)*16);
s_c = SAM_obj.sammap.get('cM');
s_c = s_c(1:trim(1)*16, 1:trim(2)*16);
s_a = abs(SAM_obj.sammap.get('alphadB'));
s_a = s_a(1:trim(1)*16, 1:trim(2)*16);
s_z = SAM_obj.sammap.get('ZM');
s_z = s_z(1:trim(1)*16, 1:trim(2)*16);

% plot images to match
figure(1)
subplot(2, 1, 1)
imagesc(flipud(h_img))
daspect([1, 1, 1])
xticklabels = 0:2:h_dx_dy(1)*size(h_img, 2)*1e-3;
xticks = 0:size(h_img, 2)/size(xticklabels, 2):size(h_img, 2)-1;
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('mm')
yticklabels = 0:2:h_dx_dy(2)*size(h_img, 1)*1e-3;
xticks = 0:size(h_img, 1)/size(yticklabels, 2):size(h_img, 1)-1;
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'YDir', 'normal')
ylabel('mm')
title('Histology')
subplot(2, 1, 2)
imagesc(flipud(s_amp))
daspect([1, 1, 1])
colormap 'gray'
    set(gca, 'XTick', 0:250:size(s_amp, 2), 'XTickLabel', ...
    0:0.5:2*size(s_amp, 2)*1e-3)
xlabel('mm')
set(gca, 'YTick', 0:250:size(s_amp, 1), 'YTickLabel', ...
    0:0.5:2*size(s_amp, 1)*1e-3, 'YDir', 'normal')
ylabel('mm')
title('SAM Amplitude')
savefig(sprintf('match_reg_iter=%d.fig', i))

% overlay images
figure(3)
N = 11;
sigma = 50;
% design gaussian filter
ind_f = -floor(N/2) : floor(N/2);
[X, Y] = meshgrid(ind_f, ind_f);
h = exp(-(X.^2 + Y.^2)/(2*sigma*sigma));
H = h / sum(h(:));
h = uint8(conv2(double(hist(:, :, 1)), double(H)));
h = h(5:end-6, 5:end-6);
s = s_amp;
s(s > 210) = 255;
s(s < 210) = 15*(log2(double(s(s < 210))));
h = hist(:, :, 1);
h(h > 225) = 255;
h(h < 225) = 15*(log2(double(h(h < 225))));
img = double(cat(3, h, h, s))/255;
imshow(img)
savefig(sprintf('overlay_reg_iter=%d.fig', i))

% show parameter images
figure(2)
subplot(2, 2, 1)
imagesc(flipud(s_c), [1400, 1700])
daspect([1, 1, 1])
set(gca, 'XTick', 0:250:size(s_c, 2), 'XTickLabel', ...
    0:0.5:2*size(s_c, 2)*1e-3)
xlabel('mm')
set(gca, 'YTick', 0:250:size(s_c, 1), 'YTickLabel', ...
    0:0.5:2*size(s_c, 1)*1e-3, 'YDir', 'normal')
ylabel('mm')
cb = colorbar
title(cb, 'm/s')
title('Speed of Sound')
subplot(2, 2, 2)
imagesc(flipud(s_a), [3, 8])
daspect([1, 1, 1])
set(gca, 'XTick', 0:250:size(s_c, 2), 'XTickLabel', ...
    0:0.5:2*size(s_c, 2)*1e-3)
xlabel('mm')
set(gca, 'YTick', 0:250:size(s_c, 1), 'YTickLabel', ...
    0:0.5:2*size(s_c, 1)*1e-3, 'YDir', 'normal')
ylabel('mm')
cb = colorbar
title(cb, 'dB/MHz/cm')
title('Attenuation')
subplot(2, 2, 3)
imagesc(flipud(s_z), [1.4, 1.7])
daspect([1, 1, 1])
set(gca, 'XTick', 0:250:size(s_c, 2), 'XTickLabel', ...
    0:0.5:2*size(s_c, 2)*1e-3)
xlabel('mm')
set(gca, 'YTick', 0:250:size(s_c, 1), 'YTickLabel', ...
    0:0.5:2*size(s_c, 1)*1e-3, 'YDir', 'normal')
ylabel('mm')
cb = colorbar
title(cb, 'Mrayl')
title('Acoustic Impedance')
subplot(2, 2, 4)
imagesc(flipud(hist))
daspect([1, 1, 1])
set(gca, 'XTick', 0:250:size(s_c, 2), 'XTickLabel', ...
    0:0.5:2*size(s_c, 2)*1e-3)
xlabel('mm')
set(gca, 'YTick', 0:250:size(s_c, 1), 'YTickLabel', ...
    0:0.5:2*size(s_c, 1)*1e-3, 'YDir', 'normal')
ylabel('mm')
title('Histology')
savefig(sprintf('parameter_reg_iter=%d.fig', i))

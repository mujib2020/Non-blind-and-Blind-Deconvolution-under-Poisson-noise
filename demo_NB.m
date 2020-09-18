clear; close all; clc;
addpath('.\utilities')

%% Data Generation
% Ground Truth
I = imread('.\data\galaxy256.png');
I = double(I);

% Peak
peak = 25.5;
I = I/max(I(:))*peak;

% Adding blur
H = fspecial('gaussian',9,sqrt(3)); % blurring kernel
I_blurry = imfilter(I,H,'circular');

% Adding Poisson noise
f = poissrnd(I_blurry);

%% TV
% Paremeter setting
pm.beta = 20;
pm.alpha = 1;
pm.mu1 = 0.1;
pm.mu2 = 20;
pm.maxit = 300;

[u_TV,output_TV] = FOTVDeblur_NB(f,H,pm);

%% FOTV
% Paremeter setting
pm.alpha = 1.6; 
pm.beta = 20;
pm.mu1 = 0.1;
pm.mu2 = 20;
pm.maxit = 300;

[u_FOTV,output_FOTV] = FOTVDeblur_NB(f,H,pm);

%% visualize the results
figure;
subplot(221); imshow(f,[0,peak]);title(['Input, PSNR=', num2str(PSNR(I,f))]);
subplot(222); imshow(u_TV,[0,peak]);title(['TV, PSNR=', num2str(PSNR(I,u_TV))]);
subplot(223); imshow(u_FOTV,[0,peak]);title(['FOTV, PSNR=', num2str(PSNR(I,u_FOTV))]);
subplot(224); 
plot(output_FOTV.cpu,output_FOTV.X,'LineWidth',1.2); axis square;  xlabel('CPU Time'); title('Energy');




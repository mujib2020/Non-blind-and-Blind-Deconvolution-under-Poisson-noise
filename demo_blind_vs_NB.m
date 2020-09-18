clear; close all; clc;
addpath('.\utilities')

%% Data Generation
% Ground Truth
I = imread('.\data\sate128.jpg');
I = double(I); 

%Peak
peak = 255; 
I = I/max(I(:))*peak;

% Adding blur
Hg = fspecial('gaussian',9,sqrt(3)); % blurring kernel
[Mh,Nh] = size(Hg); % size of the kernel
I_blurry = conv2(I,Hg,'valid');

% Adding Poisson noise
f = poissrnd(I_blurry);

%% Non-blind deconvolution
% Padding boundary condition for non-blind
f1 = padarray(f,[floor(Mh/2) floor(Nh/2)]);

% % parameter setting
pm.beta = 110;
pm.mu1 = 0.1;
pm.mu2 = 1;
pm.maxit = 50;
pm.alpha = 1;

[u_NB,output_NB] = FOTVDeblur_NB(f1,Hg,pm);

%% Blind deconvolution
% EM
pm.maxit = 20;
[u_EM,H_EM,output_EM] = EM_Blind_Deconv(f,Mh,Nh,pm);

% FOTV 
% % parameter setting
pm.alpha = 1;
pm.mu1 = 1e-1;
pm.mu2 = 1;
pm.maxit = 150;
pm.beta = 110;

[u_blind,H_blind,output_blind] = FOTV_deconv_blind(f,Mh, Nh,pm);


%% visualize the results
figure;
subplot(221); imshow(f1,[0,peak]);title(['Input, PSNR=', num2str(PSNR(I,f1))]);
subplot(222); imshow(u_EM,[0,peak]);title(['EM, PSNR=', num2str(PSNR(I,u_EM))]);
subplot(223); imshow(u_NB,[0,peak]);title(['Non-blind, PSNR=', num2str(PSNR(I,u_NB))]);
subplot(224); imshow(u_blind,[0,peak]);title(['Blind, PSNR=', num2str(PSNR(I,u_blind))]);



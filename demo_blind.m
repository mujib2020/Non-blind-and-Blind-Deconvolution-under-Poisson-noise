clear; close all; clc;
addpath('.\utilities')

%% Data Generation
% Ground Truth
I = imread('.\data\sate128.jpg');
I = double(I);

% Peak
peak = 1e3;
I = I/max(I(:))*peak;

% Adding blur
Hg = fspecial('motion',11,45); %blurring kernel
[Mh,Nh] = size(Hg); % size of the kernel
I_blurry = conv2(I,Hg,'valid');

% Adding Poisson noise
f = poissrnd(I_blurry);

%% EM
pm.maxit = 20;
[u_EM,H_EM,output_EM] = EM_Blind_Deconv(f, Mh, Nh, pm);

%% FOTV
% % parameter setting
% pm.beta = 200;
% pm.alpha = 1.1;
% pm.mu1 = 1e-3;
% pm.mu2 = 0.1;
% pm.maxit = 300;

[u_FOTV,H_FOTV,output] = FOTV_deconv_blind(f,Mh, Nh, pm);

%% visualize the results
f1 = padarray(f,[floor(Mh/2) floor(Nh/2)]);

% plot the blurring kernel
figure;
subplot(131);imshow(Hg,[],'InitialMagnification','fit');title('Original')
subplot(132);imshow(H_EM,[],'InitialMagnification','fit'); title('EM')
subplot(133);imshow(H_FOTV,[],'InitialMagnification','fit'); title('FOTV')

% plot the image
figure;
subplot(221); imshow(f1,[0,peak]);title(['Input, PSNR=', num2str(PSNR(I,f1))]);
subplot(222); imshow(u_EM,[0,peak]);title(['EM, PSNR=', num2str(PSNR(I,u_EM))]);
subplot(223); imshow(u_FOTV,[0,peak]);title(['FOTV, PSNR=', num2str(PSNR(I,u_FOTV))]);
subplot(224); plot(output.cpu,output.X,'LineWidth',1.2); axis square;  xlabel('CPU Time'); title('Energy');




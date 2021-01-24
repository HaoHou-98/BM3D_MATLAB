
clear all
close all

    original_image = im2double(imread('house.png')); %% Input neat image and put it in intensity range [0,1].

    sigma = 15;   %% Standard deviation of the AWGN, one can changes this value, the larger the number, the stronger the noise. 
   
    randn('seed', 0);   %% Generate random number seed.
    
    noisy_image = original_image + (sigma/255)*randn(size(original_image)); %% Create a noisy image. 
 
tic,

[basic_estimation] = BM3D_matlab(noisy_image,sigma); %% BM3D image denoising,basic estimation, hard thresholding stage.

[denoised_image] = BM3D_matlab_wiener(noisy_image,basic_estimation,sigma); %% BM3D image denoising,final estimation, Wiener filtering stage.
% [denoised_image] = BM3D_matlab02(noisy_image,basic_estimation,sigma); %% BM3D image denoising,final estimation, Wiener filtering stage.

toc,

figure,imshow(original_image); title('original neat image');
figure,imshow(noisy_image); title('noisy image');
figure,imshow(basic_estimation); title('basic estimation image');
figure,imshow(denoised_image);title('denoised image');

PSNR = psnr(original_image,denoised_image)  %% PSNR value between original neat image and denoised image by BM3D.
SSIM = ssim(original_image,denoised_image)  %% SSIM value between original neat image and denoised image by BM3D.
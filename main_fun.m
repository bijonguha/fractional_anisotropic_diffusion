display('Code Start')
clear all;
close all;

ig = imread('lena.png');

if ndims(ig)==3
  ig = rgb2gray(ig);
end

im = imnoise(ig,'gaussian',0,0.01);
[psnr(uint8(im),uint8(ig))]

[wA,wH,wV,wD] = dwt2(im,'db4');
wA = im2double(wA);
wH = im2double(wH);
wV = im2double(wV);
wD = im2double(wD);

imgHi = idwt2(wA,wH,zeros(size(wV)),zeros(size(wD)),'db4');
imgH = wave_denoising(imgHi,'db4','s',1);
imgVi = idwt2(wA,zeros(size(wH)),wV,zeros(size(wD)),'db4');
imgV = wave_denoising(imgVi,'db4','s',1);
imgDi = idwt2(wA,zeros(size(wH)),zeros(size(wV)),wD,'db4');
imgD = wave_denoising(imgDi,'db4','s',1);

resALL = medfilt2(fanisodiff(im,14,20,0.30,2,1.1,'AL4'),[5,5]);
resH = medfilt2(fanisodiff(imgH,10,20,0.30,2,1.1,'NS'),[5,1]);
resV = medfilt2(fanisodiff(imgV,10,20,0.30,2,1.1,'EW'),[1,5]);
resD = medfilt2(fanisodiff(imgD,12,20,0.30,2,1,'D'),[5,5]);

img = (resALL+0.4*resH+0.4*resV+0.2*resD)/2;
img = fanisodiff(img,2,20,0.20,1,1,'AL4');

figure()
subplot(2,2,2)
imshow(uint8(img));title('Proposed');
subplot(2,2,3)
imshow(uint8(anisodiff(im,18,25,0.30,2)));title('P-M Model');
subplot(2,2,4)
imshow(uint8(wave_denoising(im,'haar',2)));title('Wavelet Thresholding');
subplot(2,2,1)
imshow(uint8(im));title('Nosiy Image');

[psnr(uint8(img),uint8(ig)),psnr(uint8(anisodiff(im,13,25,0.30,2)),uint8(ig)),...
    psnr(uint8(wave_denoising(im,'haar',2)),uint8(ig))]
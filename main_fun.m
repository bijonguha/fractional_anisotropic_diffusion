display('Code Start')
clear all;
close all;

ig = imread('lena.png');

if ndims(ig)==3
  ig = rgb2gray(ig);
end

imN = imnoise(ig,'gaussian',0,0.03);
[psnr(uint8(imN),uint8(ig))]

im = wave_denoising(imN,'sym4','s',1);

[wA,wH,wV,wD] = dwt2(im,'sym4');
wA = im2double(wA);
wH = im2double(wH);
wV = im2double(wV);
wD = im2double(wD);

imgH = idwt2(wA,wH,zeros(size(wV)),zeros(size(wD)),'sym4');
%imgH = wave_denoising(imgHi,'db4','s',1);
imgV = idwt2(wA,zeros(size(wH)),wV,zeros(size(wD)),'sym4');
%imgV = wave_denoising(imgVi,'db4','s',1);
imgD = idwt2(wA,zeros(size(wH)),zeros(size(wV)),wD,'sym4');
%imgD = wave_denoising(imgDi,'db4','s',1);

resALL = wave_denoising(fanisodiff(im,8,20,0.30,2,1.1,'AL4'),'sym4','s',4);
resH = wave_denoising(fanisodiff(imgH,8,20,0.30,2,1.1,'NS'),'sym4','s',1);
resV = wave_denoising(fanisodiff(imgV,8,20,0.30,2,1.1,'EW'),'sym4','s',1);
resD = wave_denoising(fanisodiff(imgD,5,20,0.30,2,1,'D'),'sym4','s',1);

img = (resALL+0.4*resH+0.4*resV+0.2*resD)/2;
img = fanisodiff(img,5,20,0.20,1,1,'AL4');

figure()
subplot(2,2,2)
imshow(uint8(img));title('Proposed');
subplot(2,2,3)
imshow(uint8(anisodiff(imN,16,25,0.25,2)));title('P-M Model');
subplot(2,2,4)
imshow(uint8(wave_denoising(imN,'db4','s',2)));title('Wavelet Thresholding');
subplot(2,2,1)
imshow(uint8(imN));title('Nosiy Image');

[psnr(uint8(img),uint8(ig)),psnr(uint8(anisodiff(imN,16,25,0.25,2)),uint8(ig)),...
    psnr(uint8(wave_denoising(imN,'db4','s',2)),uint8(ig))]
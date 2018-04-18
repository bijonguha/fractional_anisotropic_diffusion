display('Code Start')
clear all;
close all;

ig = imread('wallpaint.png');
ig = imresize(ig,[1024,1024]);


if ndims(ig)==3
  ig = rgb2gray(ig);
end

ig1 = im2double(ig);
imN = imnoise(ig1,'gaussian',0,0.02);
imN = imN*255;

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

resALL = medfilt2(fanisodiff(im,10,20,0.30,2,1.1,'AL4'),[5,5]);
resH = medfilt2(fanisodiff(imgH,8,20,0.30,2,1.1,'NS'),[5,1]);
resV = medfilt2(fanisodiff(imgV,8,20,0.30,2,1.1,'EW'),[1,5]);
resD = medfilt2(fanisodiff(imgD,8,20,0.30,2,1.1,'D'),[5,5]);

img = (resALL+0.2*resH+0.2*resV+0.6*resD)/2;
%img = fanisodiff(img,5,20,0.20,1,1,'AL4');

andif = anisodiff(imN,13,25,0.25,2);
nshrink = NeighShrinkSUREdenoise(imN);

figure(1)
subplot(2,2,2)
imshow(uint8(img));title('Proposed');
subplot(2,2,3)
imshow(uint8(andif));title('P-M Model');
subplot(2,2,4)
imshow(uint8(wave_denoising(imN,'db4','s',2)));title('Wavelet Thresholding');
subplot(2,2,1)
imshow(uint8(imN));title('Noisy Image');
figure(2)
subplot(2,2,1)
imshow(uint8(nshrink));title('Neighshrink');
['Proposed ||', ' PM Model ||', ' Wavelet thresholding ||', 'NeighShrink']
[psnr(uint8(img),ig),psnr(uint8(andif),ig),...
    psnr(uint8(wave_denoising(imN,'db4','s',2)),ig),...
    psnr(uint8(nshrink),ig)]
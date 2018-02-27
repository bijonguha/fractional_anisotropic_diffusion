%calling anisodiff
clear all;
close all;

'PERONA MALIK simple'
z = imread('lena.jpg');
zn = imnoise(z,'gaussian',0,0.03);
anisodiff(zn,15,25,0.25,2,z);
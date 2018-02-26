%calling anisodiff
clear all;
close all;

'PERONA MALIK simple'
z = imread('body.jpg');
zn = imnoise(z,'gaussian',0,0.01);
anisodiff(zn,10,25,0.25,2,z);
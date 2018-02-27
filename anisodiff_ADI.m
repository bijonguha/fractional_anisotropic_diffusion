% ANISODIFF_mod - Anisotropic diffusion.
%
% 
% diff = fradiff(im, niter, kappa, lambda, option, v)
%
% 
%         im     - input image
%         niter  - number of iterations.
%         kappa  - conduction coefficient 20-100 ?
%         lambda - max value of .25 for stability
%         option - 1 Perona Malik diffusion equation No 1
%                  2 Perona Malik diffusion equation No 2
%         v      - order of fractional derivative ~ [(1,2)]
%
% Return
%         diff   - diffused image.
%
% kappa controls conduction as a function of gradient.  If kappa is low
% then mall intensity gradients are able to block conduction and hence diffusion
% across step edges.  A large value reduces the influence of intensity
% gradients on conduction.
%
% lambda controls speed of diffusion (you usually want it at a maximum of
% 0.25)
%
% Diffusion equation 1 preserve high contrast edges over low contrast ones.
% Diffusion equation 2 favours wide regions over smaller ones.


%function diff = fradiff(im, niter, kappa, lambda, option, v)

close all;
clear all;

%----Manual Running----
im = imread('body.jpg');
im_good = im;

kappa = 25;
lambda = 0.19;
option = 2;
niter = 15;
v = 1.1;
func = 1;
beta = 1;
gv = 0.6;

im = imnoise(im,'gaussian',0,0.03);
psnr(im, im_good)

if ndims(im)==3
  im = rgb2gray(im);
  im_good = rgb2gray(im_good);
end

im = double(im);
[rows,cols] = size(im);
diff = im;


for i = 1:niter
    
 % fprintf('\rIteration %d',i);

  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl = zeros(rows+4, cols+4);
  diffl(3:rows+2, 3:cols+2) = diff;

  UxM1y = diffl(3:rows+2 , 2:cols+1);
  UxM2y = diffl(3:rows+2 , 1:cols  );
  UxP1y = diffl(3:rows+2 , 4:cols+3);
  UxP2y = diffl(3:rows+2 , 5:cols+4);
 
  
  % North, South differences
  %P ~ +(plus), M ~ -(minus)
  deltaN = v*(1-v)/2 * UxP2y + v*UxP1y - diff;
  deltaS = -1*(diff - v * UxM1y - v*(1-v)/2 * UxM2y);
  
  
  %new(-----
  
  diff_b = imgaussfilt(diff, gv);
  
  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl_b = zeros(rows+4, cols+4);
  diffl_b(3:rows+2, 3:cols+2) = diff_b;

  UxM1y_b = diffl_b(3:rows+2 , 2:cols+1);
  UxM2y_b = diffl_b(3:rows+2 , 1:cols  );
  UxP1y_b = diffl_b(3:rows+2 , 4:cols+3);
  UxP2y_b = diffl_b(3:rows+2 , 5:cols+4);
  
  %trial
  %North, South differences
  
  delN = UxP1y_b - diff_b;
  delS = -1*(diff_b - UxM1y_b );
  
  %new)------
  
  % Conduction
  
if func  == 1
  if option == 1
    cN = exp(-(delN/kappa).^2);
    cS = exp(-(delS/kappa).^2);
    
  elseif option == 2
    cN = 1./(1 + (delN/kappa).^2);
    cS = 1./(1 + (delS/kappa).^2);
  end
end

if func == 2
  if option == 1
    cN = exp(-(deltaN/kappa).^2)/2.5;
    cS = exp(-(deltaS/kappa).^2)/2.5;
    
  elseif option == 2
    cN = 1./(1 + (deltaN/kappa).^2)/2.5;
    cS = 1./(1 + (deltaS/kappa).^2)/2.5;
  end
end

% APPLYING FOUR-POINT-TEMPLETE FOR numerical solution of DIFFUSION P.D.E.
  
diff = diff + lambda*(cN.*deltaN + cS.*deltaS );
  
%------------------------------------------
%step 2
  
diffl = zeros(rows+4, cols+4);
diffl(3:rows+2, 3:cols+2) = diff;
  
UxyM1 = diffl(2:rows+1 , 3:cols+2);
UxyM2 = diffl(1:rows   , 3:cols+2);
UxyP1 = diffl(4:rows+3 , 3:cols+2);
UxyP2 = diffl(5:rows+4 , 3:cols+2);
  
%East and West differences
%P ~ +(plus), M ~ -(minus)
  
deltaW = v*(1-v)/2 * UxyP2 + v*UxyP1 - diff;
  
deltaE = -1*(diff - v * UxyM1 - v*(1-v)/2 * UxyM2);
  
%trial
% North, South, East and West differences
  
%new(-----
  
diff_b = imgaussfilt(diff, gv);
  
% Construct diffl which is the same as diff but
% has an extra padding of zeros around it.
diffl_b = zeros(rows+4, cols+4);
diffl(3:rows+2, 3:cols+2) = diff_b;
  
UxyM1_b = diffl(2:rows+1 , 3:cols+2);
UxyM2_b = diffl(1:rows   , 3:cols+2);
UxyP1_b = diffl(4:rows+3 , 3:cols+2);
UxyP2_b = diffl(5:rows+4 , 3:cols+2);
  
  
  
delW = UxyP1_b - diff_b;
delE = -1*(diff_b - UxyM1_b );
  
% Conduction
  
if func  == 1
  if option == 1
    cE = exp(-(delE/kappa).^2);
    cW = exp(-(delW/kappa).^2);
    
  elseif option == 2
    cE = 1./(1 + (delE/kappa).^2);
    cW = 1./(1 + (delW/kappa).^2);
    
  end
end

if func == 2
  if option == 1
    cE = exp(-(delE/kappa).^2)/2.5;
    cW = exp(-(delW/kappa).^2)/2.5;
    
  elseif option == 2
    cE = 1./(1 + (delE/kappa).^2)/2.5;
    cW = 1./(1 + (delW/kappa).^2)/2.5;
  end
end

% APPLYING FOUR-POINT-TEMPLETE FOR numerical solution of DIFFUSION P.D.E.
  
diff = diff + lambda*(cE.*deltaE + cW.*deltaW);
  
diff_q = uint8(diff);
  
diff_q = imhistmatch(diff_q , im_good);
  
[ i, psnr(diff_q, im_good), mse(diff_q, im_good) ]
  
%figure();
% Uncomment the following to see a progression of images
%subplot(ceil(sqrt(niter)),ceil(sqrt(niter)), i)
% figure();
% pause(0)
% imagesc(diff), colormap(gray), axis image

end
% fprintf('\n');
% 
figure()
imhist(cN)
figure()
imhist(cS)
figure()
imhist(cW)
figure()
imhist(cE)


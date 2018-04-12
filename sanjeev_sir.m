ig = imread('cameraman.jpg');
ig = rgb2gray(ig);
[wA1,wH1,wV1,wD1] = dwt2(ig,'db4');

im = imnoise(ig,'gaussian',0,0.02);
niter = 15;
kappa = 25;
option = 2;
lambda = 0.25;

[wA,wH,wV,wD] = dwt2(im,'db4');
wA = double(wA);
wH = double(wH);
wV = double(wV);
wD = double(wD);

[rows,cols] = size(wA);
wA = anisodiff(wA,niter,kappa,lambda,option);

diff_H = wH;
for i = 1:niter
  
  diffl = zeros(rows+2, cols+2);
  diffl(2:rows+1, 2:cols+1) = diff_H;

  % North, South diffusion
  deltaN = diffl(1:rows,2:cols+1)   - diff_H;
  deltaS = diffl(3:rows+2,2:cols+1) - diff_H;
  
  if option == 1
    cN = exp(-(deltaN/kappa).^2);
    cS = exp(-(deltaS/kappa).^2);
    
  elseif option == 2
      
    cN = 1./(1 + (deltaN/kappa).^2);
    cS = 1./(1 + (deltaS/kappa).^2);
     
  end
  
  diff_H = diff_H + lambda*(cN.*deltaN + cS.*deltaS);

  
end

diff_V = wV;

for i = 1:niter

  diffl = zeros(rows+2, cols+2);
  diffl(2:rows+1, 2:cols+1) = diff_V;
  
  
  deltaE = diffl(2:rows+1,3:cols+2) - diff_V;
  deltaW = diffl(2:rows+1,1:cols)   - diff_V;
  
  if option == 1
    cE = exp(-(deltaE/kappa).^2);
    cW = exp(-(deltaW/kappa).^2);
    
  elseif option == 2
      
    cE = 1./(1 + (deltaE/kappa).^2);
    cW = 1./(1 + (deltaW/kappa).^2);
     
  end
  
  diff_V = diff_V + lambda*(cE.*deltaE + cW.*deltaW);
  
end

img = idwt2(wA,diff_H,diff_V,medfilt2(wD,[3,3]),'db4');

figure()
imshow(uint8(img))
figure()
imshow(uint8(im))

psnr(uint8(im),uint8(ig))
psnr(uint8(img),uint8(ig))
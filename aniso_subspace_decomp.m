'New'
ig = imread('lena.jpg');

if ndims(ig)==3
  ig = rgb2gray(ig);
end

[wA1,wH1,wV1,wD1] = dwt2(ig,'db4');

im = imnoise(ig,'gaussian',0,0.01);
[psnr(uint8(im),uint8(ig))]
niter = 25;
kappa = 20;
option = 2;
lambda = 0.50;

[rows,cols] = size(ig);

[wA,wH,wV,wD] = dwt2(im,'db4');
wA = double(wA);
wH = double(wH);
wV = double(wV);
wD = double(wD);
 
% wA = medfilt2(wA,[3,3]);
% wH = medfilt2(wH,[3,1]);
% wV = medfilt2(wV,[1,3]);
% wD = medfilt2(wD,[2,2]);

imgH = idwt2(wA,wH,zeros(size(wV)),zeros(size(wD)),'db4');
imgV = idwt2(wA,zeros(size(wH)),wV,zeros(size(wD)),'db4');
imgD = idwt2(wA,zeros(size(wH)),zeros(size(wV)),wD,'db4');


for j = 1:niter
    diff_H = imgH;
    
for i = 1:j
  
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

diff_V = imgV;

for k = 1:j

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

diff_D = imgD;

for k = 1:j
        dx = 1;
        dy = 1;
        dd = sqrt(2);
        
        hNE = [0 0 1; 0 -1 0; 0 0 0];
        hSE = [0 0 0; 0 -1 0; 0 0 1];
        hSW = [0 0 0; 0 -1 0; 1 0 0];
        hNW = [1 0 0; 0 -1 0; 0 0 0];
    
        deltaNE = imfilter(diff_D,hNE,'conv');
        deltaSE = imfilter(diff_D,hSE,'conv');   
        deltaSW = imfilter(diff_D,hSW,'conv');
        deltaNW = imfilter(diff_D,hNW,'conv'); 
        
        % Diffusion function.
        if option == 1
            cNE = exp(-(deltaNE/kappa).^2);
            cSE = exp(-(deltaSE/kappa).^2);
            cSW = exp(-(deltaSW/kappa).^2);
            cNW = exp(-(deltaNW/kappa).^2);
        elseif option == 2
            cNE = 1./(1 + (deltaNE/kappa).^2);
            cSE = 1./(1 + (deltaSE/kappa).^2);
            cSW = 1./(1 + (deltaSW/kappa).^2);
            cNW = 1./(1 + (deltaNW/kappa).^2);
        end

        % Discrete PDE solution.
        diff_D = diff_D + ...
                  lambda*(...
                  (1/(dd^2))*cNE.*deltaNE + (1/(dd^2))*cSE.*deltaSE + ...
                  (1/(dd^2))*cSW.*deltaSW + (1/(dd^2))*cNW.*deltaNW );

end

img_Rec = (0.8*diff_H +0.8* diff_V +1.4* diff_D)/3;
[psnr(uint8(anisodiff(im,j,kappa,lambda,option)),uint8(ig)),psnr(uint8(img_Rec),uint8(ig))] 

end

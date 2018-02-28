close all;
clear all;

'Fractional perona malik'

im = imread('boat.jpg');
im_good = im;

kappa = 25;
lambda = 0.19;
option = 2;
niter = 10;
v = 1.1;
func = 1;
beta = 1;
gv = 0.6;
mean = 0;
var = 0.01;

KO = 25; K1 = 5;

im = imnoise(im,'gaussian',mean,var);
psnr(im, im_good)

if ndims(im)==3
  im = rgb2gray(im);
  im_good = rgb2gray(im_good);
end

im = double(im);
[rows,cols] = size(im);
diff = im;

for i = 1:niter
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
    
% delN = v*(1-v)/2 * UxP2y_b + v*UxP1y_b - diff_b;
% delS = -1*(diff_b - v * UxM1y_b - v*(1-v)/2 * UxM2y_b);

  if func == 3
      % MAD implementation
      deln_med = medfilt2(delN,[3,3],'symmetric');       
      deln_abs = abs(delN - deln_med);
      delN_M = 1.4826* medfilt2(deln_abs,[3,3],'symmetric') * sqrt(5);

      dels_med = medfilt2(delS,[3,3],'symmetric');
      dels_abs = abs(delS - dels_med);
      delS_M = 1.4826* medfilt2(dels_abs,[3,3],'symmetric') * sqrt(5);
  end
  
%--Under testing
if func == 4
% Median K implementation
% North

dnp = padarray(delN, [1,1]);
delNmax = delN;
delNmin = delN;
delNmed = medfilt2(delN,[3,3],'symmetric');

for ij=2:rows+1
for j=2:cols+1
    maxn = dnp(ij,j);
    minn = dnp(ij,j);
    for k=-1:1
        for l=-1:1
            if dnp(ij+k,j+l) >= maxn
                maxn = dnp(ij+k,j+l);
            end
            if dnp(ij+k,j+l) <= minn
                minn = dnp(ij+k,j+l);
            end
        end
    end
    delNmax(ij-1,j-1) = maxn;
    delNmin(ij-1,j-1) = minn;
end
end

delNd = delNmax - delNmin;
delNn = 1+abs(delNmed - delN);

%South
dns = padarray(delS, [1,1]);
delSmax = delS;
delSmin = delS;
delSmed = medfilt2(delS,[3,3],'symmetric');

for ik=2:rows+1
for j=2:cols+1
        maxs = dns(ik,j);
        mins = dns(ik,j);
        for k=-1:1
            for l=-1:1
                if dns(ik+k,j+l) >= maxs
                    maxs = dns(ik+k,j+l);
                end
                if dns(ik+k,j+l) <= mins
                    mins = dns(ik+k,j+l);
                end
            end
        end
        delSmax(ik-1,j-1) = maxs;
        delSmin(ik-1,j-1) = mins;
end
end

delSd = delSmax - delSmin;
delSn = 1+abs(delSmed - delS);

end
%-----
  
  
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
  
  if func == 3
      indN  = ( abs(delN) < delN_M );
      indS  = ( abs(delS) < delS_M );
    
      cN = zeros(size(delN));
      cS = zeros(size(delS));
    
      cN(indN) =  1*(1 - ((delN(indN)./delN_M(indN)).^2)).^2;
      cS(indS) =  1*(1 - ((delS(indS)./delS_M(indS)).^2)).^2;
  end
  
    
  if func == 4
      
      cN = 1./(1 + (delN./(KO + K1 .* exp(-delNn./delNd))).^2);
      cS = 1./(1 + (delS./(KO + K1 .* exp(-delSn./delSd))).^2);
      
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
  
  % delW = v*(1-v)/2 * UxyP2_b + v*UxyP1_b - diff_b;
  % delE = -1*(diff_b - v * UxyM1_b - v*(1-v)/2 * UxyM2_b);
  
  if func == 3

      % MAD implementation
      delw_med = medfilt2(delW,[3,3],'symmetric');
      delw_abs = abs(delW - delw_med);
      delW_M = 1.4826* medfilt2(delw_abs,[3,3],'symmetric') * sqrt(5);

      dele_med = medfilt2(delE,[3,3],'symmetric');
      dele_abs = abs(delE - dele_med);
      delE_M = 4.4826* medfilt2(dele_abs,[3,3],'symmetric') * sqrt(5);
   
  end
%-----Under testing

if func ==4
% Median K implementation
% West
dnw = padarray(delW, [1,1]);
delWmax = delW;
delWmin = delW;
delWmed = medfilt2(delW,[3,3],'symmetric');

for il=2:rows+1
for j=2:cols+1
    maxw = dnw(il,j);
    minw = dnw(il,j);
    for k=-1:1
        for l=-1:1
            if dnw(il+k,j+l) >= maxw
                maxw = dnw(il+k,j+l);
            end
            if dnw(il+k,j+l) <= minw
                minw = dnw(il+k,j+l);
            end
        end
    end
    delWmax(il-1,j-1) = maxw;
    delWmin(il-1,j-1) = minw;
end
end

delWd = delWmax - delWmin;
delWn = 1+abs(delWmed - delW);

%East
dne = padarray(delE, [1,1]);
delEmax = delE;
delEmin = delE;
delEmed = medfilt2(delE,[3,3],'symmetric');

for im=2:rows+1
    for j=2:cols+1
        maxe = dne(im,j);
        mine = dne(im,j);
        for k=-1:1
            for l=-1:1
                if dne(im+k,j+l) >= maxe
                    maxe = dne(im+k,j+l);
                end
                if dne(im+k,j+l) <= mine
                    mine = dne(im+k,j+l);
                end
            end
        end
        delEmax(im-1,j-1) = maxe;
        delEmin(im-1,j-1) = mine;
    end
end

delEd = delEmax - delEmin;
delEn = 1+abs(delEmed - delE);

end
%--------

  
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
  
  
  if func == 3
      indW  = ( abs(delW) < delW_M );
      indE  = ( abs(delE) < delE_M );
    
      cW = zeros(size(delW));
      cE = zeros(size(delE));
    
      cW(indW) = 1 * (1 - ((delW(indW)./delW_M(indW)).^2)).^2;
      cE(indE) = 1 * (1 - ((delE(indE)./delE_M(indE)).^2)).^2;
  end
  
  if func == 4
      
      cW = 1./(1 + (delW./(KO + K1 .* exp(-delWn./delWd))).^2);
      cE = 1./(1 + (delE./(KO + K1 .* exp(-delEn./delEd))).^2);
  
  end

  
  % APPLYING Alternate Direction now in X, FOR numerical solution of 
  % DIFFUSION P.D.E.
  
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

figure()
imhist(cN)
figure()
imhist(cS)
% figure()
% imhist(cW)
% figure()
% imhist(cE)
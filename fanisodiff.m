function f_diff = fanisodiff(im, niter, kappa, lambda, option,v,dir)

%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               NUM_ITER - number of iterations. 
%               DELTA_T  - integration constant (0 <= delta_t <= 1/7).
%                          Usually, due to numerical stability this 
%                          parameter is set to its maximum value.
%               KAPPA    - gradient modulus threshold that controls the conduction.
%               OPTION   - conduction coefficient functions proposed by Perona & Malik:
%                          1 - c(x,y,t) = exp(-(nablaI/kappa).^2),
%                              privileges high-contrast edges over low-contrast ones. 
%                          2 - c(x,y,t) = 1./(1 + (nablaI/kappa).^2),
%                              privileges wide regions over smaller ones. 
%               DIR      - Direction in which anisotropic diffusion has to
%                          be performed
%                          1-['AL4']- All four 2-['EW']-East-West 3-['NS']-North-South 4-['D']-Both
%                          diagonals 5- Individual ['E','W','N','S']
%               V        - Order of fractional derivative( [-1,-2) )
%
%       OUTPUT DESCRIPTION:
%                DIFF_IM - (diffused) image with the largest scale-space parameter.


% Convert input image to double.
im = double(im);

% PDE (partial differential equation) initial condition.
f_diff = im;

% Center pixel distances.
dx = 1;
dy = 1;
dd = sqrt(2);

% 2D convolution masks - integral finite differences.

hN = [0 1 0; 0 -1 0; 0 0 0];
hS = [0 0 0; 0 -1 0; 0 1 0];
hE = [0 0 0; 0 -1 1; 0 0 0];
hW = [0 0 0; 1 -1 0; 0 0 0];
hNE = [0 0 1; 0 -1 0; 0 0 0];
hSE = [0 0 0; 0 -1 0; 0 0 1];
hSW = [0 0 0; 0 -1 0; 1 0 0];
hNW = [1 0 0; 0 -1 0; 0 0 0];

% 2D convolution masks - fractional finite differences.

om_2 = v*(v-1)/2;
om_1 = v;
normalize = ((om_1+om_2+1)/2)^10;

hNf = zeros(5,5);
hSf = zeros(5,5);
hEf = zeros(5,5);
hWf = zeros(5,5);

hNf(3,3) =-1; hNf(2,3)=om_1; hNf(1,3) = om_2;
hSf(3,3) =-1; hSf(4,3)=om_1; hSf(5,3) = om_2;
hEf(3,3) =-1; hEf(3,4)=om_1; hEf(3,5) = om_2;
hWf(3,3) =-1; hWf(3,2)=om_1; hEf(3,1) = om_2;
 

switch dir
    case 'AL4'
        for t = 1:2
            
        nablaN = imfilter(f_diff,hNf,'conv');
        nablaS = imfilter(f_diff,hSf,'conv');   
        nablaW = imfilter(f_diff,hWf,'conv');
        nablaE = imfilter(f_diff,hEf,'conv');   

        if option == 1
            cN = exp(-(nablaN/(normalize*kappa)).^2);
            cS = exp(-(nablaS/(normalize*kappa)).^2);
            cW = exp(-(nablaW/(normalize*kappa)).^2);
            cE = exp(-(nablaE/(normalize*kappa)).^2);

        elseif option == 2
            cN = 1./(1 + (nablaN/(1/normalize*kappa)).^2);
            cS = 1./(1 + (nablaS/(1/normalize*kappa)).^2);
            cW = 1./(1 + (nablaW/(1/normalize*kappa)).^2);
            cE = 1./(1 + (nablaE/(1/normalize*kappa)).^2);

        end

        % Discrete PDE solution.

        f_diff = f_diff + ...
                  lambda*(...
                  (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
                  (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE);
 
        end
        if niter > 2
            for t = 2:niter
                
                nablaN = imfilter(f_diff,hN,'conv');
                nablaS = imfilter(f_diff,hS,'conv');
                nablaW = imfilter(f_diff,hW,'conv');
                nablaE = imfilter(f_diff,hE,'conv');
                
                if option == 1
                    cN = exp(-(nablaN/kappa).^2);
                    cS = exp(-(nablaS/kappa).^2);
                    cW = exp(-(nablaW/kappa).^2);
                    cE = exp(-(nablaE/kappa).^2);
                    
                elseif option == 2
                    cN = 1./(1 + (nablaN/kappa).^2);
                    cS = 1./(1 + (nablaS/kappa).^2);
                    cW = 1./(1 + (nablaW/kappa).^2);
                    cE = 1./(1 + (nablaE/kappa).^2);
                    
                end
                
                % Discrete PDE solution.
                
                f_diff = f_diff + ...
                    lambda*(...
                    (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
                    (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE);
                
            end
        end
        
    case 'EW'
        for t = 1:2
        nablaW = imfilter(f_diff,hWf,'conv');
        nablaE = imfilter(f_diff,hEf,'conv');
         
         if option == 1
            cW = exp(-(nablaW/(normalize*kappa)).^2);
            cE = exp(-(nablaE/(normalize*kappa)).^2);
         elseif option == 2
            cW = 1./(1 + (nablaW/(1/normalize*kappa)).^2);
            cE = 1./(1 + (nablaE/(1/normalize*kappa)).^2);
         end
         
         f_diff = f_diff + ...
                  lambda*((1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE );
        end
        
        if niter > 2
            for t = 2:niter
                nablaW = imfilter(f_diff,hW,'conv');
                nablaE = imfilter(f_diff,hE,'conv');
                
                if option == 1
                    cW = exp(-(nablaW/kappa).^2);
                    cE = exp(-(nablaE/kappa).^2);
                elseif option == 2
                    cW = 1./(1 + (nablaW/kappa).^2);
                    cE = 1./(1 + (nablaE/kappa).^2);
                end
                
                f_diff = f_diff + ...
                    lambda*((1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE );
            end
        end
        
    case 'NS'
        
        for t = 1:2
        nablaN = imfilter(f_diff,hNf,'conv');
        nablaS = imfilter(f_diff,hSf,'conv');
        
         if option == 1
          cN = exp(-(nablaN/(normalize*kappa)).^2);
          cS = exp(-(nablaS/(normalize*kappa)).^2);
         elseif option == 2
          cN = 1./(1 + (nablaN/(1/normalize*kappa)).^2);
          cS = 1./(1 + (nablaS/(1/normalize*kappa)).^2);
         end
        
         f_diff = f_diff + ...
                  lambda*((1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS);
         
        end
        if niter > 2
            for t = 2:niter
                nablaN = imfilter(f_diff,hN,'conv');
                nablaS = imfilter(f_diff,hS,'conv');
                
                if option == 1
                    cS = exp(-(nablaS/kappa).^2);
                    cN = exp(-(nablaN/kappa).^2);
                elseif option == 2
                    cS = 1./(1 + (nablaS/kappa).^2);
                    cN = 1./(1 + (nablaN/kappa).^2);
                end
                
                f_diff = f_diff + ...
                    lambda*((1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS);
                
            end
        end
        
    case 'D'
        
        for t=1:niter
         nablaNE = imfilter(f_diff,hNE,'conv');
         nablaSE = imfilter(f_diff,hSE,'conv');   
         nablaSW = imfilter(f_diff,hSW,'conv');
         nablaNW = imfilter(f_diff,hNW,'conv');
         
         if option == 1
            cNE = exp(-(nablaNE/kappa).^2);
            cSE = exp(-(nablaSE/kappa).^2);
            cSW = exp(-(nablaSW/kappa).^2);
            cNW = exp(-(nablaNW/kappa).^2);
            
         elseif option == 2
            cNE = 1./(1 + (nablaNE/kappa).^2);
            cSE = 1./(1 + (nablaSE/kappa).^2);
            cSW = 1./(1 + (nablaSW/kappa).^2);
            cNW = 1./(1 + (nablaNW/kappa).^2);
         
         end
         
         f_diff = f_diff + ...
                  lambda*((1/(dd^2))*cNE.*nablaNE + ... 
                  (1/(dd^2))*cSE.*nablaSE + ...
                  (1/(dd^2))*cSW.*nablaSW + (1/(dd^2))*cNW.*nablaNW );
              
        end
    
end
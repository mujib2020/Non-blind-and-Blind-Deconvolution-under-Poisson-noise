function [u,output] = FOTVDeblur_NB(f,H,pm)
% Fractional-order TV regularized method

% If you use this code or part of this code please cite the following references:
%
% REFERENCES: 1) Chowdhury, M. R. and Qin, J. and Lou, Y.; Non-blind and Blind Deconvolution under Poisson Noise using Fractional-order Total Variation
%                Journal of Mathematical Imaging and Vision; 2020
%             2) Chowdhury, M. R. and Zhang, J. and Qin, J. and Lou, Y.; Poisson image denoising based on fractional-order total variation
%                Inverse Problem and Imaging Volume 14, No. 1, 2020, 77-96
%                doi: 10.3934/ipi.2019064

%     Input :
%             f                   : noisy-blurry image
%             H                   : kernel/psf
%             pm.beta             : balancing parameter for regularization and data fitting term
%             pm.mu1 and pm.mu2   : penalty parameter
%             pm.alpha            : order of the derivative
%             pm.maxit            : maximum iteration number
%
%     Output:
%             u          :    recovered/denoised image
%             output.X   :    energy
%             output.Err :    error between two consecutive iterates
%             output.cpu :    computation time


%% parameter setting
alpha = 1.6; 
beta = 1;
mu1 = 0.1;
mu2 = 0.1;
maxit = 300;

if isfield(pm,'alpha'); alpha = pm.alpha; end
if isfield(pm,'beta'); beta = pm.beta; end
if isfield(pm,'mu1'); mu1 = pm.mu1; end
if isfield(pm,'mu2'); mu2 = pm.mu2; end
if isfield(pm,'maxit'); maxit = pm.maxit; end

%% initizing image u
u = f;
[m,n] = size(u);

%% initialization of the auxiliary and dual variables
lam_1x = zeros(m,n);
lam_1y = lam_1x;
lam2 = lam_1x;
z1 =lam_1x;
z2 = lam_1x;

%% initialization of energy and break
output.X = [];
output.Err = [];

%% FOTV setup
K = 15; w = zeros(1,K); w(K) = 1;
for i = 1:(K-1)
    w(K-i) = (-1)^i*gamma(alpha+1)/gamma(alpha-i+1)/factorial(i);
end
[D,Dt] = defDDt(w,alpha); % compute fractional order derivatives

% for the denominator of u
D1 = abs(psf2otf(w,[m,n])).^2; % column
D2 = abs(psf2otf(w',[m,n])).^2; % row
C = D1+D2;

%% Deconvolution
U = u;

% Find the convolution matrix A
A = psf2otf(H,[m n]);
AAt = A.^2; % A'*A
g = real(ifft2(A.*fft2(u)));

tstart = tic;
for i = 1:maxit
    
    %%%%%%%%%% u-subproblem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % terms of numerator
    u1 = mu2*real(ifft2(conj(A).*fft2(g)));
    u2 = real(ifft2(conj(A).*fft2(lam2)));
    u3 = Dt(lam_1x,lam_1y);
    u4 = mu1*Dt(z1,z2);
    
    % terms of denominator
    u5 = mu2 * AAt;
    uu = fft2(u1+u2+u3+u4)./(mu1*C + u5);
    
    u = real(ifft2(uu));
    
    
    %%%%%%%%%%% z - subprblem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Dux,Duy]= D(u);
    Zx = Dux-lam_1x/mu1;
    Zy = Duy-lam_1y/mu1;
    
    W = sqrt (Zx.^2+Zy.^2);
    Max = max(W - 1/mu1,0);
    W (W == 0) = 1;
    
    z1 = Max.*Zx./W;
    z2 = Max.*Zy./W;
        
    
    %%%%%%%%%%%% g-subproblem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1st term
    g1 = real(ifft2(A.*fft2(u)));
    g2 = (beta+lam2)/mu2;
    g3 = (g1-g2)/2;
    
    g = g3 + sqrt (g3.^2 + beta * f /mu2);
    
    
    %%%%%%%%%%%%Lambda update%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lam_1x = lam_1x + mu1*(z1-Dux);
    lam_1y = lam_1y + mu1*(z2-Duy);
    
    lam2 = lam2+ mu2 * (g-g1); % q1 = Au
    
    output.cpu(i) = toc(tstart);
    
    %%%%%%%%%%%%Energy curve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    Unorm = sqrt(Dux.^2+Duy.^2);
    
    xTV = sum(sum(Unorm));
    xData = sum(sum((  real(ifft2(A.*fft2(u)))  -f.*log( abs( real(ifft2(A.*fft2(u)))  ))) ));
    
    x1 = xTV + beta * xData;
    output.X = [output.X x1];
    
    
    %%%%%%%%%%%%Error curve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    err = norm(u(:)-U(:))/norm(U(:));
    output.Err = [output.Err err];
    
    if err <= 1e-5
        break;
    end
    U = u;
    
end


end

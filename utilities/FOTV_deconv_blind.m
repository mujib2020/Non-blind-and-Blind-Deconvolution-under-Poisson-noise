function [u,H,output] = FOTV_deconv_blind(f,Mh, Nh,pm)
% Fractional-order TV regularized method

% If you use this code or part of this code please cite the following references:
%
% REFERENCES: 1) Chowdhury, M. R. and Qin, J. and Lou, Y.; Non-blind and Blind Deconvolution under Poisson Noise using Fractional-order Total Variation
%                Journal of Mathematical Imaging and Vision; 2020
%             2) Chowdhury, M. R. and Zhang, J. and Qin, J. and Lou, Y.; Poisson image denoising based on fractional-order total variation
%                Inverse Problem and Imaging Volume 14, No. 1, 2020, 77-96
%                doi: 10.3934/ipi.2019064
%
%     Input :
%             f                   : noisy-blurry image
%             pm.beta             : balancing parameter for regularization and data fitting term
%             pm.mu1 and pm.mu2   : penalty parameter
%             pm.alpha            : order of the derivative
%             pm.maxit            : maximum iteration number
%
%     Output:
%             u          :    recovered/denoised image
%             H          :    recovered kernel/psf
%             output.X   :    energy
%             output.Err :    error between two consecutive iterates
%             output.cpu :    computation time


%% parameter setting
alpha = 1;
beta = 1;
mu1 = 0.1;
mu2 = 0.1;
maxit = 100;

if isfield(pm,'alpha'); alpha = pm.alpha; end
if isfield(pm,'beta'); beta = pm.beta; end
if isfield(pm,'mu1'); mu1 = pm.mu1; end
if isfield(pm,'mu2'); mu2 = pm.mu2; end
if isfield(pm,'maxit'); maxit = pm.maxit; end

%% initizing image u and kernerl H
u = padarray(f,[floor(Mh/2) floor(Nh/2)],'replicate');
H = ones(Mh,Nh)/Mh/Nh;
[m,n] = size(u);

%% initialization of the auxiliary and dual variables
Z = zeros(m,n);
lam_1x = zeros(m,n);
lam_1y = lam_1x;
lam2 = lam_1x;
u_const = ones(size(u));

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

tstart = tic;
for  it = 1:maxit
    
    %%%%%%%%%% update sharp image%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%% EM step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i2 = 1:1
        f1a = u./conv2( u_const,rot90(H,2), 'same');
        
        f1b_1 = conv2(u, H, 'valid');
        f1b_2 = f./max(1e-3,f1b_1);
        f1b = conv2(f1b_2, rot90(H,2), 'full') ;
        f1 = f1a.*f1b; % u_half+k
    end
    
    
    %%%%%%%%%% FOTV step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%% u-subproblem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    c = conv2( u_const,rot90(H,2), 'same');
    u1 = beta*c - mu2*Z-lam2;
    u = -u1/(2*mu2) + (sqrt(u1.^2 + 4*mu2*beta*f1.*c))./(2*mu2);
    
    
    %%%%%%%%%% v-subproblem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Dux,Duy]= D(u);
    Vx = Dux-lam_1x/mu1;
    Vy = Duy-lam_1y/mu1;
    W = sqrt (Vx.^2+Vy.^2);
    Max = max(W - 1/mu1,0);
    W (W == 0) = 1;
    
    v1 = Max.*Vx./W;
    v2 = Max.*Vy./W;
        
    
    %%%%%%%%%% Z-subproblem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    z1 = Dt(lam_1x,lam_1y);
    z2 = mu1*Dt(v1,v2);
    z3 = mu2 * u - lam2;
    
    z4 = mu1*C + mu2;
    
    ZZ = fft2(z1+z2+z3)./z4;
    Z = real(ifft2(ZZ));
    
    
    %%%%%%%%%%%%Lambda update%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lam_1x = lam_1x + mu1*(v1-Dux);
    lam_1y = lam_1y + mu1*(v2-Duy);
    lam2 = lam2 + mu2*(Z-u);
    
       
    %%%%%%%%%%%%update blur%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gradk  = zeros(Mh,Nh);
    h1 = H./conv2fft(rot90(u, 2),u_const,'valid');
    h2 = conv2fft(rot90(u, 2), f./max(conv2( u,H, 'valid'),1e-3),'valid');
    H = gradk + h1.*h2;
    
    % Kernel projection
    H = H.*(H>0);
    H = H/sum(H(:));
    
    
    output.cpu(it) = toc(tstart);
    
    %%%%%%%%%%%%Energy curve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [Dux,Duy]= D(u);
    Unorm = sqrt(Dux.^2+Duy.^2);
    
    xTv = sum(sum( Unorm ));
    
    Data_1 = conv2(u,H,'valid');
    xData = beta*sum(sum( Data_1- f.* log(abs(Data_1)+1e-14)));
    
    x1 = xTv + xData;
    output.X = [output.X x1];
    
    
    %%%%%%%%%%%%Error curve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    err = norm(u(:)-U(:))/norm(U(:));
    output.Err = [output.Err err];
    
    if err <= 1e-5
        break;
    end
    U = u;
    
end






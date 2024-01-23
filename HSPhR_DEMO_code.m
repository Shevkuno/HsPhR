%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Demo-code for hyperspectral quantitative phase 
%  imaging from spectrally multiplexed observations,
%  I. Shevkunov, V. Katkovnik,  K.Eguazarian, 2023
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2020-2023 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only
% be used for nonprofit noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes
% is prohibited.
%
% AUTHORS:
%     I. Shevkunov, V. Katkovnik, K.Eguazarian
%     email: igor.shevkunov@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

addpath('.\AUX files')
fprintf ('\n Initialization started \n')
tic
%% ---------General parameters -------------
n1 = 64;                         % image size in dimension 1
n2 = 64;                         % image size in dimension 2
distance = 0.002;                % propagation distance between object and sensor
dx = 3.45e-6;                    % pixel size of sensor
K = 6;                           % number of wavelengths
T = 18;                          % number of masks and experiments
nn = 150;                        % number of iterations

%% --------- Options for running reconstruction ------------------
do_CCF_filtering=1;              % 1 - CCF filtering is on; 0 - off
do_Lagrange=1;                   % 1 - Lagrange variable are calculated; 0 - not
betta_CCF=0.3;                   % reconstruction parameter
%% Effects of Lagrange Multipliers and CCF: demonstration
it_start1=10;                    % iteration number to turn  on Lagrange variables 
it_start2=50;                    % iteration number to turn on CCF  

%% ------------ Wavelength spectrum-----------------------
lambda1 = 550e-9;                % The first wavelength of the range
lambda_end = 950e-9;%            % the last wavelength of the range

delta_lambda = (lambda_end-lambda1)/(K-1);
lambda_set = lambda1:delta_lambda:lambda_end; % wavelengths set

%% ------------ Refractive index creation --------------------------------
glass_type=1;                    % 1 is for BK7 glass
n_ref = Refraction_Cauchy(lambda_set,glass_type); %refractive index according to Cauchy's equation

%% ------------ Noise parameters --------------------------------------------
do_poisson_noise = 0;            % 1 - is for Poisson noise; 0 - Gaussian noise
sigma = 0.001;                   % standard deviation for Gaussian noise
gamma_gauss = 20*sigma^2;        % reconstruction parameter
noise_type = 'Gaussian noise';   % prefix for figures

if do_poisson_noise
   noise_type = 'Poisson noise';      % prefix for figures
   kappa = 170000;                    % Poisson noise parameter
   gamma_poisson = 6/kappa;           % reconstruction parameter
end
%% Make Hyperspectral complex-valued object
%--------(1) upload images for amplitude and phase ------------------------
x_phase = double(imread('image_Cameraman64.bmp'))/255; % upload cameraman 64x64 image
x_ampl = double((imread('peppers_gray.png')))/255;     % upload peppers image
x_ampl = x_ampl(211+(1:n1),246+(1:n2));                % crop 64x64 part from peppers

%--------(2) create phase shifting properties of the object w.r.t. wavelenghts and refractive index
s_lambda = 0;
for lambda = lambda_set
    s_lambda = s_lambda+1;
    varphi(:,:,s_lambda) = x_phase*lambda1*2*pi/lambda*(n_ref(s_lambda)-1); % lambda1 for scaling
    fprintf ('.')
end
%--------(3) Make a Hyperspectral object ----------------------------------
x = x_ampl.*exp(1j*varphi);

%--------(4) Make zero-padding for angular spectrum proper propagation ---------------------
Nzp  =  ceil(max(lambda_set)*distance/(dx)^2);   % zero-padding size
if rem(Nzp,2)~= 0;Nzp = Nzp+1;end % zero-padding size needed to be even
x = padarray(x,[Nzp/2 Nzp/2]);
object_coord_1=Nzp/2+(1:n1);
object_coord_2=Nzp/2+(1:n2);
%% Masks creation
Masks_set = zeros(Nzp+n1,Nzp+n2,K,T,'single');       % preallocation
s_lambda = 0;
for lambda = lambda_set
    s_lambda = s_lambda+1;
    rng(44)
    for t = 1:T %% Masks depending on wavelength and refractive index
        temp = randsrc(n1,n2,[0 pi/2 -pi/2 pi/4  -pi/4])/lambda*lambda1*(n_ref(s_lambda)-1);% lambda1 for scaling
        temp = exp(1j*temp);
        Masks_set(:,:,s_lambda,t) = padarray(temp,[round(Nzp/2) round(Nzp/2)],1); %used masks set
        fprintf ('.')
    end
end
%% Angular Spectrum (AS) Propagation operator
%--------- (1) AS Transfer Function ------------------------------
[LL,KK,~] = size(x);
k = single(-KK/2:KK/2 - 1);
l = single(-LL/2:LL/2 - 1);
[k,l] = meshgrid(k,l);
AS = zeros(size(x));
s_lambda = 0;
for lambda = lambda_set
    s_lambda = s_lambda+1;
    shift = 2*pi/lambda*distance;
    U = 1 - lambda^2*((k/(dx*KK)).^2+(l/(dx*LL)).^2);
    TFunc = exp(1i*2*pi/lambda*distance*sqrt(U)); TFunc(U<0) = 0;
    AS(:,:,s_lambda) = (TFunc)*exp(-1j*shift);   %% AS Trunsfer function
end
clear k l U
%--------- (2) hyperspectral propagation with masks included ----------------
A = @(wf,Masks) ifft2(arrayIshift(AS.*(arrayshift(fft2(conj(Masks).*wf))))); % forward propagation
At = @(wf, Masks) Masks.*ifft2(arrayIshift(conj(AS).*arrayshift(fft2(wf)))); % backward propagation
%% Observations model
%-------preallocations-----------------------------------------------------
B = zeros(LL,KK,K,T,'single');
BB = zeros(LL,KK,K,T,'single');
D = zeros(LL,KK,K,T,'single');
Y = zeros(LL,KK,T,'single');
Z = zeros(LL,KK,T,'single');
%-------observations creation ---------------------------------------------
for t = 1:T %create observations for each mask
    Masks = (Masks_set(:,:,:,t));
    B(:,:,:,t) = A(x,Masks);                         % HS wavefront propagated to sensor
    Y(:,:,t) = sum((squeeze(abs(B(:,:,:,t)))).^2,3); % noiseless obesrvation summed along wavelengths
    
    if do_poisson_noise                              % Poisson noise Intensity Observations
        Z(:,:,t)  = poissrnd(Y(:,:,t)*kappa);
         SNR(t) = 10*log10(sum(kappa^2*Y(:,:,t).^2)/sum((Y(:,:,t)*kappa-Z(:,:,t)).^2));
        
    else                                             % Gaussian noise Intensity Observations
        Z(:,:,t) = Y(:,:,t)+randn(size(squeeze(Y(:,:,t))))*sigma;
        SNR(t) = 10*log10(sum(Y(:,:,t).^2)/sum((Y(:,:,t)-Z(:,:,t)).^2));
    end
    fprintf ('.')
end
%% Reconstruction Algorithm
%------------(1) Random Initialization ------------------------------------
x_phase = rand(size(x));
s_lambda = 0;
for lambda = lambda_set
    s_lambda = s_lambda+1;
    varphi_0(:,:,s_lambda) = x_phase(s_lambda)*2*pi/lambda*lambda1*(n_ref(s_lambda)-1);
end
xs = (rand(size(x))).*exp(1j*randn(size(x)));

%% ------------- figures section --------------------------------------------
fig1 = figure('Name','Relative errors','units','normalized','outerposition',...
    [0.1 0.4 0.4 0.35],'color', 'w');
Logo=axes('position',[0, 0, .12,.12,]);[logo_im, ~]=imread('Tuni logo.png');image(logo_im); set(Logo,'handlevisibility','off','visible','off'),

lambda_to_show = ceil(K/2); % wavelength number to show wavefront images for ampl and phase 
fig2 = figure('Name','Reconstructed amplitude and phase','units',...
    'normalized','outerposition',[0.5 0.4 0.22 0.4],'color', 'w');
Logo=axes('position',[0, 0, .21,.105,]);[logo_im, ~]=imread('Tuni logo.png');image(logo_im); set(Logo,'handlevisibility','off','visible','off'),

subplot 221, imshow(abs(x(Nzp/2+(1:n1),Nzp/2+(1:n2),lambda_to_show)),[]), ...
    title('True amplitude'), ...
    c = colorbar;
c.Label.String = 'Amplitude, a.u.';
subplot 222, imshow(varphi(:,:,lambda_to_show),[]),
title('True phase'),
c = colorbar;
c.Label.String = 'Phase, rad';
%--------------------------------------------------------------------------
fprintf ('\n Initialization took %g seconds', toc)
fprintf ('\n ------------------------------ \n Start of reconstruction.')
E = zeros(size(Z),'single');  %preallocation

%% Main Iterations
for s = 1:nn 
tic 
%% ---------------------(2) forward propagation--------------------------
    for t = 1:T
        B(:,:,:,t) = A(xs,Masks_set(:,:,:,t)); % HS object propagated to the sensor
        Y(:,:,t) = sum((abs(B(:,:,:,t))).^2,3);
        
        if do_poisson_noise
            E(:,:,t) = Z(:,:,t)- kappa*Y(:,:,t);  % residual estimation
        else
            E(:,:,t) = Z(:,:,t)- Y(:,:,t);              % residual estimation
        end
    end
    
    Bold = B; % temporal variable for saving B for the next iteration
    %
    clear E_solution
    %% ------------------(3) Update B by the proximal operators--------------
    for t = 1:T
        if do_poisson_noise % gives solutions for Poisson noise
            E_solution(:,:,t) = Poisson_quadratic_solution(kappa,gamma_poisson,B(:,:,:,t),D(:,:,:,t),Z(:,:,t));
             B(:,:,:,t) = (B(:,:,:,t)+D(:,:,:,t))./(1+kappa*gamma_poisson-gamma_poisson*Z(:,:,t)./(E_solution(:,:,t)+eps));
        else   % gives solutions for Gaussian noise
            E_solution(:,:,t) = analysis_cardano_solutions(sigma,gamma_gauss,B(:,:,:,t),D(:,:,:,t),Z(:,:,t));
             B(:,:,:,t) = (B(:,:,:,t)+D(:,:,:,t))./(2*gamma_gauss/sigma^2*E_solution(:,:,t)+1);
        end
    end
    %% ---------------- (4) Update Lagrange variables------------------------
    if do_Lagrange && s>=it_start1
        D = D-1*(B-Bold);
    end
    %% ---------------- (5) Backward propagation ----------------------------
    xst = zeros(size(B),'single');
    for t = 1:T
        xst(:,:,:,t) = At(B(:,:,:,t)-D(:,:,:,t), Masks_set(:,:,:,t));
    end
    xs = mean(xst,4);
    %% ---------------- (6) Update of U_{o,k} by CCF filtering --------------
    if do_CCF_filtering && s>=it_start2
        [xs_ccf,~] = CCF(xs);
         xs  = (1-betta_CCF)*xs +betta_CCF*xs_ccf; % smooth injection of filtering results 
    end
    
    %% ------------- Relative errors estimation ----------------------------
    s_lambda = 0;
    for lambda = lambda_set
        s_lambda = s_lambda+1;
        Relerrs(s_lambda,s)  =  norm(x(object_coord_1,object_coord_2,s_lambda)...
            - exp(-1i*angle(trace(x(object_coord_2,object_coord_2,s_lambda)'*xs(object_coord_1,object_coord_2,s_lambda))))...
            * xs(object_coord_1,object_coord_2,s_lambda), 'fro')/norm(x(object_coord_1,object_coord_2,s_lambda),'fro'); % rel. error
    end
    %% ------------ Draw figures ------------------------------------------
    if (rem(s,10) == 0)||s == 2
        figure(fig1);
            semilogy(squeeze(Relerrs(:,:))','lineWidth',1), grid on  %semilogy(1:numel(Relerrs),Relerrs)
            xlabel('Iteration'), ylabel('Relative error '),% legend(num2str(lambda_set'*1e9),'NumColumns',2),... 'Orientation','horizontal'
            title(['Relative errors vs. iteration count, T=',num2str(T),', K=',num2str(K) '. ' noise_type ', SNR=',num2str(mean(SNR),3) ' dB'])
            set(gca,'FontSize',12)
            lgnd = legend( num2str(lambda_set'*1e9,3) ,'NumColumns',1,'FontSize',8,'location','bestoutside'); % ,'Orientation','horizontal'
            title(lgnd,'\lambda, nm');
        
        figure(fig2);
        subplot 223, imshow(abs(xs(Nzp/2+(1:n1),Nzp/2+(1:n2),lambda_to_show)),[]), ...
            title('Reconstructed amplitude'), ...
            c = colorbar;
            c.Label.String = 'Amplitude, a.u.';
        subplot 224, imshow(angle(xs(Nzp/2+(1:n1),Nzp/2+(1:n2),lambda_to_show)),[]),
            title('Reconstructed phase'),
            c = colorbar;
            c.Label.String = 'Phase, rad';
            sgt = sgtitle([noise_type ' SNR=' num2str(mean(SNR),3)  ' dB, \lambda=' num2str(lambda_set(lambda_to_show)*1e9,3) ' nm, {\it ERROR_{rel}}=' num2str(Relerrs(lambda_to_show,s),2), ', iter=' num2str(s)]);
            sgt.FontSize = 8;
        drawnow
    end
  fprintf ('\n iteration %u, %2.1f sec., reler=%1.4f', s, toc, mean(squeeze(Relerrs(:,s))))
end

% figure(3), sliceViewer(abs(xs)), title('Reconstructed amplitude'),
% figure(4),sliceViewer(angle(xs)), title('Reconstructed phase')
fprintf ('\n End of reconstructions, mean relative error = %1.3g \n', mean(Relerrs(:,s)))

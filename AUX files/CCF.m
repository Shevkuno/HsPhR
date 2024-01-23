function [Output,kf]=CCF(Z)

%% CCF.m , Complex-domain Cube Filtering from paper "Hyperspectral phase imaging based on denoising in
% complex-valued eigensubspace",I.Shevkunov, V.Katkovnik, D.Claus, G.Pedrini, N.V.Petrov, K.Egiazarian.
% 
% I.Shevkunov and V.Katkovnik,  2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCF is a complex domain version of the algorithm presented in
% L. Zhuang and J. Bioucas-Dias, “Fast hyperspertral image denoising and inpainting", 
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote
% Sensing vol. 11, no. 3, pp. 730-742, 2018;
% with MATLAB codes in http://www.lx.it.pt/~bioucas/publications.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Inputs:
% Z-Noisy Hyperspectral Cube, Nx x Ny x Nz, Nz is a number of slices
%
% Outputs:
% Output - complex domain filtering of Z
% kf - size of used eigensubspace
%--------------------------------------------------------------------------

%% %%%%%%% Generation of Eigen Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Lines, Columns, B] = size(Z);
N=Lines*Columns;

temp = transpose(reshape(Z, N, B));   %% Line images
clear Z

%% estimation of random CD noise and correlation matrix %%%%%%%%%%%%%%%%%%%
[w, Rw] = estNoise_CD(temp,'additive');
[kf, E]=hysime_CD(temp,w,Rw); %%

eigen_im=E'*temp;   %% Eigen Images Lines
k_subspace=kf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% YY - Eigen images of Original shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YY = zeros(Lines,Columns,k_subspace,'single');

for ss = 1:k_subspace
    YY(:,:,ss) = single( reshape(eigen_im(ss,:), Lines, Columns));
end

clear eigen_im
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Complex Domain BM3D Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eigen_im_filtered = zeros(Lines,Columns,k_subspace,'single');
eigen_im_filtered = zeros(Lines,Columns,k_subspace,'single');
for ss = 1: k_subspace 
    
    sigma = function_stdEst([real(YY(:,:,ss)); imag(YY(:,:,ss))]); %% STD of the additive Gaussian noise
    
     eigen_im_filtered(:,:,ss)  = CDID(YY(:,:,ss), sigma, 'analysis', 1); %% Thresholding filtering
%     eigen_im_filtered(:,:,ss) = CDID(YY(:,:,ss), sigma, 'analysis', 3,'thresholding',2,'basicestimate',eigen_im_filtered_temp ); %% Wiener filtering
    
end
clear eigen_im_filtered_temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return to original space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = transpose(reshape(eigen_im_filtered, N, k_subspace));  %% to the spectral lines
temp=E*temp;

Output=zeros(Lines, Columns, B,'single');

for ss=1:B
    Output(:,:,ss) = reshape(temp(ss,:), Lines, Columns);  %% filtered cube Z
end


end

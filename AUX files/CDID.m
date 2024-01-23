function CDBM3D_hat = CDID(Z, sigma, varargin)
%
% COMPLEX DOMAIN IMAGE DENOISING TOOLBOX 
%
% CDBM3D_hat = CDID(Z, sigma, varargin);
%
% The function "CDID" performs higher-order singular value decomposition (HOSVD) 
% based denoising of complex-valued image Z for various settings.
% 
% Inputs
% ---------
%
% Required:
%
%   Z - noisy complex-valued 2D image
%
%   sigma - standard deviation of noise in Real and Imaginary components of Z
%
% Optional:
%
%   'block matching': 1 (default value) - in complex domain, 2 - in phase domain, 3 - in amplitude domain;
%      
%   '2d or 3d': 2 - pseudo 2D blocks, 3 - 3D (or 4D) blocks. Default value is 3;
%
%   'analysis': 1 (default value) - analysis of complex values, 2 - joint analysis of Phase and 
%               Amplitude components, 3 - Joint analysis of Real and Imaginary components, 4 - Separate 
%               analysis of Phase and Amplitude components;
%
%   'thresholding': 0 (default value) - hard thresholding, 1 - soft thresholding, 2 - Wiener filtering;
%      
%   'cd thresholding': 0 - joint thresholding, 1 - separate thresholding of Im and Re components.
%                      The option is applied only for thresholding of complex values.
%                      Default value is 0.
%      
%   'nm1': nm1 x nm1 is the block (patch) size used for the filtering (default value is 8)
%
%   'nm2': maximum number of similar blocks (maximum size of the 3rd dimension of a 3D array);
%          default value calculates as floor(nmfact*nm1*nm1/2);  
%
%   'nmfact': corfficient used for calculation of default value of nm2 (default value of nmfact is 1)
%      
%   'ns': length of the side of the search neighborhood for full-search block-matching (BM), 
%         must be odd. Default value of ns is 39.
%      
%   'nstep': sliding step to process every next reference block (default value is 3)
%      
%   'lambda3d': factor for hard and soft thresholding (default value is 1)
%      
%   'twosided': 1 - to use two-sided confidence interval for grouping of similar patches;
%               0 - to uses only upper bound;
%               Default value is 0.  
%      
%   'lambda_conf_int': parameter for the two sides symmetric confidence interval for grouping;
%                      Default value is 2;         
%      
%   'lambda_conf_int_factor': parameter for the upper bound in the two sides non-symmetric 
%                             confidence interval for grouping. Default value is 10.
%
%   'lambda_conf1': lower bound of confidence interval for grouping of similar patches.
%                   Default value is calculated using 'lambda_conf_int' parameter.
%
%   'lambda_conf2': upper bound of confidence interval for grouping of similar patches.
%                   Default value is calculated using 'lambda_conf_int' and 'lambda_conf_int_factor' 
%                   parameters.
%
%   'puma usage': 0 - PUMA script is not applied to noisy phase, 1 - PUMA script is applied to noise phase.
%                 Default values is 0.
%
%   'basicestimate': Basic estimate for Wiener filtering. Is not needed for hard and soft thresholding.
%      
%   'comments':  to show or not to show comments (default value is false)
%
% Outputs
% ---------
%   CDBM3D_hat - filtered complex-valued 2D image
%
% Examples of usage
% -----------------
%
% C = CDID(Z, 0.2);
%
% C = CDID(Z, 0.1, 'lambda3d', 0.9, 'analysis', 3);

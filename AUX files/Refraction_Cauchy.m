function [n,glass] =Refraction_Cauchy(lambda,glass_type)

%% calculate refractive indeces for lambda_0, lambda and focal distance for lambda

%% B and C depend on material
% Material                   B       C (mkm2)
% Fused silica              1.4580  0.00354
% Borosilicate glass BK7  1.5046  0.00420
% Hard crown glass K5       1.5220  0.00459
% Barium crown glass BaK4   1.5690  0.00531
% Barium flint glass BaF10  1.6700  0.00743
% Dense flint glass SF10  1.7280  0.01342

%% Cauchy formula
%% Take
switch glass_type
    case 1 %BK7
        glass='BK7';
B=  1.5046; C=0.00420*1e-12;
% n_0=B+C/lambda_0^2;
    case 2 %fused silica
        glass='fused silica';
        B=  1.4580; C=0.00354*1e-12;
    case 3 %
        glass='Hard crown glass K5';
        B=  1.5220; C=0.00459*1e-12;
    case 4 %
        glass='Barium crown glass BaK4';
        B=  1.5690; C=0.00531*1e-12;
    case 5 %
        glass='Barium flint glass BaF10';
        B=  1.6700 ; C=0.00743*1e-12;
    case 6 %
        glass='Dense flint glass SF10';
        B=  1.7280; C=0.01342*1e-12;
    case 7 %
        glass='approach to LCD';
        B=  1.5; C=0.05*1e-12;
end

%%
n=B+C./(lambda.^2);
% lambda_0 - reference wavelength
% lambda
% f=f_0*(n_0-1)/(n-1);


end


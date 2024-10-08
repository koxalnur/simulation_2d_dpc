%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_dpc -                                                              %
% Main file for DPC absorption and phase recovery                         %
%                                                                         %                                        %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

set(0, 'DefaultFigureWindowStyle', 'normal');
addpath('dpc_functions_sim');
F             = @(x) fft2(x);
IF            = @(x) ifft2(x);

%% load data
%load DPC images, subtract DC, and normalized by the total energy
 
IDPC_a1 = 1; 
IDPC_p1 = mat2gray(double(imread('cameraman.tif'))/255);
groundtruth_phase = IDPC_p1;
IDPC_p1 = imresize(IDPC_p1,[256,256])./max(max(IDPC_p1));
IDPC_i1 = IDPC_a1.*exp(1i.*IDPC_p1);
subplot(121);imshow(angle(IDPC_i1),[]);title('phase object');colorbar;caxis([0,1]);
subplot(122);imshow(abs(IDPC_i1),[]);title('amplitude object');colorbar;caxis([0.9,1.1]);

%% system parameters
dim           = [size(IDPC_i1, 1), size(IDPC_i1, 2)]; % image size
sigma         = 1.25;                                 % partial coherence factor
na            = 0.40;                                 % numerical aperture of the imaging system
na_illum      = sigma*na;                             % numerical aperture of the illumination
magnification = 20;                                   % magnification of the imaging system
lambda1       = 0.521;                                % wavelength in micron
ps            = 6.5/magnification;                    % pixel size in micron
wavenumber    = 2*pi/lambda1;                         % wave number
illu_rotation1 = [0, 180, 90, 270];                   % orientation of the DPC half circular patterns           
num_rotation1  = numel(illu_rotation1);               % number of illumination used in DPC 
na_inner      = [0, 0, 0, 0];                         % if annular illumination is used,
                                                      % set the na corresponds to the inner radius     
SystemSetupSim();

%% generate illumination sources
%CHANGE THE SOURCE PATTERNS BASED ON YOUR ILLUMINATION SETTINGS
titles1={'Bottom','Top','Left','Right'};
source1      = zeros(dim(1), dim(2), num_rotation1); 
for source_index = 1:num_rotation1
    source1(:, :, source_index) = genSourceSim(illu_rotation1(source_index), na_illum, na_inner(source_index),...
                                                                                                lambda1, Fx, Fy);
end

figure('Name', 'LED Illumination Patterns', 'NumberTitle', 'off');
fig_rows    = floor(sqrt(num_rotation1));
fig_cols    = floor((num_rotation1)/fig_rows);
for fig_col = 1:num_rotation1
    fig_index = fig_col;
    ax        = subplot(fig_rows, fig_cols, fig_index);
    imagesc(fftshift(fx), fftshift(fy), fftshift(source1(:, :, fig_col)));axis image; colorbar;axis off; 
    title([titles1(fig_col)]); colormap(ax, 'gray'); caxis([0, 0.1]);
end
drawnow;

%% Cameraman images before Tikhonov  
% show simulations
tables1   = {'Top','Bottom','Right','Left'};
Hi_new1   = zeros(dim(1), dim(2), num_rotation1);
Hr_new1   = zeros(dim(1), dim(2), num_rotation1);

% pupil function set by the numerical aperture of the imaging system
pupil1    = (Fx.^2+Fy.^2<=(na/lambda1)^2);

fIDPC_i1=F(IDPC_i1);

for source_index = 1:num_rotation1
    [Hi_new1(:,:,source_index), Hr_new1(:,:,source_index)] = genTransferFunctionSim(source1(:, :, source_index), pupil1);
end


figure('Name', 'DPC Images', 'NumberTitle', 'off');
for source_index=1:num_rotation1
    subplot(2, 2, source_index);
    IDPC1(:,:,source_index)   = imag(IF(Hi_new1(:, :, source_index).*(fIDPC_i1)));
    imagesc(fftshift(fx),fftshift(fy),(IDPC1(:,:,source_index)));axis image;colormap('gray');colorbar;axis off; 
    title([titles1(source_index)], 'FontSize', 12);
end


%% DPC amplitude and phase recovery with Tikhonov regularization
% calculate frequency spectrum of the measurements
fIDPC1           = F(IDPC1);
% parameters for Tikhonov regurlarization [absorption, phase] ((need to tune this based on SNR)
reg_tik1         = 1.0*[1e-1, 5e-3];
% DPC deconvolution to solve the absorption and phase
[absorption1,...
 phase1]         = DPC_tik_sim(fIDPC1, source1, pupil1, reg_tik1);

figure('Name', 'DPC Recovery Results', 'NumberTitle', 'off');
ax1 = subplot(1, 2, 1);
imagesc(x, y, (absorption1)); axis image;  colormap(ax1, 'gray'); colorbar; axis off; caxis([0.95, 1.15]); 
title('recovered \alpha','FontSize', 14);
ax2 = subplot(1, 2, 2);
imagesc(x, y, (phase1)); axis image; colorbar; colormap(ax2, 'gray'); axis off; caxis([-1.5, 1.5]);
title('recovered \phi', 'FontSize', 14);
linkaxes([ax1, ax2]);


%% Forward PTF models
figure('Name', 'forward PTF models', 'NumberTitle', 'off');
tables1   = {'Top','Bottom','Right','Left'};
Hi_final1 = zeros(dim(1), dim(2), num_rotation1);
Hi_new1   = zeros(dim(1), dim(2), num_rotation1);
Hr_new1   = zeros(dim(1), dim(2), num_rotation1);

for source_index = 1:num_rotation1
    [Hi_new1(:,:,source_index), Hr_new1] = genTransferFunctionSim(source1(:, :, source_index), pupil1);
end

for source_index = 1:num_rotation1
    subplot(2, 2, source_index);
    Hi_final1(:, :, source_index)  =  fftshift(imag(Hi_new1(:, :, source_index)));      
    imagesc(fftshift(fx),fftshift(fy), Hi_final1(:, :, source_index)); 
    colorbar; colormap('jet'); axis image; caxis([-1 1]); axis off;  
    title([tables1(source_index)], 'FontSize', 14); 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_ColorDPC_v1.m            
%                                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
set(0,'DefaultFigureWindowStyle','normal');
addpath('Code');

F = @(x) fftshift(fft2(ifftshift(x)));
IF = @(x) fftshift(ifft2(ifftshift(x)));

%% load data
%load DPC images, subtract DC, and normalized by the total energy
IDPC_a2 = 1; %imresize(double(imread('westconcordorthophoto.png'))/255, [256,256]);
IDPC_p2 = mat2gray(double(imread('cameraman.tif'))/255);
IDPC_p2 = imresize(IDPC_p2,[256,256])./max(max(IDPC_p2));
IDPC_i2 = IDPC_a2.*exp(1i.*IDPC_p2);
subplot(121);imshow(angle(IDPC_i2),[]); title('phase object');colorbar;%caxis([0,3]);
subplot(122);imshow(abs(IDPC_i2),[]); title('amplitude object');colorbar;caxis([0.9,1.1]);

%% parameters
sigma = 1.25;                            % partial coherence factor
NA = 0.4;                             % system numerical aperture
NA_illu = sigma*NA;                   % illumination numerical aperture
Magnification = 20;                 % total magnification
rotation2= [0, 180, 90, 270];        % rotation of patterns for r,g,b channels

lambda2 = [0.621 0.466 0.621 0.466];       % wavelengths in micron
ps = 6.5/Magnification;              % pixel size in micron
na_inner  = [0, 0, 0, 0]; 
num_rotation2  = numel(rotation2);     % number of illumination used in DPC 

deltaX = 256;
deltaY = 256;
yMin   = 1;
xMin   = 1;
yMax   = yMin+deltaY-1;
xMax   = xMin+deltaX-1;

Intensity2 = IDPC_i1(yMin:yMax,xMin:xMax,:);

num_wavelength2 = length(lambda2);  
nangle2 = length(rotation2);

rows2 = size(Intensity2,1);  % dimension of the image
cols2 = size(Intensity2,2);           

% Cartesian coordinate setting
x2 = -(cols2-mod(cols2,2))/2:1:(cols2-mod(cols2,2))/2-(mod(cols2,2)==0);
x2 = ps*x2;
y2 = -(rows2-mod(rows2,2))/2:1:(rows2-mod(rows2,2))/2-(mod(rows2,2)==0);
y2 = ps*y2;
[X2, Y2] = meshgrid(x2,y2);

dfx2 = 1/cols2/ps; dfy2 = 1/rows2/ps;
fx2 = dfx2*(-(cols2-mod(cols2,2))/2:1:(cols2-mod(cols2,2))/2-(mod(cols2,2)==0));
fy2 = dfy2*(-(rows2-mod(rows2,2))/2:1:(rows2-mod(rows2,2))/2-(mod(rows2,2)==0));
[Fx2, Fy2] = meshgrid(fx2,fy2);


%% Source and transfer functions
Source2        = zeros(rows2, cols2, num_wavelength2); 

for n = 1:num_wavelength2
    Source2(:,:,n)  = SourceComputeSim(rotation2(n), NA_illu, lambda2(n), Fx2, Fy2);
    Pupil2(:,:,n)    = ((Fx2.^2+Fy2.^2)*(0.521)^2 < NA^2);
    [Hi2(:,:,n) Hr2(:,:,n)] = TransferFunction_ColorDPC_Sim(lambda2(n), Source2(:,:,n), Pupil2(:,:,n)); 

end

%% Generate illumination sources

Source_red1    = zeros(rows2, cols2, 3);
Source_blue1   = zeros(rows2, cols2, 3);
Source_red2    = zeros(rows2, cols2, 3);
Source_blue2   = zeros(rows2, cols2, 3);

Source_red1(:,:,1)   =  Source2(:,:,1);
Source_blue1(:,:,3)  =  Source2(:,:,2);
Source_RGB1          =  Source_red1 + Source_blue1;

Source_red2(:,:,1)   =  Source2(:,:,3);
Source_blue2(:,:,3)  =  Source2(:,:,4);
Source_RGB2          =  Source_red2 + Source_blue2;

figure('Name', 'LED Illumination Patterns for TB', 'NumberTitle', 'off');
subplot(221);
imshow(Source_red1); 
title('Source red bottom'); 

subplot(222);
imshow(Source_blue1); 
title('Source blue top'); 

subplot(223);
imshow(Source_red2); 
title('Source red left'); 

subplot(224);
imshow(Source_blue2); 
title('Source blue right'); 

%% RGB Composites

figure('Name', 'RGB Composites', 'NumberTitle', 'off');
subplot(121);
imshow(Source_RGB1); 
title('Source RGB1-TB');

subplot(122);
imshow(Source_RGB2); 
title('Source RGB2-LR');

%% forward PTF models

figure('Name', 'forward PTF models', 'NumberTitle', 'off');
tables2={'Red bottom','blue top','red left','blue right'};

Hi_final2 = zeros(rows2, cols2, num_rotation2);
Hi2       = zeros(rows2, cols2, num_rotation2);
Hr2       = zeros(rows2, cols2, num_rotation2);
for n = 1:num_wavelength2
    [Hi2(:,:,n) Hr2(:,:,n)] = TransferFunction_ColorDPC_Sim(lambda2(n), Source2(:,:,n), Pupil2(:,:,n)); 
end

for n2 = 1:4
    subplot(2, 2, n2);
    Hi_final2(:, :, n2)   =  imag(Hi2(:, :, n2));      
    imagesc((fx2),(fy2),(Hi_final2(:, :, n2)));colorbar;colormap('jet');axis image;axis off;%caxis([-1 1]);   
    title([tables2(n2)], 'FontSize', 14);
end


%% DPC Images before Tikhonov-TBLR

figure('Name', 'Images before Tikhonov-TB & LR', 'NumberTitle', 'off');
titles2={'Red bottom','Blue top', 'Red left','Blue right'};

fintensity2=F(Intensity2);

for n = 1:num_wavelength2
    Source2(:,:,n)  = SourceComputeSim(rotation2(n), NA_illu, lambda2(n), Fx2, Fy2);
    Pupil2(:,:,n)    = ((Fx2.^2+Fy2.^2)*(0.521)^2 < NA^2); % pupil function set by the numerical aperture of the imaging system
    [Hi2(:,:,n) Hr2(:,:,n)] = TransferFunction_ColorDPC_Sim(lambda2(n), Source2(:,:,n), Pupil2(:,:,n)); 

end

for source_index2=1:4
    subplot(2, 2, source_index2);
    IDPC2(:,:,source_index2) = imag(IF(Hi2(:,:, source_index2).*fintensity2));
    imagesc((fx2),(fy2),IDPC2(:,:,source_index2));axis image;axis off;colormap('gray');
    title([titles2(source_index2)], 'FontSize', 24);
end
%%
IDPC2_dualTB = zeros(rows2, cols2,3);
IDPC2_dualTB(:,:,1)=IDPC2(:,:,1);
IDPC2_dualTB(:,:,3)=IDPC2(:,:,2);

figure;
subplot(121);imshow(IDPC2_dualTB);

IDPC2_dualLR = zeros(rows2, cols2,3);
IDPC2_dualLR(:,:,1)=IDPC2(:,:,3);
IDPC2_dualLR(:,:,3)=IDPC2(:,:,4);

subplot(122);imshow(IDPC2_dualLR);

%% Phase and amplitude retrieval
% calculate frequency spectrum of the measurements
fIDPC2      = F(IDPC2);
reg_amp2    = 1e-2;                % regularization for amplitude
reg_phase2  = 1e-1;                % regularization for phase

[amplitude2, phase2] = ColorDPC_L2_Sim(IDPC2, Hr2, Hi2, reg_amp2, reg_phase2, lambda2(2));

figure; subplot(121);imagesc(x2,y2,amplitude2); axis off; axis image; colorbar; colormap gray;caxis([0.95 1.15]);
title('Amplitude2','fontsize',20);
subplot(122);imagesc(x2,y2,phase2);axis off;axis image;colorbar;colormap gray;caxis([-1.5 1.5]);
title('Phase2','fontsize',20);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_ColorDPC_v1.m               
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
set(0,'DefaultFigureWindowStyle','normal');
F = @(x) fftshift(fft2(ifftshift(x)));
IF = @(x) fftshift(ifft2(ifftshift(x)));

%% load data
% load DPC images, subtract DC, and normalized by the total energy
IDPC_a3 = 1; %imresize(double(imread('westconcordorthophoto.png'))/255, [256,256]);
IDPC_p3 = mat2gray(double(imread('cameraman.tif'))/255);
IDPC_p3 = imresize(IDPC_p3,[256,256])./max(max(IDPC_p3));
IDPC_i3 = IDPC_a3.*exp(1i.*IDPC_p3);
subplot(121);imshow(angle(IDPC_i3),[]); title('phase object');colorbar;%caxis([0,3]);
subplot(122);imshow(abs(IDPC_i3),[]); title('amplitude object');colorbar;caxis([0.9,1.1]);


%% parameters
sigma = 1.25;                         % partial coherence factor
NA = 0.4;                             % system numerical aperture
NA_illu = sigma*NA;                   % illumination numerical aperture
Magnification = 20;                   % total magnification
rotation3 = [0, 120, 240];            % rotation of patterns for r,g,b channels

lambda3 = [0.621 0.521 0.466];        % wavelengths in micron
ps = 6.5/Magnification;              % pixel size in micron
na_inner      = [0, 0, 0, 0]; 
num_rotation3  = numel(rotation3);    % number of illumination used in DPC 

num_wavelength3 = length(lambda3);  
nangle3 = length(rotation3);

rows3 = size(IDPC_i3,1);  % dimension of the image
cols3 = size(IDPC_i3,2);           

% Cartesian coordinate setting
x3 = -(cols3-mod(cols3,2))/2:1:(cols3-mod(cols3,2))/2-(mod(cols3,2)==0);
x3 = ps*x3;
y3 = -(rows3-mod(rows3,2))/2:1:(rows3-mod(rows3,2))/2-(mod(rows3,2)==0);
y3 = ps*y3;
[X3, Y3] = meshgrid(x3,y3);

dfx3 = 1/cols3/ps; dfy3 = 1/rows3/ps;
fx3 = dfx3*(-(cols3-mod(cols3,2))/2:1:(cols3-mod(cols3,2))/2-(mod(cols3,2)==0));
fy3 = dfy3*(-(rows3-mod(rows3,2))/2:1:(rows3-mod(rows3,2))/2-(mod(rows3,2)==0));
[Fx3, Fy3] = meshgrid(fx3,fy3);

%% Source and transfer functions
Source3       = zeros(rows3, cols3, num_wavelength3); 
Source_red    = zeros(rows3, cols3, num_wavelength3);
Source_green  = zeros(rows3, cols3, num_wavelength3);
Source_blue   = zeros(rows3, cols3, num_wavelength3);

for n = 1:num_wavelength3
    Source3(:,:,n)  = SourceCompute120Sim(rotation3(n), NA_illu, lambda3(n), Fx3, Fy3);
    Pupil3(:,:,n)    = ((Fx3.^2+Fy3.^2)*(0.521)^2 < NA^2);
    [Hi3(:,:,n) Hr3(:,:,n)] = TransferFunction_ColorDPC_120Sim(lambda3(n), Source3(:,:,n), Pupil3(:,:,n)); 

end

%% Generate illumination sources

Source_red(:,:,1)   =  Source3(:,:,1);
Source_green(:,:,2) =  Source3(:,:,2);
Source_blue(:,:,3)  =  Source3(:,:,3);
Source_RGB          =  Source3;

figure('Name', 'LED Illumination Patterns for TB', 'NumberTitle', 'off');
subplot(131);
imshow(Source_red); 
title('Source red bottom'); 

subplot(132);
imshow(Source_green); 
title('Source green'); 

subplot(133);
imshow(Source_blue); 
title('Source blue'); 

%% RGB Composites

figure('Name', 'RGB Composites', 'NumberTitle', 'off');
imshow(Source3); 
title('Source RGB');


%% Phase transfer function

figure('Name', 'PTF models', 'NumberTitle', 'off');
subplot(131);    
imagesc(imag(Hi3(:,:,1))); colorbar; colormap('jet'); axis off;  axis image; 
title('PTF Red'); caxis([-2,2]);
subplot(132);    
imagesc(imag(Hi3(:,:,2))); colorbar; colormap('jet'); axis off;  axis image; 
title('PTF Green'); caxis([-2,2]);
subplot(133);    
imagesc(imag(Hi3(:,:,3))); colorbar; colormap('jet'); axis off;  axis image; 
title('PTF Blue'); caxis([-2,2]);

%% 3-axis intensity of phase transfer function
figure;
Hi_final3 = zeros(rows3, cols3);
for i=1:num_wavelength3
    [Hi3(:,:,i) Hr3(:,:,i)] =TransferFunction_ColorDPC_120Sim(lambda3(i), Source3(:,:,i), Pupil3(:,:,i)); 
     Hi_final3 = Hi_final3  + (abs(imag(Hi3(:,:,i)))).^2 ;
end

imagesc((Hi_final3)); colorbar; colormap('jet'); axis image; axis off; 
title('3-axis PTF', 'FontSize', 24);

%% 3-axis intensity of amplitude transfer function
figure;
Hr_final3 = zeros(rows3, cols3);
for i=1:num_wavelength3
    [Hi3(:,:,n) Hr3(:,:,n)] = TransferFunction_ColorDPC_120Sim(lambda3(n), Source3(:,:,n), Pupil3(:,:,n)); 
    Hr_final3 = Hr_final3  + (Hr3(:,:,n)) ;
end
imagesc(Hr_final3); colorbar; colormap('jet'); axis image; axis off;
title('3-axis ATF', 'FontSize', 24);

%% DPC Images before Tikhonov

figure('Name', 'Images before Tikhonov-TB', 'NumberTitle', 'off');
titles3={'Red','Green', 'Blue'};
fintensity3=F(IDPC_i3);

for source_index3=1:num_wavelength3
    subplot(1, 3, source_index3);
    IDPC3(:,:,source_index3) = imag(IF(Hi3(:,:, source_index3).*fintensity3));
    imagesc((fx3),(fy3),(IDPC3(:,:,source_index3))); axis image; axis off; colormap('gray');colorbar;
    title([titles3(source_index3)], 'FontSize', 24);
end


%%
figure;
imshow(IDPC3);axis image; axis off; colormap('jet');colorbar;
title("IDPC image of 3 colors");

%% Phase and amplitude retrieval
fIDPC3      = F(IDPC3);
reg_amp3    = 1e-2;                % regularization for amplitude
reg_phase3  = 1e-1;                % regularization for phase

[amplitude3, phase3] = ColorDPC_L2_Sim(IDPC3, Hr3, Hi3, reg_amp3, reg_phase3, lambda3(2));

figure; subplot(121);imagesc(amplitude3); axis off; axis image; colorbar; colormap gray;
title('Amplitude','fontsize',20);caxis([0.95 1.15]);
subplot(122);imagesc(phase3); axis off; axis image; colorbar; colormap gray;caxis([-1.5 1.5]);
title('Phase','fontsize',20);











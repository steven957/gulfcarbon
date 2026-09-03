% ed_calculation_daily_extrapolation_with_satellite_input.m - Program to compute PAR and
% PUR profiles from satellite-derived PAR and Kd490 data
%
% Syntax:  ed_calculation_daily_extrapolation_with_satellite_input.m
%
% Inputs:
%    1) Folder locations for ed0 surface spectral irradiance file, phyto
%    shape vector files, and mean kd files
%    2) Satellite image output from 'read_pace_kd490_PAR_aph442_netcdf_single.m' script
%
% Outputs:
%    1) ed_z (spectral downwelling irradiance as a function of depth);
%    2) PAR_z (photosynthetically available irradiance as a function of
%    depth;
%    3) PUR_z (photosynthetically utilizable irradiance as a function of
%    depth
%   
% Other m-files required: 
%   1) read_pace_kd490_PAR_aph442_netcdf.m - reads in data from satellite
%   image file
%
% MAT-files required: 
%    1) '*.mat' files for filterpad shape vectors generated using
%    plot_pad_phyto_size_vectors_seabass.m ('fmicro_aph_dat.mat','fpico_aph_dat.mat'
%    2) '*.mat' files with surface spectral irradiance shape files
%    generated using plot_ed0_mean.m ('ed0_es_shape_vector.mat')
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 7 Aug 2026

%% ------------- BEGIN CODE --------------%

%close all
% clc
% clearvars

%% Load data files
ed0path = '/run/media/slohrenz/PrimaryDrive/ocean_color/pace/Ancillary/';
padpath = '/run/media/slohrenz/PrimaryDrive/ocean_color/pace/Ancillary/';
kdpath = '/run/media/slohrenz/PrimaryDrive/ocean_color/pace/Ancillary/';

% Read in satellite image data (par, kd490, aph442, lat, lon)
disp('Processing satellite image data'); % 

% Load data files 
load([ed0path,'ed0_es_shape_vector.mat']); % 'ed0_shape_mean','es_shape_mean','ed_lambda'
% load([padpath,'fmicro_aph_dat.mat']);  % 'pad_dat_fmicro','pad_lambda_fmicro'
% load([padpath,'fpico_aph_dat.mat']);  % 'pad_dat_fpico','pad_lambda_fpico'
load([kdpath,'kd_cluster_averages.mat']); % 'plt_lambda','meankdest','meankdin','meankdmid','meankdout'

% Reduce wavelength resolution if needed to reduce memory usage 
lmbd_incr = 10;
new_lambda = 400:lmbd_incr:700;
lmbd_n = length(new_lambda);
ed0_shape_mean = interp1(plt_lambda(~isnan(ed0_shape_mean))',ed0_shape_mean(~isnan(ed0_shape_mean)),...
    new_lambda,'nearest','extrap');
meankdest =  interp1(plt_lambda,meankdest,new_lambda,"nearest","extrap");
meankdin =  interp1(plt_lambda,meankdin,new_lambda);
meankdmid =  interp1(plt_lambda,meankdmid,new_lambda);
meankdout =  interp1(plt_lambda,meankdout,new_lambda);

%% Setup input variables

% Find wavelength index for 490 nm
lmbd_490_indx = find(new_lambda==490);
% lmbd_490_indx = find(plt_lambda>492 & plt_lambda<492.5);  % For full spectral resolution

% Subsample kd and aph to new_lambda
% ensure column vectors
wv_kd = wv_kd(:);
wv_aph = wv_aph(:);
new_lambda = new_lambda(:);
K = 1; % change to >1 for k-nearest

% kd
[idx1, dists1] = knnsearch(wv_kd, new_lambda, 'K', K);
  % idx is N x K (N = numel(new_lambda)), dists same size
nearest_wv_kd = wv_kd(idx1);  % nearest values (N x K)
kd_intrp = kd(:,:,idx1);

% aph
[idx2, dists2] = knnsearch(wv_aph, new_lambda, 'K', K);
  % idx is N x K (N = numel(new_lambda)), dists same size
nearest_wv_aph = wv_aph(idx2);  % nearest values (N x K)
aph_intrp = aph(:,:,idx2);
aph_intrp(:,:,new_lambda>510) = nan;

% Get image dimensions
[imagem,imagen] = size(par);

% Reshape variable to one dimensional vectors for calculations
par_rshp = reshape(par,imagem.*imagen,1);   % par reshaped to vector (micromol quanta m-2 s-1)
kd_rshp = reshape(kd490,imagem.*imagen,1);  % kd490 reshaped to vector
aph_rshp = reshape(aph440,imagem.*imagen,1);   % aph440 reshaped to vector
lat_rshp = repmat(lat,imagem,1);   % lat replicated to match number of elements
lon_rshp = repmat(lon,imagen,1);   % lon replicated to match number of elements
chl_rshp = reshape(chl,imagem.*imagen,1); 
kd_intrp_rshp = reshape(kd_intrp,imagem.*imagen,size(kd_intrp,3));  
aph_intrp_rshp = reshape(aph_intrp,imagem.*imagen,size(kd_intrp,3));  
coeffs = zeros(size(aph_intrp_rshp,1),2)';

% Predictors (e.g., lambda)
wvfit_indx = find(new_lambda(new_lambda<=510 & new_lambda>=440));
wvl_pred = new_lambda(wvfit_indx);
X_reg = [ones(size(wvl_pred,1),1) wvl_pred];
aphfit_indx = find(~isnan(aph_intrp_rshp(:,1)));

% Perform regression
coeffs(:,aphfit_indx) = X_reg \ aph_intrp_rshp(aphfit_indx,wvfit_indx)'; % Linear fit for every pixel

% Get max aph for each pixel
% Calculate maximum aph for each pixel
% [maxAph, idx] = max(aph_intrp_rshp, [], 2);  % Maximum aph across wavelengths for each pixel
            
% lat_rshp = reshape(lat,imagem.*imagen,1);   % lat reshaped to vector
% lon_rshp = reshape(lon,imagem.*imagen,1);   % lon reshaped to vector

% Calculate euphotic depth for PP calculations
zeu = -log(0.01)./kd_rshp;   % Estimated euphotic depth from Chakraborty et al. (2017) and Lehrter et al. (2009) 

% Generate 10 depth increments based on zeu
z_array=zeros(size(zeu,1),lmbd_n,10);
z_incr=zeros(size(zeu,1),1);

for iz = 1:10
    z_incr = zeu./10;
    z_array(:,:,iz) = repmat(z_incr.*iz,1,lmbd_n);
end

% if zeu > 10
%     z_incr = 2;
% else
%     z_incr = 1;
% end
% z = 0:z_incr:round(zeu,0);  

disp('Variables reshaped. Beginning surface irradiance and kd calculations...')

%% Surface irradiance calculations
% Convert satellite PAR (micromol quanta m-2 s-1) to surface spectral irradiance using ed0 shape vector:
%    m x n array with pixels as rows and wavelengths as columns
par_molQ_perm2_perh = 3600.*par_rshp./1000000; % mol Q or mol photons m-2 h-1
ed0_molQ_perm2_perh = (par_molQ_perm2_perh.*ed0_shape_mean)./sum(ed0_shape_mean,'omitnan');  % mol Q or mol photons m-2 h-1

% Divide by average cosine to get scalar irradiance
e0_molQ_perm2_perh = ed0_molQ_perm2_perh./avg_cos;

if strcmp(sensor_type,'MODIS')
    % Convert kd490 to spectral k using appropriate shape vector
    kd_out_indx = find(kd_rshp<0.07);
    kd_mid_indx = find(kd_rshp>=0.07 & kd_rshp<0.2);
    kd_in_indx = find(kd_rshp>=0.2 & kd_rshp<1.2);
    kd_est_indx = find(kd_rshp>=1.2);

    % Preallocate spectral k variable
    kd_spectral = single(zeros(size(kd_rshp,1),length(new_lambda)));

    kd_spectral(kd_out_indx,:) = kd_rshp(kd_out_indx).*meankdout./meankdout(lmbd_490_indx);
    kd_spectral(kd_mid_indx,:) = kd_rshp(kd_mid_indx).*meankdmid./meankdmid(lmbd_490_indx);
    kd_spectral(kd_in_indx,:) = kd_rshp(kd_in_indx).*meankdin./meankdin(lmbd_490_indx);
    kd_spectral(kd_est_indx,:) = kd_rshp(kd_est_indx).*meankdest./meankdest(lmbd_490_indx);
elseif strcmp(sensor_type,'PACE')    % Alternatively use satellite-derived kd
    kd_spectral = kd_intrp_rshp;
end

disp('Completed surface irradiance and kd calculations. Loading and interpolating shape vectors...')

%% Prepare aph shape vectors
% Load aph data and calculate aph440 normalized shape vector
% 
%  fmicro shape vector
% load([padpath,'fmicro_aph_dat.mat']);
% fmicro_aph = pad_dat_fmicro;
% aph440_micro = fmicro_aph(pad_lambda_fmicro(:,1)>=439.9 & pad_lambda_fmicro(:,1)<=440.1,:);
% aph_micro_440_norm = pad_dat_fmicro./aph440_micro; % aph normalized to aph440
% fmicro_aph_shape = mean(aph_micro_440_norm,2); %The fmicro shape vector is the mean aph_micro_440_norm
%             %  of all the fmicro-dominated stations
% 
% % Interpolate aph and aph shape vector to same wavelengths as ed_z
% fmicro_aph_shape_interp = interp1(pad_lambda_fmicro(end:-1:1,1),fmicro_aph_shape(end:-1:1),plt_lambda);
% fmicro_aph_shape_interp = interp1(pad_lambda_fmicro(:,1),fmicro_aph_shape,new_lambda);
% 
% %  fpico shape vector 
% load([padpath,'fpico_aph_dat.mat']);
% fpico_aph = pad_dat_fpico;
% aph440_pico = fpico_aph(pad_lambda_fpico(:,1)>=439.9 & pad_lambda_fpico(:,1)<=440.1,:);
% aph_pico_440_norm = pad_dat_fpico./aph440_pico;
% fpico_aph_shape = mean(aph_pico_440_norm,2); %The fpico shape vector is the mean aph_micro_440_norm
%             %  of all the fpico-dominated stations
% fpico_aph_shape_interp = interp1(pad_lambda_fpico(:,1),fpico_aph_shape,new_lambda);
% 
% % Save shape vector file for later use
% save([padpath,'GC2_aph_shape_vectors.mat'],'fmicro_aph_shape_interp','fpico_aph_shape_interp');

load([padpath,'GC2_aph_shape_vectors.mat']);  % wavelength-interpolated aph
% shape vectors normalized to aph440

disp('Loading shape vectors completed. Beginning irradiance and PUR vs. depth calculations...')

%% Irradiance and PUR calculations
% Compute irradiance as a function of depth, looping through depths and 
%    calculation of PAR as a function of depth by integration of ed_z over
%    wavelength:  edz is a m x n x p x q array with pixel, wavelength,
%    depth, and hourly cos_solzen as dimensions; all others are m x p x q arrays with pixel, depth
%    and cos_solzen as dimensions

% Preallocate ed_z, PAR_z, and PUR_z variables
ed_z = single(zeros(size(kd_spectral,1),size(kd_spectral,2),size(z_array,3),size(cos_indx,1)));
PAR_z = single(zeros(size(kd_spectral,1),size(z_array,3),size(cos_indx,1)));
PUR_z_micro = single(zeros(size(kd_spectral,1),size(z_array,3),size(cos_indx,1)));
PUR_z_piconano = single(zeros(size(kd_spectral,1),size(z_array,3),size(cos_indx,1)));

% for idep = 1:size(z,2)
%     ed_z(:,:,idep) = ed0_molQ_perm2_perh.*exp(-kd_spectral.*z(:,idep)); % mol Q or mol photons m-2 h-1 
%     PAR_z(:,idep) = sum(ed_z(:,:,idep),2,'omitnan');  % units: mol Q or mol photons m-2 h-1; 
%     % Calculate PUR (integration of ed_z across wavelength spectrum
%     PUR_z_micro(:,idep) = sum(ed_z(:,:,idep).*fmicro_aph_shape_interp,2,'omitnan'); % mol Q or mol photons m-2 h-1
%     PUR_z_pico(:,idep) = sum(ed_z(:,:,idep).*fpico_aph_shape_interp,2,'omitnan'); % mol Q or mol photons m-2 h-1
% end

for iedz = 1:length(cos_indx)
    % for idep = 1:size(z,2)
    % 
    %     % Compute irradiance as a function of depth, adjusting surface
    %     %   irradiance to account for cosine of solar zenith angle at time of
    %     %   profile by multiplying times ed0
    %     ed_z(:,:,idep,iedz) = cos_solzen_ratio(iedz).*e0_molQ_perm2_perh.*exp(-kd_spectral.*z(:,idep)); % mol Q or mol photons m-2 h-1 nm-1
    % 
    %     % Sensitivity analysis
    %     % ed_z = ed0_molQ_perm2_perh(:,1).*exp(-(1.5.*kd_calc_mean(1,:)').*z); % mol Q or mol photons m-2 h-1 nm-1
    %     % ed_z = ed0_molQ_perm2_perh(:,1).*exp(-(0.5.*kd_calc_mean(1,:)').*z); % mol Q or mol photons m-2 h-1 nm-1
    % 
    % end

    % Compute irradiance as a function of depth, adjusting surface
    %   irradiance to account for cosine of solar zenith angle at time of
    %   profile by multiplying times ed0
    ed_z(:,:,:,iedz) = cos_solzen_ratio(iedz).*e0_molQ_perm2_perh.*exp(-kd_spectral.*z_array); % mol Q or mol photons m-2 h-1 nm-1

    % Calculation of PAR as a function of depth by integration of ed_z over wavelength
    % Limit spectral range to 400-700 nm (PAR)
    % lambda_400_700_indx = plt_lambda>399.8 & plt_lambda<700.2;
    PAR_z(:,:,iedz) = sum(ed_z(:,:,:,iedz),2,'omitnan');  % units: mol Q or mol photons m-2 h-1
    
    % Second method to calculate PAR using constant value of Kd
    % Estimate KPAR for upper water column
    % dep_rng = 10;
    % par_mdl=polyfit(z(:,z<=dep_rng),log(PAR_z(z<=dep_rng,iedz)),1);
    % KPAR = -par_mdl(1);
    % PAR_0 = exp(par_mdl(2));
    
    % PAR_z_constK(:,iedz) = PAR_z(1,iedz).*exp(-KPAR.*z);
    
    % PUR Calculation
    
    % Calculate PUR (integration of ed_z across wavelength spectrum with 3.34 nm scaling factor
    %   to account for finite bandwidth of HyperPro measurement)
    
    PUR_z_micro(:,:,iedz) = sum(ed_z(:,:,:,iedz).*fmicro_aph_shape_interp,2,'omitnan');
    % mol Q or mol photons m-2 h-1
    PUR_z_piconano(:,:,iedz) = sum(ed_z(:,:,:,iedz).*fpico_aph_shape_interp,2,'omitnan');
    % mol Q or mol photons m-2 h-1
end

profile_depth = z_array(:,1,:);
z = reshape(profile_depth,size(profile_depth,1),size(profile_depth,3));

disp('Completed irradiance and PUR calculations...');



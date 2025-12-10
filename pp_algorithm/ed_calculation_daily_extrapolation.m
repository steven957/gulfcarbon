% ed_calculation_daily_extrapolation.m - Program to compute PAR and PUR profile
% from HyperPro-derived ed and kd values and filterpad aph values. This
% program is called by pp_daily_algorithm.m script.
%
% Syntax:  ed_calculation_daily_extrapolation.m
%
% Inputs:
%    1) '*.mat' files for HyperPro estimate near-surface kd and ed values
%    generated using hyperpro_kd_averages.m
%    2) '*.mat' files for filterpad shape vectors generated using
%    plot_pad_phyto_size_vectors_seabass.m
%
% Outputs:
%    1) '*.mat' files with PAR and PUR profiles for primary production
%    calculations with units converted to mol Q or mol photons m-2 h-1
%   
% Other m-files required: None 
%
% MAT-files required: 
%  1) *.mat' files for HyperPro estimate near-surface kd and ed values
%    generated using hyperpro_kd_averages_all_Apr2009.m
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 9 Dec 2025

%% ------------- BEGIN CODE --------------%

% close all
% clc
% clearvars

% Comment this statement out for pp algorithm (this program is called by 
%   pp_daily_algorithm.m script)

% cruisename = 'GC5';

switch cruisename
    case 'GC2'
        propath = '\GC2_hyperPro_data\';
        padpath = '\Apr2009\Filterpad\';
    case 'GC3'
        propath = '\GC3_hypersensors\HyperPro_GC3\';
        padpath = '\Jul2009\Filterpad\';
    case 'GC4'
        propath = '\GC4_hyper\HyperPro_data_GC4\';
        padpath = '\Nov2009\Filterpad\';
    case 'GC5'
        propath = '\GC5_hyperPro\';
        padpath = '\Mar2010\Filterpad\';
end
 
% Comment these statements out for pp algorithm
% sta_txt = {'g3'};
% istn = 1;

if exist([propath,sta_txt{istn},'_hyperpro_profile.mat'],'file')
    load([propath,sta_txt{istn},'_hyperpro_profile.mat']);
else
    kd490 = nan;
    return
end

disp(['Processing station ',sta_txt{istn}]); % Create a depth variable
% 
% z = 0:1:100;
% z = 0:2:38;

% Convert irradiance units from uW cm-2 nm-1 to mol quanta m-2 s-1 nm-1
% (divide by hplanck * c_light/lambda_m, which is the energy associated with a given wavelength)

hplanck = 6.62607e-034;  %J Hz-1 (Planck's constant)
c_light = 3.0e008;  %Speed of light in m/s
lambda_m = plt_lambda./10^9;  % Convert nm to m
Av_num = 6.022e023;  % Avogadro's number

% Retrieve kd490 from HyperPro profiles 
kd490_indx = find(plt_lambda>488 & plt_lambda<491);
kd490 = kd_calc_mean(kd490_indx);

% Calculate euphotic depth for PP calculations
zeu = -log(0.01)./kd490;   % Estimated euphotic depth from Chakraborty et al. (2017) and Lehrter et al. (2009) 
if zeu > 10
    z_incr = 1;
else
    z_incr = 0.25;
end
z = 0:0.25:round(zeu,0);  

%% Prepare aph shape vectors
% Load aph data and calculate aph440 normalized shape vector
% 
%  fmicro shape vector
load([padpath,'fmicro_aph_dat.mat']);
fmicro_aph = pad_dat_fmicro;
aph440_micro = fmicro_aph(pad_lambda_fmicro(:,1)>=439.9 & pad_lambda_fmicro(:,1)<=440.1,:);
aph_micro_440_norm = pad_dat_fmicro./aph440_micro; % aph normalized to aph440
fmicro_aph_shape = mean(aph_micro_440_norm,2); %The fmicro shape vector 
%  is the mean aph_440-normalized spectrum of all the fmicro-dominated stations

% Interpolate aph and aph shape vector to same wavelengths as ed_z
fmicro_aph_shape_interp = interp1(pad_lambda_fmicro(:,1),fmicro_aph_shape,plt_lambda);

%  fpiconano shape vector 
load([padpath,'fpiconano_aph_dat.mat']);
fpiconano_aph = pad_dat_fpiconano;
aph440_piconano = fpiconano_aph(pad_lambda_fpiconano(:,1)>=439.9 & pad_lambda_fpiconano(:,1)<=440.1,:);
aph_piconano_440_norm = pad_dat_fpiconano./aph440_piconano;
fpiconano_aph_shape = mean(aph_piconano_440_norm,2); % The fpiconano shape 
%   vector is the mean aph_440-normalized spectrum of all the fpiconano-dominated stations

% Interpolate aph and aph shape vector to same wavelengths as ed_z
fpiconano_aph_shape_interp = interp1(pad_lambda_fpiconano(:,1),fpiconano_aph_shape,plt_lambda);

%% Irradiance calculations

% Convert irradiance units from uW cm-2 nm-1 to mol quanta m-2 s-1 nm-1 
%  Unit conversion:
%    1 W cm-2 = 1 J s-1 cm-2 
%    Conversion of ed_0 to units of photons per s per cm2
% 
ed0_J_pers_percm2=ed0_calc_mean'./1000000;  % converts uW to W      
% 
hc_over_lambda = hplanck.*c_light./lambda_m; % Planck's computation of energy per quanta in J s-1 cm-2
ed0_photons_pers_percm2 = ed0_J_pers_percm2./hc_over_lambda;  % Quanta or photons s-1 cm-2 nm-1

ed0_molQ_perm2_perh = ed0_photons_pers_percm2.*3600.*1.0e004/Av_num;  % mol Q or mol photons m-2 h-1 nm-1

% Divide by average cosine to get scalar irradiance
e0_molQ_perm2_perh = ed0_molQ_perm2_perh./avg_cos;

% Compute light levels for different hour increments through photoperiod

PAR_z = zeros(length(z),length(cos_indx));
PAR_z_constK = zeros(length(z),length(cos_indx));
PUR_z_micro = zeros(length(z),length(cos_indx));
PUR_z_piconano = zeros(length(z),length(cos_indx));

for iedz = 1:length(cos_indx)

    % Compute irradiance as a function of depth, adjusting surface
    %   irradiance to account for cosine of solar zenith angle at time of
    %   profile and dividing z by cosine solar zenith angle for each hour
    %   increment to account for longer effective pathlength
    ed_z = e0_molQ_perm2_perh(:,1).*exp(-kd_calc_mean(1,:)'.*z./cos_solzen(cos_indx(iedz))); % mol Q or mol photons m-2 h-1 nm-1
    
    % For sensitivity analysis (adjust kd by +/-50%)
    % ed_z = ed0_molQ_perm2_perh(:,1).*exp(-(1.5.*kd_calc_mean(1,:)').*z); % mol Q or mol photons m-2 h-1 nm-1
    % ed_z = ed0_molQ_perm2_perh(:,1).*exp(-(0.5.*kd_calc_mean(1,:)').*z); % mol Q or mol photons m-2 h-1 nm-1
    
    % Calculation of PAR as a function of depth by integration of ed_z over wavelength
    % Limit spectral range to 400-700 nm (PAR)
    lambda_400_700_indx = plt_lambda>399.8 & plt_lambda<700.2;
    PAR_z(:,iedz) = sum(ed_z(lambda_400_700_indx,:).*3.34,1,'omitnan');  % Factor of 3.34 is wavelength interval; units: mol Q or mol photons m-2 h-1
    
    % Second method to calculate PAR using constant value of Kd
    % Estimate KPAR for upper water column
    dep_rng = 10;
    par_mdl=polyfit(z(z<=dep_rng),log(PAR_z(z<=dep_rng,iedz)),1);
    KPAR = -par_mdl(1);
    PAR_0 = exp(par_mdl(2));
    
    PAR_z_constK(:,iedz) = PAR_z(1,iedz).*exp(-KPAR.*z);
    
    % PUR Calculation
    
    % Calculate PUR (integration of ed_z across wavelength spectrum with 3.34 nm scaling factor
    %   to account for finite bandwidth of HyperPro measurement)
    
    PUR_z_micro(:,iedz) = sum(ed_z(lambda_400_700_indx,:)'.*3.34.*fmicro_aph_shape_interp(lambda_400_700_indx)',2,'omitnan');
    % mol Q or mol photons m-2 h-1
    PUR_z_piconano(:,iedz) = sum(ed_z(lambda_400_700_indx,:)'.*3.34.*fpiconano_aph_shape_interp(lambda_400_700_indx)',2,'omitnan');
    % mol Q or mol photons m-2 h-1
end

% Save values for PP calculation
% 
profile_depth = z;
% save([padpath,sta_txt{istn},'_PUR_aph_dat.mat'],'profile_depth','PAR_z','KPAR','PAR_z_constK','PUR_z_micro','PUR_z_pico',...
%     'aph440_micro','aph440_pico','fmicro_aph_shape_interp','fmicro_aph_interp',...
%     'fpiconano_aph_shape_interp','fpiconano_aph_interp');

disp('Completed');



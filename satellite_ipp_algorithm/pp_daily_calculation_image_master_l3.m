% Plot pp_daily_calculation_image_master.m
% Syntax: pp_calculation_image_master
%
% This script calls another script ('ed_calculation_with_satellite_par_input_single_file.m')
% to generate spectral irradiance and PUR depth profiles for each pixel.
% It then generates PP estimates for each pixel vs. depth and integrates to
% produce a water column primary production values for each pixel.  
%
% Inputs:
%   1) Folder location with 'ed_calculation_with_satellite_par_input_single_file.m' script
%   2) Folder location with 'read_l3_kd490_PAR_aph440_netcdf_single_file' script
%    
% Outputs:
%   1) Figures with plotted results
%  
% Other m-files required: 
%   1) read_pace_kd490_PAR_aph442_netcdf_single_file.m - reads in data from satellite
%   image file
%   2) ed_calculation_with_satellite_par_input_single_file.m - calculates spectral
%   downwelling irradiance and PAR and PUR as a function of depth
% 
% MAT-files required (files required for ed_calculation script): 
%   1) 'ed0_es_shape_vector.mat' - shape vector for surface spectral irradiance
%   2) 'fmicro_aph_dat.mat' - shape vector for microphytoplankton absorption
%   3) 'fpiconano_aph_dat.mat' - shape vector for pico- and nanophytoplankton absorption
%   4) 'kd_cluster_averages.mat' - mean diffuse spectral attenuation
%          coefficients for downwelling irradiance for different water mass
%          types (see Mili and Lohrenz, 2026)
%
% Author: Steven E. Lohrenz
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 5 Sep 2026

%% ------------- BEGIN CODE --------------%% 

clc
clearvars

% Load file information
% codepath = '/home/user/matlab/'; % Path to matlab code
codepath = '/home/slohrenz/ocean_color/matlab/';

% Specify number of depth intervals for calculations
% z_incr = 2;  % Depth interval
% z = 0:z_incr:50;
% % z = 0:2:38;
depth_n = 10; %length(z);

% Get number of files to be read
%labeltext=cellstr(input_file(1:12));

sensor_type = 'PACE'; % 'MODIS'; % 'PACE'

monthrange = 1:12;

for mnthn = monthrange

    % List of input files for IOP and PAR
    switch sensor_type
        case 'PACE'
            prod_flag = {'kd','aph','chlor_par'};
            iop_list = {'kd','aph','chlor_par'};
        case 'MODIS'
            prod_flag = {'OC','aph'};
        end
    
    nprod = length(prod_flag);
    
    for iprod = 1:nprod
        prod_type = prod_flag{iprod};
        sensor_type = sensor_type; % 'MODIS'; % 'PACE'
        
        % Get list of input files for IOP and PAR and read in data
        switch sensor_type
            case 'PACE'
                inputfolder=['/home/user/pace/','L3_',prod_type,'/']; % Path to imagery data products
                file_lst=dir([inputfolder,'L3_',prod_type,'/','PACE_OCI.2025',sprintf('%02d', mnthn),'*_',prod_type,'.L3.nc']);
                iop_list = {prod_type};
                % Load data
                fname=char(file_lst.name);

                filepath = [inputfolder,'L3_',prod_type,'/',fname];
                finfo = ncinfo(filepath);

                disp(['Reading image data from ',fname]);

                for ivar = 1:length(iop_list)
                    iop = iop_list{ivar};
                    run([codepath,'read_l3_kd490_PAR_aph440_netcdf_single_file']);
                end

            case 'MODIS'                
                if strcmp(prod_type,'OC') 
                    iop_list = {'kd','chlor_par'};
                elseif strcmp(prod_type,'aph')
                    iop_list = {'aph'};
                end% Transpose to match dimensions of PUR

                inputfolder='/home/user/modis/';  % Path to imagery data products
                file_lst=dir([inputfolder,'A2025',sprintf('%02d', mnthn),'*_',prod_type,'.L3.nc']);

                % Load data
                fname=char(file_lst.name);

                filepath=[inputfolder,fname];
                finfo = ncinfo(filepath);

                disp(['Reading image data from ',fname]);

                for ivar = 1:length(iop_list)
                    iop = iop_list{ivar};
                    run([codepath,'read_l3_kd490_PAR_aph440_netcdf_single_file']);
                end
        end
    end
    
    %% Get solar elevation and azimuth for photoperiod
    
    % Define lat_center and lon_center for sun angle information -
    %  based on attributes from image
    lat_center = (lat_south + lat_north)./2;
    lon_center = (lon_east + lon_west)./2;
    
    % Time increment in hours 
    hr_incr = 2;
    mDatetm = datetime(dt,'Format','uuuu/MM/dd HH:mm:ss') + days((0:hr_incr:23)'./24);  % Increment days in hour intervals
    
    [sAz,sEl] = SolarAzEl(mDatetm,zeros(size(mDatetm,1),1)+lat_center,zeros(size(mDatetm,1),1)+lon_center,zeros(size(mDatetm,1),1));
    
    % Find indices and times of sunrise and sunset
    sol_zen1 = 90 - sEl;
    cos_solzen_air = cosd(sol_zen1);
    utc_h_from_midnght = mod((0:23.98),24);
    cos_indx = find(cos_solzen_air>0);
    
    solzen2_test = asind(sind(acosd(cos_solzen_air))./1.34); 
    
    % Apply Snell's law to get transmitted angle across air-water interface and compute cosine (cos_solzen_wtr)
    sol_zen2 = asind(sind(sol_zen1(cos_indx))./1.34);   % Divide sin of theta by water refractive index (1.34)
    cos_solzen_wtr = cosd(sol_zen2);
    
    % Optional to plot cos_solzen_air vs utc_h_from_midnght
    % clf
    % plot(utc_h_from_midnght(cos_solzen_air>0),cos_solzen_wtr,'bo')
    
    % Adjust irradiance values by cos_solzen_wtr, referencing time surface
    % irradiance measurement was acquired
    
    % Find time closest to irradiance profile
    [~,~,~,h,mi,s] = datevec(utc);
    pro_hr = h; % + mi/60 + s/3600;
    
    % Get cos_solzen_wtr for irradiance profile time
    cos_solzen_pro_tm = cos_solzen_wtr(hour(mDatetm(cos_indx))==pro_hr);
    cos_solzen_ratio = cos_solzen_wtr./cos_solzen_pro_tm;
    
    %% Calculate light profiles
    
    % Approximate value for Ed/Eo
    avg_cos = 0.8; % Assumed value for converting planar to scalar irradiance
    
    % Call script to compute irradiance profile data
    ed_calculation_daily_extrapolation_with_satellite_input_git;
    
    % Get aph slope and determine mean 
    % 
    % if log(abs(aph_slope_var))>-8.3  % Approach based on Verma et al. (2021)
    %     PUR_z=PUR_z_micro; % mol Q or mol photons m-2 h-1
    % else based
    %     PUR_z=PUR_z_piconano; % mol Q or mol photons m-2 h-1
    % end
    
    %% Productivity Calculations
    % Get indices for micro- and piconano-dominated stations (based on aph440 threshold; see Hirata et al., 2008)
    % micro_indx = find(aph_rshp>0.1);
    % piconano_indx = find(aph_rshp<=0.1);
    
    % Alternatively use aph slope from PACE or MODIS QAA aph product
    micro_indx = find(log(abs(coeffs(2,:)))>-8.3);
    piconano_indx = find(log(abs(coeffs(2,:)))<=-8.3);
    
    % Estimate aph440 normalized pmax (pmax_aph440_spec_est) and quantum
    % yield (Based on results from Mili and Lohrenz, 2026; see pmax_box_plot.m)
    pmax_aph440_spec_est = 0.0217*(1 - exp(-5.69*kd_rshp));   % mol C m-2 h-1

    % Estimated phi based on correlation with pmax_aph440_spec ((Based on results from Mili and Lohrenz, 2026; see  pmax_aph_spec_vs_phi.m)
    phi_est = -60.3.*pmax_aph440_spec_est.^2 + 5.17.*pmax_aph440_spec_est - 0.00354;
    
    avg_cos = 0.8;  % Estimated value based on literature
    
    % Estimate standard P-E par using relationships in Mili and Lohrenz(2026)
    alpha_est = zeros(size(aph_rshp));
    pbmax_est = pmax_aph440_spec_est.*aph_rshp./chl_rshp;  % mol C m-2 h-1 mg Chl-1
    alpha_est(micro_indx) = phi_est(micro_indx).*mean((aph_rshp(micro_indx)./chl_rshp(micro_indx)).*fmicro_aph_shape_interp,2);
    alpha_est(piconano_indx) = phi_est(piconano_indx).*mean((aph_rshp(piconano_indx)./chl_rshp(piconano_indx)).*fpiconano_aph_shape_interp,2);
    
    disp('Beginning productivity calculations...');
    
    % Run pp algorithm on profile for each hour, adjusting irradiance
    % for cosine of solar zenith angle (Based on Mili and Lohrenz, 2026)
    
    clear P_z_est_day P_z_PAR_day
    
    % Pre-allocate P_est
    P_z_est = zeros(size(aph_rshp,1),10);
    P_z_std = zeros(size(aph_rshp,1),10);
    
    for icos = 1:length(cos_indx)
    
        % Calculate PP for micro-dominated stations (mol C m-3 h-1)
        P_z_est(micro_indx,:) = aph_rshp(micro_indx).*pmax_aph440_spec_est(micro_indx).*...
            (1 - exp(-phi_est(micro_indx).*(PUR_z_micro(micro_indx,:,icos)./avg_cos)./(pmax_aph440_spec_est(micro_indx)))); % mol C m-3 h-1
        
        % Calculate PP for piconano-dominated stations (mol C m-3 h-1)
        P_z_est(piconano_indx,:) = aph_rshp(piconano_indx).*pmax_aph440_spec_est(piconano_indx).*...
            (1 - exp(-phi_est(piconano_indx).*(PUR_z_piconano(piconano_indx,:,icos)./avg_cos)./(pmax_aph440_spec_est(piconano_indx)))); % mol C m-3 h-1
    
        % Standard wavelength-integrated, PAR-based algorithm for comparison
        % P_z_std = chl_rshp.*pbmax_est.*(1 - exp(-alpha_est.*(PAR_z(:,:,icos)./avg_cos)./pbmax_est))./12000; % mol C m-3 h-1
    
        if exist('P_z_est_day','var')
            P_z_est_day = P_z_est_day + P_z_est.*hr_incr;
            % P_z_std_day = P_z_std_day + P_z_std;
        else
            P_z_est_day = P_z_est.*hr_incr;
            % P_z_std_day = P_z_std;
        end
    
    end
    
    % Calculate water column integrated primary production (mol C m-2 h-1)
    PP_est_int = sum(P_z_est_day.*z_incr,2);
    PP_est_int_array = reshape(PP_est_int,imagem,imagen);
    lat_rshp = reshape(lat_rshp,imagen,imagem)';
    lon_rshp = reshape(lon_rshp,imagem,imagen);
        
    disp('Completed productivity calculations. Saving data...');
    
    % Save output
    % Prepare filename for saving
    file_name = strrep(strrep(fname,'_aph',''),'.nc','');
    file_name = [strrep(strrep(file_name,'_chlor_par',''),'.nc',''),'.PP.mat'];
    save([inputfolder,file_name],'P_z_est','z','PP_est_int','PP_est_int_array','lat_rshp','lon_rshp','lat','lon');
    
    disp(['Saving ',file_name]);
end

disp(' ');
disp('Completed');

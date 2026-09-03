% Plot pp_daily_calculation_image_master.m
% Syntax: pp_calculation_image_master
%
% This script calls another script ('ed_calculation_with_satellite_par_input_single_file.m')
% to generate spectral irradiance and PUR depth profiles for each pixel.
% It then generates PP estimates for each pixel vs. depth and integrates to
% produce a water column primary production values for each pixel. 
%
% Inputs:
%   1) Folder location with 'ed_calculation_with_satellite_par_input_single_file.m'
%   2) Folder location with 'read_pace_kd490_PAR_aph442_netcdf_single_file.m' script
%           or 'read_pace_l3_kd490_PAR_aph442_netcdf_single_file.m' script
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
% MAT-files required: None
%
% Author: Steven E. Lohrenz
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 13 Aug 2026

%% ------------- BEGIN CODE --------------%% 

clc
clearvars

% Load file information
codepath = '/home/slohrenz/ocean_color/matlab/';

% Specify number of depth intervals for calculations
% z_incr = 2;  % Depth interval
% z = 0:z_incr:50;
% % z = 0:2:38;
depth_n = 10; %length(z);

% Get number of files to be read

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
        
        % List of input files for IOP and PAR
        switch sensor_type
            case 'PACE'
                inputfolder=['/run/media/slohrenz/PrimaryDrive/ocean_color/pace/','L3_',prod_type,'/'];
                file_lst=dir([inputfolder,'PACE_OCI.2025',sprintf('%02d', mnthn),'*_',prod_type,'.L3.nc']);
                iop_list = {prod_type};
                % Load data
                fname=char(file_lst.name);

                %labeltext=cellstr(input_file(1:12));
                filepath=[inputfolder,fname];
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
                end
                inputfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/modis/';
                file_lst=dir([inputfolder,'A2025',sprintf('%02d', mnthn),'*_',prod_type,'.L3.nc']);

                % Load data
                fname=char(file_lst.name);

                %labeltext=cellstr(input_file(1:12));
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
    % sun_indx = find(cos_sEl > 0);
    
    solzen2_test = asind(sind(acosd(cos_solzen_air))./1.34); 
    
    % Apply Snell's law to get transmitted angle across air-water
    % interface
    sol_zen2 = asind(sind(sol_zen1(cos_indx))./1.34);   % Divide sin of theta by water refractive index (1.34)
    cos_solzen_wtr = cosd(sol_zen2);
    
    % day_lngth = utc_h_from_midnght(sun_indx(1))-utc_h_from_midnght(sun_indx(2));
    % sunset = utc_h_from_midnght(sun_indx(1));
    % sunrise = utc_h_from_midnght(sun_indx(2));
    
    % Optional to plot cos_sEl vs utc_h_from_midnght
    % clf
    % plot(utc_h_from_midnght(cos_solzen_air>0),cos_solzen_wtr,'bo')
    
    % Adjust irradiance values by cos_sE1, referencing time surface
    % irradiance measurement was acquired
    
    % Find time closest to irradiance profile
    [~,~,~,h,mi,s] = datevec(utc);
    pro_hr = h; % + mi/60 + s/3600;
    % pro_tm = floor(pro_hr);
    % if isnan(pro_tm)
    %     continue;
    % end
    
    % Get cos_sEl for irradiance profile time
    cos_solzen_pro_tm = cos_solzen_wtr(hour(mDatetm(cos_indx))==pro_hr);
    cos_solzen_ratio = cos_solzen_wtr./cos_solzen_pro_tm;
    
    %% Calculate light profiles
    
    % Approximate value for Ed/Eo
    avg_cos = 0.8;
    
    % Load light profile data
    ed_calculation_daily_extrapolation_with_satellite_input;
    
    % % Get aph slope?
    % 
    % if log(abs(aph_slope_var))>-8.3
    %     PUR_z=PUR_z_micro; % mol Q or mol photons m-2 h-1
    %     aph440_mean = mean(aph440_micro);
    % else 
    %     PUR_z=PUR_z_pico; % mol Q or mol photons m-2 h-1
    %     aph440_mean = mean(aph440_pico);
    % end
    
    %% Productivity Calculations
    % Get indices for micro- and pico-dominated stations (based on aph440
    % threshold)
    % micro_indx = find(aph_rshp>0.1);
    % pico_indx = find(aph_rshp<=0.1);
    
    % Alternatively use aph slope from PACE QAA aph product
    micro_indx = find(log(abs(coeffs(2,:)))>-8.3);
    pico_indx = find(log(abs(coeffs(2,:)))<=-8.3);
    
    % Estimate aph440 normalized pmax (pmax_aph440_spec_est) and quantum yield
    % (phi) (Based on results from pmax_box_plot.m)
    pmax_aph440_spec_est = 0.0217*(1 - exp(-5.69*kd_rshp));   % mol C m-2 h-1

    % Estimated phi based on correlation with pmax_aph440_spec (see pmax_aph_spec_vs_phi.m)
    phi_est = -60.3.*pmax_aph440_spec_est.^2 + 5.17.*pmax_aph440_spec_est - 0.00354;
    % (Based on correlation between pmax_aph440_spec and phi in
    %           'PEcurves_GC-phyto_allvar_revised_051624.xlsx')
    
    avg_cos = 0.8;  % Estimated value based on literature
    
    % Estimate standard P-E parameters
    alpha_est = zeros(size(aph_rshp));
    pbmax_est = pmax_aph440_spec_est.*aph_rshp./chl_rshp;  % mol C m-2 h-1 mg Chl-1
    alpha_est(micro_indx) = phi_est(micro_indx).*mean((aph_rshp(micro_indx)./chl_rshp(micro_indx)).*fmicro_aph_shape_interp,2);
    alpha_est(pico_indx) = phi_est(pico_indx).*mean((aph_rshp(pico_indx)./chl_rshp(pico_indx)).*fpico_aph_shape_interp,2);
    
    % P_z = aph_440.*pmax_aph440_spec.*(1 - exp(-aph440_mean.*phi.*(PUR_z./avg_cos)./(aph_440.*pmax_aph440_spec))); % mol C m-3 h-1
    % P_z_PAR = chl.*pbmax.*(1 - exp(-alpha.*(PAR_z_constK./avg_cos)'./pbmax))./12000; % mol C m-3 h-1
    
    disp('Beginning productivity calculations...');
    
    % Run pp algorithm on profile for each hour, adjusting irradiance
    % for cosine of solar zenith angle
    
    clear P_z_est_day P_z_PAR_day
    
    % Pre-allocate P_est
    P_z_est = zeros(size(aph_rshp,1),10);
    P_z_std = zeros(size(aph_rshp,1),10);
    
    for icos = 1:length(cos_indx)
    
        % Calculate PP for micro-dominated stations (mol C m-3 h-1)
        P_z_est(micro_indx,:) = aph_rshp(micro_indx).*pmax_aph440_spec_est(micro_indx).*...
            (1 - exp(-phi_est(micro_indx).*(PUR_z_micro(micro_indx,:,icos)./avg_cos)./(pmax_aph440_spec_est(micro_indx)))); % mol C m-3 h-1
        
        % Calculate PP for pico-dominated stations (mol C m-3 h-1)
        P_z_est(pico_indx,:) = aph_rshp(pico_indx).*pmax_aph440_spec_est(pico_indx).*...
            (1 - exp(-phi_est(pico_indx).*(PUR_z_piconano(pico_indx,:,icos)./avg_cos)./(pmax_aph440_spec_est(pico_indx)))); % mol C m-3 h-1
    
        P_z_std = chl_rshp.*pbmax_est.*(1 - exp(-alpha_est.*(PAR_z(:,:,icos)./avg_cos)./pbmax_est))./12000; % mol C m-3 h-1
    
        if exist('P_z_est_day','var')
            P_z_est_day = P_z_est_day + P_z_est.*hr_incr;
            % P_z_PAR_day = P_z_PAR_day + P_z_PAR;
        else
            P_z_est_day = P_z_est.*hr_incr;
            % P_z_PAR_day = P_z_PAR;
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

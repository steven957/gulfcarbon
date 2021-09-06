% hypersas_data_read.m - Script to read hypersas data from individual run files for further
%  processing; called 
%
% Syntax: hypersas_data_read
%
% Inputs:
%    1) 'file_name' variable from main program (see below)
%
% Outputs:
%    Variables used in processing HyperSAS reflectance observations
%   
% Other m-files required: 
%  1) 'plot_hypersas_hyperpro_rrs_Apr2009_wqaa_output_v8c.m'
%      or 'hypersas_rrs_output_to_seabass_v4.m', which provide 'file_name'
%      variable
%
% MAT-files required: none
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 6 Sep 2021
%
%% ------------- BEGIN CODE --------------%

disp(text_msg);
load(file_name);
disp('Loaded.')

%Populate variables
rnum=rad_rrs_output.run_no==run_no;
if iscell(rad_rrs_output.es_data)  
    es_data=rad_rrs_output.es_data{rnum};
    es_data_interp=rad_rrs_output.es_data_interp{rnum};
    gps_time=rad_rrs_output.gps_time{rnum};
    gps_datenum=rad_rrs_output.gps_datenum{rnum};
    gps_lat=rad_rrs_output.gps_lat{rnum};
    gps_lon=rad_rrs_output.gps_lon{rnum};
    gps_course=rad_rrs_output.gps_course{rnum};
    gps_speed=rad_rrs_output.gps_speed{rnum};
    es_datenum=rad_rrs_output.es_datenum{rnum};
    es_timer=rad_rrs_output.es_timer{rnum};
    es_time=rad_rrs_output.es_time{rnum};
    es_lambda=rad_rrs_output.lambda_es;
    lt_data=rad_rrs_output.lt_data{rnum};
    lw_data=rad_rrs_output.lt_data_ref_corr{rnum};
    lt_int_time=rad_rrs_output.lt_int_time{rnum};
    alt_lt_data=rad_rrs_output.alt_lt_data{rnum};
    lt_timer=rad_rrs_output.lt_timer{rnum};
    lt_time=rad_rrs_output.lt_time{rnum};
    lt_datenum=rad_rrs_output.lt_datenum{rnum};
    lt_lambda=rad_rrs_output.lambda_lt;
    li_int_time=rad_rrs_output.li_int_time{rnum};
    li_dat_interp=rad_rrs_output.li_dat_interp{rnum};
    li_data=rad_rrs_output.li_data{rnum};
    li_timer=rad_rrs_output.li_timer{rnum};
    li_time=rad_rrs_output.li_time{rnum};
    rrs=rad_rrs_output.rrs{rnum};
    alt_rrs=rad_rrs_output.alt_rrs{rnum};
    drrs=rad_rrs_output.drrs{rnum};
    alt_drrs=rad_rrs_output.alt_drrs{rnum};
    des_datenum_interp=rad_rrs_output.des_datenum_interp{rnum};
    des_time=rad_rrs_output.des_time{rnum};
    des_timer=rad_rrs_output.des_timer{rnum};
    des_data_interp=rad_rrs_output.des_data_interp{rnum};
    drrs_lambda=rad_rrs_output.lambda_dlt;
    dlt_data=rad_rrs_output.dlt_data{rnum};
    dli_data=rad_rrs_output.dli_data{rnum};
    dlw_data=rad_rrs_output.dlt_data_ref_corr{rnum};
    alt_dlt_data=rad_rrs_output.alt_dlt_data{rnum};
    dlt_timer=rad_rrs_output.dlt_timer{rnum};
    dlt_time=rad_rrs_output.dlt_time{rnum};
    rrs_lat=rad_rrs_output.rrs_lat{rnum};
    rrs_lon=rad_rrs_output.rrs_lon{rnum};
    ths_timer=rad_rrs_output.ths_timer{rnum};
    ths_roll=rad_rrs_output.ths_roll{rnum};
    ths_pitch=rad_rrs_output.ths_pitch{rnum};
    solazi=rad_rrs_output.solazi{rnum};  
    solelev=rad_rrs_output.solelev{rnum};  
    azimuth_viewing_angle=rad_rrs_output.azimuth_viewing_angle{rnum};
    wind_speed=rad_rrs_output.wind_speed{rnum};
    nir_rho=rad_rrs_output.nir_rho{rnum};
    lnir=rad_rrs_output.lnir{rnum};
    rho=rad_rrs_output.rho{rnum};
    cloud_cover=rad_rrs_output.ship_cloud_cover(rnum);
    wave_h=rad_rrs_output.ship_waveh(rnum);
    ship_ws=rad_rrs_output.ship_ws(rnum);
else
    es_data=rad_rrs_output.es_data;   %For instances with only one run per day
    es_data_interp=rad_rrs_output.es_data_interp;
    gps_time=rad_rrs_output.gps_time;
    gps_datenum=rad_rrs_output.gps_datenum;
    gps_lat=rad_rrs_output.gps_lat;
    gps_lon=rad_rrs_output.gps_lon;
    gps_course=rad_rrs_output.gps_course;
    gps_speed=rad_rrs_output.gps_speed;
    es_datenum=rad_rrs_output.es_datenum;
    es_timer=rad_rrs_output.es_timer;
    es_time=rad_rrs_output.es_time;
    es_lambda=rad_rrs_output.lambda_es;
    lt_data=rad_rrs_output.lt_data;
    lw_data=rad_rrs_output.lt_data_ref_corr;
    alt_lt_data=rad_rrs_output.alt_lt_data;
    lt_timer=rad_rrs_output.lt_timer;
    lt_time=rad_rrs_output.lt_time;
    lt_int_time=rad_rrs_output.lt_int_time;
    lt_datenum=rad_rrs_output.lt_datenum;
    lt_lambda=rad_rrs_output.lambda_lt;
    li_int_time=rad_rrs_output.li_int_time;
    li_dat_interp=rad_rrs_output.li_dat_interp;
    li_data=rad_rrs_output.li_data;
    li_timer=rad_rrs_output.li_timer;
    li_time=rad_rrs_output.li_time;
    rrs=rad_rrs_output.rrs;
    alt_rrs=rad_rrs_output.alt_rrs;
    drrs=rad_rrs_output.drrs;
    alt_drrs=rad_rrs_output.alt_drrs;
    des_datenum_interp=rad_rrs_output.des_datenum_interp;
    des_time=rad_rrs_output.des_time;
    des_timer=rad_rrs_output.des_timer;
    des_data_interp=rad_rrs_output.des_data_interp;
    drrs_lambda=rad_rrs_output.lambda_dlt;
    dlt_data=rad_rrs_output.dlt_data;
    dli_data=rad_rrs_output.dli_data;
    dlw_data=rad_rrs_output.dlt_data_ref_corr;
    alt_dlt_data=rad_rrs_output.alt_dlt_data;
    dlt_timer=rad_rrs_output.dlt_timer;
    dlt_time=rad_rrs_output.dlt_time;
    rrs_lat=rad_rrs_output.rrs_lat;
    rrs_lon=rad_rrs_output.rrs_lon;
    ths_timer=rad_rrs_output.ths_timer;
    ths_roll=rad_rrs_output.ths_roll;
    ths_pitch=rad_rrs_output.ths_pitch;
    solazi=rad_rrs_output.solazi;  
    solelev=rad_rrs_output.solelev;
    azimuth_viewing_angle=rad_rrs_output.azimuth_viewing_angle;
    wind_speed=rad_rrs_output.wind_speed;
    nir_rho=rad_rrs_output.nir_rho;
    lnir=rad_rrs_output.lnir;
    rho=rad_rrs_output.rho;
    cloud_cover=rad_rrs_output.ship_cloud_cover;
    wave_h=rad_rrs_output.ship_wave_height;
    ship_ws=rad_rrs_output.ship_ws;
end

% ------------- END OF CODE --------------
% Please send suggestions for improvement
%  to Steven Lohrenz at this email address: slohrenz@umassd.edu   

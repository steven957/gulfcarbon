% hypersas_data_filter.m - Matlab script to apply hypersas data filter and 
%   baseline corrections to HyperSAS data
%
% Syntax: hypersas_data_filter
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

%Assign indices for baseline correction 
nir_index=find(lt_lambda>775 & lt_lambda<785);
index_780=find(lt_lambda>778 & lt_lambda<782);
index_720=find(lt_lambda>719 & lt_lambda<721);
nir_index_single=nir_index(end);  %Near infrared wavelength index
nir_index_lange=find(lt_lambda>750 & lt_lambda<800);

vis_index=find(lt_lambda>400 & lt_lambda<700);
deriv_index=find(lt_lambda>610 & lt_lambda<660);
blue_index=find(lt_lambda>553 & lt_lambda<556);
uv_range_index=find(lt_lambda>349 & lt_lambda<451);
uv_index=find(lt_lambda>379 & lt_lambda<381);  %uv wavelength index
index_555=find(lt_lambda>555 & lt_lambda<557);
index_600=find(lt_lambda>598 & lt_lambda<601);

%% Interpolate pitch and roll and solar zenith values to match
% lt data using either lt_timer (pitch/roll) or lt_time (solar elev) 
pitch_intrp=interp1(ths_timer,ths_pitch,lt_timer);
roll_intrp=interp1(ths_timer,ths_roll,lt_timer);
zenith_intrp=interp1(datenum(gps_time),90-solelev,datenum(lt_time));
az_view_intrp=interp1(datenum(gps_time),azimuth_viewing_angle,datenum(lt_time));
solazi_intrp=interp1(datenum(gps_time),solazi,datenum(lt_time));
date_intrp=interp1(datenum(gps_time),gps_datenum,datenum(lt_time));
gps_course_intrp=interp1(datenum(gps_time),gps_course,datenum(lt_time));
%Convert from knots to m/s
gps_speed=gps_speed./1.944;
gps_speed_intrp=interp1(datenum(gps_time),gps_speed,datenum(lt_time));
wind_speed_intrp=interp1(datenum(gps_time),wind_speed,datenum(lt_time));

%Similarly for the discrete sensors
pitch_intrp2=interp1(ths_timer,ths_pitch,dlt_timer);
roll_intrp2=interp1(ths_timer,ths_roll,dlt_timer);
%des_time=datestr(des_datenum_interp,'HH:MM:SS');  %These steps are needed to remove inconsistencies in numeric date
zenith_intrp2=interp1(datenum(gps_time),90-solelev,datenum(dlt_time));
az_view_intrp2=interp1(datenum(gps_time),azimuth_viewing_angle,datenum(dlt_time));

%% Find indices for times within specified range and omit pitch and
% roll values >=5 degrees and relative solar azimuth angles between 50 and 170 degrees
% and solar zenith angles >10 and <80 degrees (Lange et al., 2020)
%Also eliminate extreme rrs outliers

good_lt_index=find(datenum(lt_time)>datenum(time_str1) & datenum(lt_time)<datenum(time_str2) &...
    abs(pitch_intrp)<5 & abs(roll_intrp)<5 & zenith_intrp>10 & zenith_intrp<80 &...
    az_view_intrp>10 & az_view_intrp<180 & min(rrs,[],2)>-.1 & max(rrs,[],2)<0.08);

good_dlt_index=find(datenum(dlt_time)>datenum(time_str1) & datenum(dlt_time)<datenum(time_str2) &...
    abs(pitch_intrp2)<5 & abs(roll_intrp2)<5 & zenith_intrp2>10 & zenith_intrp2<80 &...
    az_view_intrp2>10 & az_view_intrp2<180 & min(drrs,[],2)>-.1 & max(drrs,[],2)<0.08);

%Move on to next station if no acceptable data
if isempty(good_lt_index)==1 || isempty(good_dlt_index)==1
    return
end

%% Omit outliers based on point to point variations
%[rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'median','ThresholdFactor',1);  %This is a median filter
% excluding values differing by one scaled median absolute deviations (MAD)

%This can be modified to do a moving median window of specified size as below:
t=datetime(2009,4,30,0,0,0)+hours(lt_timer(good_lt_index)./3600);
%Different window length and threshold depending on whether this is for a station/sample comparison or a continuous run
if strncmp(sta_txt,'Run',3)==1
    [rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'movmedian',minutes(5),'SamplePoints',t,'ThresholdFactor',2);  %5 min moving window, threshold 2 median abs dev
else
    [rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'movmedian',minutes(20),'SamplePoints',t,'ThresholdFactor',1);  %20 min moving window
end

%This is a median filter based on the number of observations rather than time
%median_window=100;
%[rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'movmedian',median_window,'ThresholdFactor',1);  %This is a median filter

rrs_new=rrs(good_lt_index(~irrs),:);

%% Popluate output variables with consistent indexing
rrs_filt=rrs_new;
rrs_lat_filt=rrs_lat(good_lt_index(~irrs),:);
rrs_lon_filt=rrs_lon(good_lt_index(~irrs),:);
lt_data_filt=lt_data(good_lt_index(~irrs),:);
lt_time_filt=lt_time(good_lt_index(~irrs));
li_data_filt=li_dat_interp(good_lt_index(~irrs),:);
lw_data_filt=lw_data(good_lt_index(~irrs),:);
es_data_filt=es_data_interp(good_lt_index(~irrs),:);

%% Interpolate ancillary variables to match lt sampling frequency
pitch_intrp_filt=interp1(ths_timer,ths_pitch,lt_timer(good_lt_index(~irrs)));
roll_intrp_filt=interp1(ths_timer,ths_roll,lt_timer(good_lt_index(~irrs)));
zenith_intrp_filt=interp1(datenum(gps_time),90-solelev,datenum(lt_time(good_lt_index(~irrs))));
az_view_intrp_filt=interp1(datenum(gps_time),azimuth_viewing_angle,datenum(lt_time(good_lt_index(~irrs))));
solazi_intrp_filt=interp1(datenum(gps_time),solazi,datenum(lt_time(good_lt_index(~irrs))));
date_intrp_filt=interp1(datenum(gps_time),gps_datenum,datenum(lt_time(good_lt_index(~irrs))));
gps_lat_intrp_filt=interp1(datenum(gps_time),gps_lat,datenum(lt_time(good_lt_index(~irrs))));
gps_lon_intrp_filt=interp1(datenum(gps_time),gps_lon,datenum(lt_time(good_lt_index(~irrs))));
gps_course_intrp_filt=interp1(datenum(gps_time),gps_course,datenum(lt_time(good_lt_index(~irrs))));
gps_speed_intrp_filt=interp1(datenum(gps_time),gps_speed,datenum(lt_time(good_lt_index(~irrs))));
wind_speed_intrp_filt=interp1(datenum(gps_time),wind_speed,datenum(lt_time(good_lt_index(~irrs))));
rho_intrp_filt=interp1(datenum(gps_time),rho,datenum(lt_time(good_lt_index(~irrs))));

%% Omit discrete rrs outliers based on point to point variations
%[drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,4),'median','ThresholdFactor',1);  %This is a median filter
% excluding values differing by one scaled median absolute deviations (MAD)
%This can be modified to do a moving median window of specified size as below:
td=datetime(2009,4,30,0,0,0)+hours(dlt_timer(good_dlt_index)./3600);
%Different window length and threshold depending on whether this is for a station/sample comparison or a continuous run
if strncmp(sta_txt,'Run',3)==1
    [drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,4),'movmedian',minutes(5),'SamplePoints',td,'ThresholdFactor',2);  %5 min moving window, threshold 2 median abs dev
else
    [drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,4),'movmedian',minutes(20),'SamplePoints',td,'ThresholdFactor',1);  %20 min moving window
end

% by time set to a 20 min moving window
%median_window=100;
%[drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,index_555),'movmedian',median_window,'ThresholdFactor',1);  %This is a median filter

%% Popluate output variables with consistent indexing
dlt_data_filt=dlt_data(good_dlt_index(~idrrs),:);
dli_data_filt=dli_data(good_dlt_index(~idrrs),:);
dlw_data_filt=dlw_data(good_dlt_index(~idrrs),:);
des_data_filt=des_data_interp(good_dlt_index(~idrrs),:);
drrs_new=drrs(good_dlt_index(~idrrs),:);

%% Interpolate variables to match lt sampling frequency
drrs_intrp_filt=interp1(dlt_timer(good_dlt_index(~idrrs),:),drrs_new,lt_timer(good_lt_index(~irrs)));
dlt_intrp_filt=interp1(dlt_timer(good_dlt_index(~idrrs),:),dlt_data_filt,lt_timer(good_lt_index(~irrs)));
dli_intrp_filt=interp1(dlt_timer(good_dlt_index(~idrrs),:),dli_data_filt,lt_timer(good_lt_index(~irrs)));
dlw_intrp_filt=interp1(dlt_timer(good_dlt_index(~idrrs),:),dlw_data_filt,lt_timer(good_lt_index(~irrs)));
des_intrp_filt=interp1(dlt_timer(good_dlt_index(~idrrs),:),des_data_filt,lt_timer(good_lt_index(~irrs)));

%% Baseline correction
if base_corr==1
    %Baseline correct
    if mean(rrs_filt(:,uv_index))<mean(rrs_filt(:,nir_index))
        %Do exponential fit to UV portion of spectrum with
        % asymptotic approach to zero
        uv_rrs=rrs_filt(:,uv_range_index);
        uv_l=lt_lambda(uv_range_index);
        rrs_baseline=zeros(size(rrs_filt,1),1);
        for iuv=1:size(uv_rrs,1)
            g = fittype('a+b*exp(c*x)');
            f0 = fit(uv_l,uv_rrs(iuv,:)',g,'StartPoint',[0;0.0001;0.01]);
            rrs_baseline(iuv)=f0(340);
        end
        rrs_corr=rrs_filt-repmat(rrs_baseline,1,size(rrs,2)); %This adjusts spectra so Rrs(780) is zero
        drrs_corr=drrs_intrp_filt-repmat(mean(rrs_baseline),size(drrs_intrp_filt)); %This adjusts the discrete Rrs by the mean of the Rrs(780) baseline for the hyper spectra
    else
        rrs_baseline=mean(rrs_filt(:,nir_index),2);
        rrs_corr=rrs_filt-repmat(rrs_baseline,1,size(rrs,2));
        drrs_adjust=mean(drrs_intrp_filt(:,4),'omitnan')-mean(rrs_corr(:,uv_index),1,'omitnan');  %Determine offset bewteen mean drrs at 380 and mean rrs at 380 nm
        drrs_corr=drrs_intrp_filt-repmat(drrs_adjust,1,size(drrs_intrp_filt,2)); %Adjust drrs spectra by offset
    end
elseif base_corr==2 %Ruddick et al., 2005,2006
    rrs_baseline=2.35.*rrs(:,index_780)-rrs(:,index_720)/(2.35-1);
    rrs_corr=rrs-repmat(rrs_baseline,1,size(rrs,2));
else
    rrs_corr=rrs_filt;
    drrs_corr=drrs_intrp_filt;            
end

%% Further filtering to remove negative spectra
include_indx2=rrs_corr(:,index_600)>0 & min(rrs_corr,[],2)>-0.01 & min(drrs_corr,[],2)>-0.01;    

if size(find(include_indx2),1)==137
    %Need to omit at least one scan (only needed for Run7)
    [min_rrs,min_indx]=min(min(rrs_corr,[],2));
    include_indx2(min_indx)=0;
end

rrs_corr=rrs_corr(include_indx2,:);
rrs_lat_filt=rrs_lat_filt(include_indx2,:);
rrs_lon_filt=rrs_lon_filt(include_indx2,:);
lt_filt=lt_data_filt(include_indx2,:);
lt_time_filt=lt_time_filt(include_indx2,:);
li_filt=li_data_filt(include_indx2,:);
lw_filt=lw_data_filt(include_indx2,:);
es_filt=es_data_filt(include_indx2,:);
pitch_intrp_filt=pitch_intrp_filt(include_indx2,:);
roll_intrp_filt=roll_intrp_filt(include_indx2,:);
zenith_intrp_filt=zenith_intrp_filt(include_indx2,:);
az_view_intrp_filt=az_view_intrp_filt(include_indx2,:);
solazi_intrp_filt=solazi_intrp_filt(include_indx2,:);
date_intrp_filt=date_intrp_filt(include_indx2,:);
gps_lat_intrp_filt=gps_lat_intrp_filt(include_indx2,:);
gps_lon_intrp_filt=gps_lon_intrp_filt(include_indx2,:);
gps_course_intrp_filt=gps_course_intrp_filt(include_indx2,:);
gps_speed_intrp_filt=gps_speed_intrp_filt(include_indx2,:);
wind_speed_intrp_filt=wind_speed_intrp_filt(include_indx2,:);
rho_intrp_filt=rho_intrp_filt(include_indx2,:);
drrs_corr=drrs_corr(include_indx2,:);
dlt_filt=dlt_intrp_filt(include_indx2,:);
dli_filt=dli_intrp_filt(include_indx2,:);
dlw_filt=dlw_intrp_filt(include_indx2,:);
des_filt=des_intrp_filt(include_indx2,:);

%% Replace NaN with -9999 (SeaBASS default for missing data)
drrs_corr(isnan(drrs_corr))=-9999;
dlt_filt(isnan(dlt_filt))=-9999;
dli_filt(isnan(dli_filt))=-9999;
dlw_filt(isnan(dlw_filt))=-9999;
des_filt(isnan(des_filt))=-9999;

% ------------- END OF CODE --------------
% Please send suggestions for improvement
%  to Steven Lohrenz at this email address: slohrenz@umassd.edu   

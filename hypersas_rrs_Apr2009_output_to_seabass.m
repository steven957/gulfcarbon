%Process radiometer data and output to SeaBASS format
% Uses files generated with
% 'generate_hypersas_reflectance_radiance_files_all_cruise_rev6.m'

clear
clc
clf reset
close all

inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\';

input_meta_file={'hypersas_seabass_header_template_file.csv'};
%Load metadata file
disp(['Reading ',char(input_meta_file)]);

fid1=fopen([inpath,char(input_meta_file)]);
    metafile_txt=textscan(fid1,...
        '%s','delimiter','\n','TreatasEmpty','EOF');
fclose(fid1);

cruise_name='gulfcarbon2';

%fields='/fields: date,time,lat,lon,heading,speed_f_w,SZA,SAZ,RelAz,wind,pitch,roll,Lt,Lsky,Lw,Es,Rrs';
%units='/units:yyyymmdd,hh:mm:ss,degrees,degrees,degrees,degrees,degrees,degrees,degrees,degrees,degrees,m/s,uW/cm^2/nm/sr,uW/cm^2/nm/sr,uW/cm^2/nm/sr,uW/cm^2/nm,1/sr';
%NOTE: Each of the last five fields has 137 wavelengths

%Spreadsheet with input parameters
infile_name=[inpath,'apr2009_hypersas_plot_qaa_parameters.xlsx'];

single_sta='Run52';  %Single station for processing just one

%Input parameters
[~,~,all_par_3]=xlsread(infile_name,1,'C3:R28');
hyper_sta=all_par_3(:,1);
hyper_doy=cell2mat(all_par_3(:,2));
hyper_date=cell2mat(all_par_3(:,16));
hyper_run=cell2mat(all_par_3(:,3));
hyper_start_time=cell2mat(all_par_3(:,4));
hyper_end_time=cell2mat(all_par_3(:,5));
sta_lab=all_par_3(:,12);
title_lab=all_par_3(:,13);
sas_basecorr_flag=cell2mat(all_par_3(:,14));

%hyper_date=datestr(datenum(hyper_date),'yyyymmdd');

station_names={'A1','A1b','A2','A3','B4','C1','C5','E0','E1','E3','D5','F4','F5','MR2','MR3',...
'G2','G3','G4','G5','H1','H2','H3','UWS1','UWS2','Run52'};

[nsta,~]=size(hyper_sta);

%hfig=figure(1);

for im=1:nsta 

    sta_txt=hyper_sta{im};  %Select station to plot

    %To do all stations, comment out the continue statement
    if strcmp(single_sta,sta_txt)~=1
        continue
    end

    %HyperSAS parameters
    file_no=hyper_doy(im);
    run_no=hyper_run(im);
    file_date=hyper_date(im,:);
    time_str1=datestr(hyper_start_time(im),'HH:MM:SS');
    time_str2=datestr(hyper_end_time(im),'HH:MM:SS');

    station_txt=sta_lab{im};
    title_txt=title_lab{im};

    %Base correction flags
    base_corr=sas_basecorr_flag(im);    %HyperSAS base correction (1=on; 0=off)

    %file_name=[inpath,'rad_rrs_',num2str(file_no),'_',num2str(run_no),'.mat'];
    file_name=[inpath,'rad_rrs_Es253_',num2str(file_no),'.mat'];
    text_msg=['loading ',file_name,'...'];

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

    %Interpolate pitch and roll and solar zenith values to match
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

    %Find indices for times within specified range and omit pitch and
    % roll values >=5 degrees and relative solar azimuth angles between 50 and 170 degrees
    % and solar zenith angles >10 and <80 degrees (Lange et al., 2020)

    good_lt_index=find(datenum(lt_time)>datenum(time_str1) & datenum(lt_time)<datenum(time_str2) &...
        abs(pitch_intrp)<5 & abs(roll_intrp)<5 & zenith_intrp>10 & zenith_intrp<80 &...
        az_view_intrp>10 & az_view_intrp<180);

    good_dlt_index=find(datenum(dlt_time)>datenum(time_str1) & datenum(dlt_time)<datenum(time_str2) &...
        abs(pitch_intrp2)<5 & abs(roll_intrp2)<5 & zenith_intrp2>10 & zenith_intrp2<80 &...
        az_view_intrp2>10 & az_view_intrp2<180);

    %Move on to next station if no acceptable data
    if isempty(good_lt_index)==1 || isempty(good_dlt_index)==1
        continue
    end

    %Omit outliers based on point to point variations
    %[rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'median','ThresholdFactor',1);  %This is a median filter
    % excluding values differing by one scaled median absolute deviations (MAD)
    %This can be modified to do a moving median window of specified size as below:
    t=datetime(2009,4,30,0,0,0)+hours(lt_timer(good_lt_index)./3600);
    [rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'movmedian',minutes(20),'SamplePoints',t,'ThresholdFactor',1);  %This is a median filter
    % by time set to a 20 min moving window
    %median_window=100;
    %[rrs_in,irrs]=rmoutliers(rrs(good_lt_index,index_555),'movmedian',median_window,'ThresholdFactor',1);  %This is a median filter
    rrs_new=rrs(good_lt_index(~irrs),:);

    %Popluate output variables with consistent indexing
    rrs_filt=rrs_new;
    rrs_lat_filt=rrs_lat(good_lt_index(~irrs),:);
    rrs_lon_filt=rrs_lon(good_lt_index(~irrs),:);
    lt_data_filt=lt_data(good_lt_index(~irrs),:);
    lt_time_filt=lt_time(good_lt_index(~irrs));
    li_data_filt=li_dat_interp(good_lt_index(~irrs),:);
    lw_data_filt=lw_data(good_lt_index(~irrs),:);
    es_data_filt=es_data(good_lt_index(~irrs),:);
    pitch_intrp_filt=interp1(ths_timer,ths_pitch,lt_timer(good_lt_index(~irrs)));
    roll_intrp_filt=interp1(ths_timer,ths_roll,lt_timer(good_lt_index(~irrs)));
    zenith_intrp_filt=interp1(datenum(gps_time),90-solelev,datenum(lt_time(good_lt_index(~irrs))));
    az_view_intrp_filt=interp1(datenum(gps_time),azimuth_viewing_angle,datenum(lt_time(good_lt_index(~irrs))));
    solazi_intrp_filt=interp1(datenum(gps_time),solazi,datenum(lt_time(good_lt_index(~irrs))));
    date_intrp_filt=interp1(datenum(gps_time),gps_datenum,datenum(lt_time(good_lt_index(~irrs))));
    gps_course_intrp_filt=interp1(datenum(gps_time),gps_course,datenum(lt_time(good_lt_index(~irrs))));
    gps_speed_intrp_filt=interp1(datenum(gps_time),gps_speed,datenum(lt_time(good_lt_index(~irrs))));
    wind_speed_intrp_filt=interp1(datenum(gps_time),wind_speed,datenum(lt_time(good_lt_index(~irrs))));
    rho_intrp_filt=interp1(datenum(gps_time),rho,datenum(lt_time(good_lt_index(~irrs))));

    %Omit discrete rrs outliers based on point to point variations
    %[drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,4),'median','ThresholdFactor',1);  %This is a median filter
    % excluding values differing by one scaled median absolute deviations (MAD)
    %This can be modified to do a moving median window of specified size as below:
    td=datetime(2009,4,30,0,0,0)+hours(dlt_timer(good_dlt_index)./3600);
    [drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,4),'movmedian',minutes(20),'SamplePoints',td,'ThresholdFactor',1);  %This is a median filter
    % by time set to a 20 min moving window
    %median_window=100;
    %[drrs_in,idrrs]=rmoutliers(drrs(good_dlt_index,index_555),'movmedian',median_window,'ThresholdFactor',1);  %This is a median filter
    drrs_new=drrs(good_dlt_index(~idrrs),:);

    drrs_filt=drrs_new;

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
            drrs_corr=drrs_filt-repmat(mean(rrs_baseline),size(drrs_filt)); %This adjusts the discrete Rrs by the mean of the Rrs(780) baseline for the hyper spectra
        else
            rrs_baseline=mean(rrs_filt(:,nir_index),2);
            rrs_corr=rrs_filt-repmat(rrs_baseline,1,size(rrs,2));
            drrs_adjust=nanmean(drrs_filt(:,4))-nanmean(rrs_corr(:,uv_index),1);  %Determine offset bewteen mean drrs at 380 and mean rrs at 380 nm
            drrs_corr=drrs_filt-repmat(drrs_adjust,1,size(drrs_filt,2)); %Adjust drrs spectra by offset
        end
    elseif base_corr==2 %Ruddick et al., 2005,2006
        rrs_baseline=2.35.*rrs(:,index_780)-rrs(:,index_720)/(2.35-1);
        rrs_corr=rrs-repmat(rrs_baseline,1,size(rrs,2));
    else
        rrs_corr=rrs_filt;
        drrs_corr=drrs_filt;            
    end

    %Further filtering to remove negative spectra
    include_indx2=rrs_corr(:,index_600)>0;

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
    gps_course_intrp_filt=gps_course_intrp_filt(include_indx2,:);
    gps_speed_intrp_filt=gps_speed_intrp_filt(include_indx2,:);
    wind_speed_intrp_filt=wind_speed_intrp_filt(include_indx2,:);
    rho_intrp_filt=rho_intrp_filt(include_indx2,:);

    %Output data to seabass format 

    %Create field headings for radiometric wavelengths
    wavel_str=num2str(lt_lambda,'%4.1f');
    Lt_fields=[repmat(',Lt',137,1),wavel_str];
    Lsky_fields=[repmat(',Lsky',137,1),wavel_str];
    Lw_fields=[repmat(',Lw',137,1),wavel_str];
    Es_fields=[repmat(',Es',137,1),wavel_str];
    Rrs_fields=[repmat(',Rrs',137,1),wavel_str];
    
    fields=['/fields=date,time,lat,lon,heading,speed_f_w,SZA,SAZ,RelAz,wind,pitch,roll',...
       reshape(Lt_fields', 1, []),reshape(Lsky_fields', 1, []),reshape(Lw_fields', 1, []),...
       reshape(Es_fields', 1, []),reshape(Rrs_fields', 1, [])];
    units=['/units=yyyymmdd,hh:mm:ss,degrees,degrees,degrees,m/s,degrees,degrees,degrees,m/s,degrees,degrees'...
        repmat(',uW/cm^2/nm/sr',1,137.*3),repmat(',uW/cm^2/nm',1,137),repmat(',1/sr',1,137)];

    %rad_data_file_name=GulfCarbon_gulfcarbon2_hypersas_St#_2009MMDD_R1.sb
    rad_data_file_name=['GulfCarbon_gulfcarbon2_hypersas_St',sta_txt,'_',file_date,'_R1.sb'];
    
    textstring=metafile_txt;
    [hdr_m]=size(char(textstring{1}));
    file_creation_date=datestr(datenum(date),'yyyymmdd');

    %Input metadata into file header
    for i=1:hdr_m
        if ~contains(char(textstring{1}(i)),'/north_latitude')~=1
            textstring{1}(i)={['/north_latitude=',num2str(max(rrs_lat_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/south_latitude'))~=1 %#ok<*STREMP>
            textstring{1}(i)={['/south_latitude=',num2str(min(rrs_lat_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/west_longitude'))~=1
            textstring{1}(i)={['/west_longitude=',num2str(min(rrs_lon_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/east_longitude'))~=1
            textstring{1}(i)={['/east_longitude=',num2str(max(rrs_lat_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/start_time'))~=1
            textstring{1}(i)={['/start_time=',time_str1,'[GMT]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/end_time'))~=1
            textstring{1}(i)={['/end_time=',time_str2,'[GMT]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/start_date'))~=1
            textstring{1}(i)={['/start_date=',file_date]};
        elseif isempty(strfind(char(textstring{1}(i)),'/end_date'))~=1
            textstring{1}(i)={['/end_date=',file_date]};
        elseif isempty(strfind(char(textstring{1}(i)),'/cruise'))~=1
            textstring{1}(i)={['/cruise=gulfcarbon2']};
        elseif isempty(strfind(char(textstring{1}(i)),'/station'))~=1
            textstring{1}(i)={['/station=',sta_txt]};
        elseif isempty(strfind(char(textstring{1}(i)),'/cloud_percent'))~=1
            textstring{1}(i)={['/cloud_percent=',num2str(cloud_cover)]};
        elseif isempty(strfind(char(textstring{1}(i)),'/waveht'))~=1
            textstring{1}(i)={['/wave_height=',num2str(wave_h,'%3.1f')]};
        elseif isempty(strfind(char(textstring{1}(i)),'/wind_speed'))~=1
            textstring{1}(i)={['/wind_speed=',num2str(ship_ws,'%3.1f')]};
        elseif isempty(strfind(char(textstring{1}(i)),'/water_depth'))~=1
            textstring{1}(i)={['/water_depth=na']};
        elseif isempty(strfind(char(textstring{1}(i)),'/measurement_depth'))~=1
            textstring{1}(i)={['/measurement_depth=na']};
        elseif isempty(strfind(char(textstring{1}(i)),'/fields'))~=1
            textstring{1}(i)={fields};
        elseif isempty(strfind(char(textstring{1}(i)),'/units'))~=1
            textstring{1}(i)={units};
        elseif isempty(strfind(char(textstring{1}(i)),'/original_file_name='))~=1
            textstring{1}(i)={['/original_file_name=','all_hyper_Apr2009_',num2str(hyper_doy(im)),'.mat']};
        elseif isempty(strfind(char(textstring{1}(i)),'/data_file_name='))~=1
            textstring{1}(i)={['/data_file_name=',rad_data_file_name]};
        elseif isempty(strfind(char(textstring{1}(i)),'!file_creation_date='))~=1
            textstring{1}(i)={['!file_creation_date=',file_creation_date]};
        else
        end
    end
    
    %fields=['/fields: date,time,lat,lon,heading,speed_f_w,SZA,SAZ,RelAz,wind,pitch,roll',...
    %   cellstr(Lt_fields)',char(cellstr(Lsky_fields)'),char(cellstr(Lw_fields)'),char(cellstr(Es_fields)'),char(cellstr(Rrs_fields)')];
    %units=['/units:yyyymmdd,hh:mm:ss,degrees,degrees,degrees,degrees,degrees,degrees,degrees,degrees,degrees,m/s',...
    %    repmat(',uW/cm^2/nm/sr',1,137.*3),repmat(',uW/cm^2/nm',1,137),repmat(',1/sr',1,13)];

    anc_data_array=[rrs_lat_filt,rrs_lon_filt,gps_course_intrp_filt,gps_speed_intrp_filt,zenith_intrp_filt,...
        solazi_intrp_filt,az_view_intrp_filt,wind_speed_intrp_filt,pitch_intrp_filt,roll_intrp_filt];
    rad_data_array=[lt_filt,li_filt,lw_filt,es_filt,rrs_corr];
    
    %Drop decimal seconds from time
    lt_time_char=char(lt_time_filt);
    lt_time_char=lt_time_char(:,1:end-4);
    
    txt_data_array=[char(datestr(date_intrp_filt,'yyyymmdd')),repmat(',',size(date_intrp_filt,1),1),lt_time_char];
    
    [rad_data_m,rad_data_n]=size(rad_data_array);
    [anc_data_m,anc_data_n]=size(anc_data_array);
    format_str=['%s,%6.4f,%6.4f,%4.1f,%4.2f,%3.1f,%4.1f,%4.1f,%3.1f,%4.2f,%4.2f,',repmat('%.2E,',1,137.*4),repmat('%.2E,',1,136),'%.2E\n'];

    disp(['Writing ',rad_data_file_name]);
    fid2 = fopen([inpath,'To_SeaBASS\',rad_data_file_name], 'w');
        for ih=1:hdr_m
            fprintf(fid2,'%s\n',char(textstring{1}(ih)));
        end
        for id=1:rad_data_m
            fprintf(fid2,format_str,txt_data_array(id,:),anc_data_array(id,:),rad_data_array(id,:));
        end
        
    fclose(fid2);

end

fclose('all');

disp('Completed')



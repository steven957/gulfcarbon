    %Program to read hypersas data from sensor dat files and plot

    clc
    clf all reset
    clear all
    
    fig=figure(1);
    scrsz = get(groot,'ScreenSize');
    set(fig,'Position',[scrsz(4).*.1 scrsz(3).*.1 scrsz(3).*.55 scrsz(4).*.65])
    set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired
    
    %Sensor ID information
    hypersensor_list={'HED188C','HED253B','HLD175C','HLD176C','HSE188C','HSE253B','HSL175C','HSL176C',};
    discretesensor_list={'DI7130D','DR7068e','DR7069d'};
    pitchroll_sensor='SATTHS0002';
    gps_sensor='GPRMC_NoMode';

    %Load data 
    homepath='C:\Users\slohrenz\Documents\Steve\DATA\';
    cruises={'January','April','July','November','March'};
    
    
    %Lookup table for reflectance correction (Mobley, 2015)
    rho_filename=[homepath,'NASA\HyperSAS GEO-CAPE\Matlab\rhoTable_AO2015.txt'];
    max_rho=0.05;  %Maximum allowable rho (to avoid overcorrection)
    
    discrete_process=1;  %0 if no discrete sensors
    alt_rrs_plot=0;  %0 to plot Rrs using Mobley(2015) reflectance correction, 1 for Lange et al. (2020) NIR reflectance correction
    
    %Cruise loop
    for cruisen=2 %1:5

        cruise_name=cruises{cruisen};
        switch cruise_name
            case 'January'
                inpath=[homepath,'NSF\GulfCarbon\Jan2009\HyperSAS\HyperSAS_data_GC1\Output\'];
                
                %Toggle to correct for local vs UTC if needed (True for Jul 2009)
                utccorr='False';  %'True';    
               
                %Met data file
                met_file_name=[inpath,'buoy42040_2009_continuous_winds.txt'];  %Jan-Jul
                
                mat_list=dir([inpath,'all_hyper_Jan2009_*.mat']);

                initial_doy=12;
                end_doy=15;
                
            case 'April'
                inpath=[homepath,'NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\'];
                
                %Toggle to correct for local vs UTC if needed (True for Jul 2009)
                utccorr='False';  %'True';    
               
                %Met data file
                met_file_name=[inpath,'buoy42040_2009_continuous_winds.txt'];  %Jan-Jul
                
                mat_list=dir([inpath,'all_hyper_Apr2009_*.mat']);

                initial_doy=110;
                end_doy=120;

                %Load sea condition information
                inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\';
                [sea_cond,raw1,raw2]=xlsread([inpath,'apr2009_hypersas_file_list.xls'],'K3:O67');
                sea_cond_run_no=sea_cond(:,1);
                sea_cond_ws=sea_cond(:,3);    %wind speed (m/s)
                sea_cond_cloud=sea_cond(:,4);  %cloud cover (%)
                sea_cond_waveh=sea_cond(:,5);  %wave height (m)

            case 'July'
                inpath=[homepath,'NSF\GulfCarbon\GC_hypersensor_data\GC3_hypersensors\HyperSAS_GC3\Output\'];
                
                %Toggle to correct for local vs UTC if needed (True for Jul 2009)
                utccorr='True';  %'True';    
               
                %Met data file
                met_file_name=[inpath,'buoy42040_2009_continuous_winds.txt'];  %Jan-Jul
                
                mat_list=dir([inpath,'all_hyper_Jul2009_*.mat']);

                initial_doy=201;
                end_doy=210;
                
            case 'November'
                inpath=[homepath,'NSF\GulfCarbon\GC_hypersensor_data\GC4_hyper\HyperSAS_Data_GC4\Output\'];
                
                %Toggle to correct for local vs UTC if needed (True for Jul 2009)
                utccorr='False';  %'True';    
               
                %Met data file
                met_file_name=[inpath,'burl1h2009.txt'];  %Jan-Jul
                
                mat_list=dir([inpath,'all_hyper_Nov2009_*.mat']);

                initial_doy=302;
                end_doy=310;

            case 'March'
                inpath=[homepath,'NSF\GulfCarbon\GC_hypersensor_data\GC5_hyperSAS\Output\'];
                
                %Toggle to correct for local vs UTC if needed (True for Jul 2009)
                utccorr='False';  %'True';    
               
                %Met data file
                met_file_name=[inpath,'burl1h2010.txt'];  %Jan-Jul
                
                mat_list=dir([inpath,'all_hyper_Mar2010_*.mat']);

                initial_doy=70;
                end_doy=79;
        end

        [filen,filem]=size(mat_list);

        %File loop for each day in cruise
        for file_no=1:filen    %% 
            file_name=[inpath,mat_list(file_no).name];
            %Specify date information
            switch cruise_name
                case 'March'
                    outfile_label=file_name(end-5:end-4);
                    curr_date=datenum('1-1-2010')+initial_doy-1;
                    yearnum=2010;
                case 'January'
                    outfile_label=file_name(end-5:end-4);
                    curr_date=datenum('1-1-2009')+initial_doy-1;
                    yearnum=2009;
                otherwise
                    outfile_label=file_name(end-6:end-4);
                    curr_date=datenum('1-1-2009')+initial_doy-1;
                    yearnum=2009;
            end
            
            %File load step
            text_msg=['loading ',file_name,'...'];
            disp(text_msg);
            load(file_name);
            disp('Loaded.')

            %Determine number of separate runs in current file
            run_range=all_hyper.hyper_run;

            [~,m]=size(all_hyper.hyper_data);

            %Test for empty data sets and omit
            good_range=zeros(m,1);
            for ik=1:m
                good_range(ik)=isempty(all_hyper.hyper_data{ik});
            end
            
            %Select run number
            min_run=min(run_range(good_range==0));
            max_run=max(run_range(good_range==0));

            %Begin processing runs (run loop)
            clear rad_rrs_output;
            for lrun=min_run:max_run
                run_no=lrun;
                
                %Ship-based sea state observations
                ship_ws=sea_cond_ws(sea_cond_run_no==run_no);
                ship_cloud_cover=sea_cond_cloud(sea_cond_run_no==run_no);
                ship_waveh=sea_cond_waveh(sea_cond_run_no==run_no);

                text_msg2=['Run no=',num2str(run_no)];
                disp(text_msg2);
                run_index=find(all_hyper.hyper_run==run_no); %Finds sensor data corresponding to run no

                [mindx,~]=size(run_index);
                
                %This tests to see if any sensors are missing data
                for fnum=1:mindx
                    if isempty(all_hyper.hyper_data{run_index(fnum)})
                        error(''); %Error capture - comment out if desired and program will skip to next run
                        disp('Exiting run due to missing or bad data...');
                        continue
                    end
                end

                %Time range
                time_str1='08:00:00.000';
                time_str2='22:30:00.000'; 
                
                %This next section gathers the sensor indices and wavelengths for the selected run

                %Get irradiance sensor (Es) index and wavelengths
                %  (NOTE: two possible irradiance sensors, 188 or 253)
                %*************************************************
                %hyper_index=strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(5))==1; %HSE188(5)
                hyper_index=strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(6))==1; %HSE253(6)
                es_run_index=run_index(hyper_index);
                disp('Getting Es indices...');
                es_index=find(strncmp(all_hyper.hyper_data{es_run_index}.hyper_hdr,'ES',2)==1);

                %Get wavelengths for irradiance sensor
                lam_hdr=all_hyper.hyper_data{es_run_index}.hyper_hdr(es_index);
                [m, ~]=size(lam_hdr);

                lambda_es=zeros(m,1);
                for j=1:m 
                    str=char(lam_hdr(j));
                    lambda_es(j)=str2double(str(4:9));
                end

                %Select wavelength indices for plotting Es
                wave_index=find(lambda_es>380 & lambda_es<383);
                refl_index=find(lambda_es>549 & lambda_es<551);

                %Get radiance sensor (Lt) index and wavelengths
                hyper_index=strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(7))==1;  %HSL175(7)
                lt_run_index=run_index(hyper_index);
                disp('Getting Lt indices...');
                lt_index=find(strncmp(all_hyper.hyper_data{lt_run_index}.hyper_hdr,'LT',2)==1);

                %Get radiance wavelengths
                lam_hdr=all_hyper.hyper_data{lt_run_index}.hyper_hdr(lt_index);
                [m,~]=size(lam_hdr);

                lambda_lt=zeros(m,1);
                for j=1:m
                    str=char(lam_hdr(j));
                    lambda_lt(j)=str2double(str(4:9));
                end
                
                %Select near infrared Lt wavelengths for filtering glint
                %nir_index=find(lambda_lt>669 & lambda_lt<672);
                nir_index=find(lambda_lt>730 & lambda_lt<740);
                %nir_index=find(lambda_lt>750 & lambda_lt<800);
                alt_nir_index=find(lambda_lt>750 & lambda_lt<800); %Alternative nir range
                
                %Get irradiance dark correction index (two possible sensors)
                %hyper_index = strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(1))==1;  %HED188(1)
                hyper_index = strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(2))==1;  %HED253B(1)
                es_dark_run_index=run_index(hyper_index);

                %Get radiance dark correction index
                hyper_index = strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(3))==1;  %HLD175(3)
                lt_dark_run_index=run_index(hyper_index);

                %Get sky radiance sensor (Li) index and wavelengths
                hyper_index = strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(8))==1;  %HSL176(8)
                li_run_index=run_index(hyper_index);
                disp('Getting Li indices...');
                li_index=strncmp(all_hyper.hyper_data{li_run_index}.hyper_hdr,'LI',2)==1;

                %Get sky radiance wavelengths
                lam_hdr=all_hyper.hyper_data{li_run_index}.hyper_hdr(li_index);

                [m, n]=size(lam_hdr);
                lambda_li=zeros(m,1);
                for j=1:m
                    str=char(lam_hdr(j));
                    lambda_li(j)=str2double(str(4:9));
                end

                %Get sky radiance dark correction index
                hyper_index=find(strcmp(all_hyper.hyper_sensor(run_index),hypersensor_list(4))==1);  %HLD176(4)
                li_dark_run_index=run_index(hyper_index);

                %***Read irradiance data****
                es_timer=all_hyper.hyper_data{es_run_index}.hyper_timer; %Timer data show duration of run
                
                %Test for no data or bad data
                if isempty(es_timer)==1
                    disp('Advancing to next run due to no data...');
                    error(''); %Error capture - comment out if desired and program will skip to next run
                    continue
                end

                es_dat=all_hyper.hyper_data{es_run_index}.hyper_rad_dat;
                es_time=all_hyper.hyper_data{es_run_index}.hyper_time; %Clock time during run
                es_int=all_hyper.hyper_data{es_run_index}.hyper_int_time; %Sampling integration time (s)
                es_samp_delay=all_hyper.hyper_data{es_run_index}.hyper_samp_delay; %Sample delay (s)
                                
                %Remove bad entries
                es_time_bad_index=find(contains(es_time,'Error')==1);
                es_time(es_time_bad_index)={'00:00:00.000'};
                formatOut='HH:MM:SS.FFF';
                if strcmp(utccorr,'True')==1   %Convert to UTC (July only)
                    es_datenum=datenum(es_time)+5./24; %Add 5 h to local time for UTC
                else
                    es_datenum=datenum(es_time);
                end
                %es_time=datestr(es_datenum,formatOut); %Convert datenum to time string
                es_test_datenum=datenum(es_time); %Generate test value for time of day only

                %Dark correction data
                es_dark_dat=all_hyper.hyper_data{es_dark_run_index}.hyper_rad_dat;
                es_dark_timer=all_hyper.hyper_data{es_dark_run_index}.hyper_timer;
                es_dark_int=all_hyper.hyper_data{es_dark_run_index}.hyper_int_time; %Sampling integration time (s)
                es_dark_samp_delay=all_hyper.hyper_data{es_dark_run_index}.hyper_samp_delay; %Sample delay (s)

                %Eliminate duplicate values
                [C,ia,ib]=unique(es_dark_timer);

                %Dark correct irradiance data
                %Check to see if there are at least two points for interpolation
                if size(es_dark_timer(ia),1)>1 && size(es_timer,1)==size(es_dat,1)
                    es_dark_interp=interp1(es_dark_timer(ia),es_dark_dat(ia,:),es_timer,'linear','extrap'); %Interpolation to match times
                else
                    error(''); %Error capture - comment out if desired and program will skip to next run
                   continue
                end
                es_dat_dcorr=es_dat-es_dark_interp;
                es_750=es_dat_dcorr(:,lambda_es>748 & lambda_es<751);

                %Find time indices within the time range specified
                time_index=find(es_test_datenum>=datenum(time_str1) & es_test_datenum<=datenum(time_str2));

                min_index=min(time_index);
                max_index=max(time_index);

                %Check that the max index es_timer value for the beginning and
                % ending times is actually greater than the min_index value. 
                % If not, set max es_timer value to the maximum in the run.
                if es_timer(max_index)<es_timer(min_index)
                    [esmx,esmaxind]=max(es_timer);
                    max_index=esmaxind;
                end
                
                %Read data for radiance
                lt_dat=all_hyper.hyper_data{lt_run_index}.hyper_rad_dat;
                lt_timer=all_hyper.hyper_data{lt_run_index}.hyper_timer;
                lt_time=all_hyper.hyper_data{lt_run_index}.hyper_time;
                %Eliminate 'Error' entries in lt_time
                lt_time(strncmp('Error',lt_time,5)==1)={'00:00:00'};
                lt_int=all_hyper.hyper_data{lt_run_index}.hyper_int_time; %Sampling integration time (s)
                lt_samp_delay=all_hyper.hyper_data{lt_run_index}.hyper_samp_delay; %Sample delay (s)

                %Dark correction data
                lt_dark_dat=all_hyper.hyper_data{lt_dark_run_index}.hyper_rad_dat;
                lt_dark_timer=all_hyper.hyper_data{lt_dark_run_index}.hyper_timer;
                lt_dark_time=all_hyper.hyper_data{lt_dark_run_index}.hyper_time;
                lt_dark_int=all_hyper.hyper_data{lt_dark_run_index}.hyper_int_time; %Sampling integration time (s)
                lt_dark_samp_delay=all_hyper.hyper_data{lt_dark_run_index}.hyper_samp_delay; %Sample delay (s)

                %Dark correct radiance data
                %Check to see if there are at least two points for interpolation
                if size(lt_dark_timer,1)>1
                    lt_dark_interp=interp1(lt_dark_timer,lt_dark_dat,lt_timer,'linear','extrap'); %Interpolation to match times
                else
                    error(''); %Error capture - comment out if desired and program will skip to next run
                    continue
                end

                lt_dat_dcorr=lt_dat-lt_dark_interp;  %Dark corrected upwelling radiance

                %Get minimum and maximum indices for lt_timer
                min_lt_index=find(lt_timer>=es_timer(min_index));
                if isempty(min_lt_index)==1  %Skip run if no matches
                    disp('No overlap in lt_timer and es_timer - advancing to next run');
                    continue
                end
                min_lt_index=min_lt_index(1);
                max_lt_index=find(lt_timer<=es_timer(max_index));
                max_lt_index=max_lt_index(end);

                %Test for no data
                if isempty(max_lt_index)==1
                    error(''); %Error capture - comment out if desired and program will skip to next run
                    disp('Advancing to next run due to no overlapping lt and es data...');
                    continue
                end

                max_lt_index=max_lt_index(end);

                %Read data for sky radiance
                li_dat=all_hyper.hyper_data{li_run_index}.hyper_rad_dat;
                li_timer=all_hyper.hyper_data{li_run_index}.hyper_timer;
                li_time=all_hyper.hyper_data{li_run_index}.hyper_time;
                li_int=all_hyper.hyper_data{li_run_index}.hyper_int_time; %Sampling integration time (s)
                li_samp_delay=all_hyper.hyper_data{li_run_index}.hyper_samp_delay; %Sample delay (s)

                %Dark correction data for sky radiance
                li_dark_dat=all_hyper.hyper_data{li_dark_run_index}.hyper_rad_dat;
                li_dark_timer=all_hyper.hyper_data{li_dark_run_index}.hyper_timer;
                li_dark_time=all_hyper.hyper_data{li_run_index}.hyper_time;
                li_dark_int=all_hyper.hyper_data{li_dark_run_index}.hyper_int_time; %Sampling integration time (s)
                li_dark_samp_delay=all_hyper.hyper_data{li_dark_run_index}.hyper_samp_delay; %Sample delay (s)

                %Test for no data
                if isempty(li_dark_dat)==1
                    error(''); %Error capture - comment out if desired and program will skip to next run
                    disp('Advancing to next run due to no data...');
                    continue
                end

                %Dark correct sky radiance data
                li_dark_interp=interp1(li_dark_timer,li_dark_dat,li_timer,'linear','extrap');  %Interpolation
                li_dat_dcorr=li_dat-li_dark_interp;
                li_750=li_dat_dcorr(:,lambda_li>748 & lambda_li<751);
                
                %Get sky radiance timer (li_timer) min and max indices in
                % relation to es_timer
                if isempty(li_timer)~=1 && isempty(find(li_timer>=es_timer(min_index),1))~=1
                    min_li_index=find(li_timer>=es_timer(min_index));
                    min_li_index=min_li_index(1);
                    max_li_index=find(li_timer<=es_timer(max_index));
                    max_li_index=max_li_index(end);
                end

                %Get GPS coordinates and times
                hyper_index=strcmp(all_hyper.hyper_sensor(run_index),gps_sensor)==1;  
                gps_run_index=run_index(hyper_index);
                gps_lat=all_hyper.hyper_data{gps_run_index}.gps_latdeg + all_hyper.hyper_data{gps_run_index}.gps_latmin./60;
                gps_lon=all_hyper.hyper_data{gps_run_index}.gps_londeg + all_hyper.hyper_data{gps_run_index}.gps_lonmin./60;
                gps_hour=all_hyper.hyper_data{gps_run_index}.gps_hour;
                gps_min=all_hyper.hyper_data{gps_run_index}.gps_min;
                gps_sec=all_hyper.hyper_data{gps_run_index}.gps_sec;
                gps_date=all_hyper.hyper_data{gps_run_index}.gps_dat{10};
                gps_datenum=all_hyper.hyper_data{gps_run_index}.gps_datenum + (gps_hour+gps_min./60+gps_sec./3600)./24;

                gps_course=all_hyper.hyper_data{gps_run_index}.gps_course;
                gps_speed_index=find(strncmp('SPEED',all_hyper.hyper_data{gps_run_index}.gps_hdr,5)==1);
                gps_speed=all_hyper.hyper_data{gps_run_index}.gps_dat{gps_speed_index};

                %Convert from knots to m/s
                gps_speed=gps_speed./1.944;
                
                gps_time=datestr([repmat([yearnum 1 1],size(gps_hour)), gps_hour, gps_min, gps_sec],13);
                %Get indices for gps times within the es_time interval
                gps_index=find(datenum(gps_time)>=datenum(es_time(min_index,:)) & datenum(gps_time)<=datenum(es_time(max_index,:)));

                %Get pitch and roll data
                hyper_index=find(strcmp(all_hyper.hyper_sensor(run_index),pitchroll_sensor)==1);
                ths_run_index=run_index(hyper_index);
                ths_pitch=all_hyper.hyper_data{ths_run_index}.ths_pitch;
                ths_roll=all_hyper.hyper_data{ths_run_index}.ths_roll;
                ths_timer=all_hyper.hyper_data{ths_run_index}.ths_timer;
                ths_time=all_hyper.hyper_data{ths_run_index}.ths_time;
                ths_comp=all_hyper.hyper_data{ths_run_index}.ths_comp; %Compass heading in degrees
                
                min_ths_index=find(ths_timer>=es_timer(min_index));
                min_ths_index=min_ths_index(1);
                max_ths_index=find(ths_timer<=es_timer(max_index));
                max_ths_index=max_ths_index(end);

                ths_index=find(ths_timer>es_timer(min_index) & ths_timer<es_timer(max_index));

                %******Correct radiance for surface reflectance - two approaches***** 
                %Mobley, C. D. (2015), Polarized reflectance and transmittance properties of windblown sea surfaces, Applied optics, 54(15), 4828-4849.
                %Lange, P. K., P. J. Werdell, Z. K. Erickson, G. Dall’Olmo, R. J. Brewin, M. V. Zubkov, G. A. Tarran, H. A. Bouman, W. H. Slade, and S. E. Craig (2020), Radiometric approach for the detection of picophytoplankton assemblages across oceanic fronts, Optics Express, 28(18), 25682-25705.

                disp('Computing surface reflectance...');

                %Get wind data for corresponding time period using read subroutine
                %[met_hdr,met_units,met_datenum,met_wspd,met_wdir]=read_ndbc_buoy_met_data(met_file_name);
                [met_hdr,met_units,met_datenum,met_wspd,met_wdir]=read_ndbc_buoy_met_data_v2(met_file_name);  %Nov

                %Interpolated windspeed for Ruddick et al. (2006) reflectance calculation
                ws_interp=interp1(met_datenum,met_wspd,gps_datenum(gps_index));

                %Find windspeed for met_datenum closest to es_datenum and
                % match to Mobley (2015) look-up values
                buoy_wspd=zeros(size(gps_index,1),1);
                for jk=1:length(gps_index)
                    [~,wsmin]=min(abs(met_datenum-gps_datenum(gps_index(1))));  %Use initial time for run - assume wind speed constant for run
                    buoy_wspd(jk)=met_wspd(wsmin);
                end
                
                possible_windspeed=[0;2;4;5;6;8;10;12;14;15];
                [~,wsmin2]=min(abs(possible_windspeed-buoy_wspd(1)));
                WS=possible_windspeed(wsmin2);    %Wind speed (m/s)

                %Compute time values, lat and lon, and solar zenith angles for the run time interval  
                DOY=gps_datenum-datenum(['1-Jan-',num2str(yearnum)])+1;
                [Y,M,D,H,MN,S] = datevec(gps_datenum(gps_index(1)));
                Lon=-gps_lon(gps_index);
                Lat=gps_lat(gps_index);

                %Matlab function to calculate solar position
                % Darin Koblick (2021). Vectorized Solar Azimuth and Elevation Estimation 
                % (https://www.mathworks.com/matlabcentral/fileexchange/23051-vectorized-solar-azimuth-and-
                % elevation-estimation), MATLAB Central File Exchange. Retrieved January 12, 2021.
                [Az,El] = SolarAzEl(datestr(gps_datenum(gps_index),'yyyy/mm/dd HH:MM:SS'),Lat,Lon,0);
                
                AZ_deg=Az;  %Solar azimuth in degress
                Theta_deg=90-El;  %Solar zenith angle in degrees
                
                %Reflectance correction using look-up table per Mobley (2015)
                polar_viewing_angle=35;  %35 degrees

                %Possible orientations for HyperSAS sensors (mounted
                % forward near bow on either side depending on sun position and 
                % directed perpendicular to the ship course track)
                
                %Here either specify orientation or assume sensors 
                % mounted on side of ship giving best relative solar azimuth angle
                % Desired optimum is between 50 and 170 degrees (Lange et al., 2020)
                
                azimuth_viewing_angle=zeros(length(gps_index),1);
                rho=zeros(length(gps_index),1);
                possible_zenith=[0:10:80,87.5];
                possible_polar_view=[0:10:80,87.5];
                [~,pmin]=min(abs(possible_polar_view-polar_viewing_angle));
                polar_view=possible_polar_view(pmin);
                possible_azimuth_view=0:15:180;
                
                %Acceptable range for azimuth_viewing angle
                max_rel_az_angle=170;
                min_rel_az_angle=50;

                %Compute li_750/es_750 ratio for Ruddick et al. (2006)
                %reflectance calculation
                %Interpolate to match gps times
                li_750_interp=interp1(datenum(li_time),li_750,datenum(gps_time(gps_index,:)));
                es_750_interp=interp1(es_test_datenum,es_750,datenum(gps_time(gps_index,:)));
                
                li_es_750_ratio=li_750_interp./es_750_interp;
                
                for ig=1:length(gps_index)
                 
                    %Azimuth view angle calculation (finds difference on scale of 0:180)
                    %IF(ABS($A2-(B2+90))>180,ABS(360-ABS($A2-(B2+90))),ABS($A2-(B2+90)))
                    %Starboard orientation
                    if abs(AZ_deg(ig)-(gps_course(gps_index(ig))+90))>180
                        starboard_azimuth_viewing_angle=abs(360-abs(AZ_deg(ig)-(gps_course(gps_index(ig))+90)));
                    else
                        starboard_azimuth_viewing_angle=abs(AZ_deg(ig)-(gps_course(gps_index(ig))+90));
                    end

                    %Port orientation (for GC2, Port orientation all days except 118)
                    %=IF(ABS($A2-(B2-90))>180,ABS(360-ABS($A2-(B2-90))),ABS($A2-(B2-90)))
                    if abs(AZ_deg(ig)-(gps_course(gps_index(ig))-90))>180
                        port_azimuth_viewing_angle=abs(360-abs(AZ_deg(ig)-(gps_course(gps_index(ig))-90)));
                    else
                        port_azimuth_viewing_angle=abs(AZ_deg(ig)-(gps_course(gps_index(ig))-90));
                    end
                    
                   if file_no==9 %file_no=9 corresponds to day 118
                       azimuth_viewing_angle(ig)=starboard_azimuth_viewing_angle;
                   else
                       azimuth_viewing_angle(ig)=port_azimuth_viewing_angle;
                   end
                    
                   %This section is commented out. This will determine which
                   % of the two orientations provides the largest angle
                   % within the max and min acceptable relative solar azimuth angles
                   %if (port_azimuth_viewing_angle<max_rel_az_angle && port_azimuth_viewing_angle>min_rel_az_angle) &&...
                   %     (starboard_azimuth_viewing_angle<max_rel_az_angle && starboard_azimuth_viewing_angle>min_rel_az_angle)
                   %         azimuth_viewing_angle(ig)=max([starboard_azimuth_viewing_angle,port_azimuth_viewing_angle]);
                   % elseif starboard_azimuth_viewing_angle<max_rel_az_angle && starboard_azimuth_viewing_angle>min_rel_az_angle
                   %     azimuth_viewing_angle(ig)=starboard_azimuth_viewing_angle;
                   % elseif port_azimuth_viewing_angle<max_rel_az_angle && port_azimuth_viewing_angle>min_rel_az_angle
                   %     azimuth_viewing_angle(ig)=port_azimuth_viewing_angle;
                   % else
                   %     A=[starboard_azimuth_viewing_angle,port_azimuth_viewing_angle];
                   %     [max_dif,xmin]=min(abs(A-max_rel_az_angle));
                   %     [min_dif,nmin]=min(abs(A-min_rel_az_angle));
                   %     if max_dif>min_dif
                   %         azimuth_viewing_angle(ig)=A(nmin);
                   %     else
                   %         azimuth_viewing_angle(ig)=A(xmin);
                   %     end
                   % end

                    %Match solar angles to possible values for Mobley (2015) look-up table
                    [~,zmin]=min(abs(possible_zenith-Theta_deg(ig)));
                    [~,amin]=min(abs(possible_azimuth_view-azimuth_viewing_angle(ig)));

                    zenith_input=possible_zenith(zmin);
                    azimuth_view=possible_azimuth_view(amin);

                    if ig==1
                        %Call Mobley subroutine to read look-up table for reflectance correction
                        rho_matrix=mobley_file_read(rho_filename);
                    end
                        
                        %Find target rho
                        rho_index=find(rho_matrix(:,1)==WS);
                        rho_indx2=find(rho_matrix(rho_index,2)==zenith_input);
                        rho_indx3=find(rho_matrix(rho_indx2,3)==polar_view);
                        rho_indx4=find(rho_matrix(rho_indx3,4)==azimuth_view);
                        mobley_rho=rho_matrix(rho_index(rho_indx2(rho_indx3(rho_indx4))),5);
                        
                        %MOBLEY CORRECTION COMMENTED OUT
                        %Set a maximum rho (to avoid overcorrection for
                        %reflectance)
                        %if mobley_rho>=max_rho
                        %    rho(ig)=max_rho;
                        %else
                        %    rho(ig)=mobley_rho; %surface reflectance of sky radiance
                        %end
                end

                %rho=0.028;  %After Mobley(1999)
                
                %THIS REFLECTANCE IS BEING USED IN PLACE OF MOBLEY(2015)
                %Test for cloudy vs clear sky conditions for Ruddick et al. (2006) correction
                rho=repmat(0.0256,length(gps_index),1);  %Cloudy
                rho(li_es_750_ratio<0.05)=0.0256+ws_interp(li_es_750_ratio<0.05).*0.00039+...
                    0.000034.*ws_interp(li_es_750_ratio<0.05).^2; %Clear sky - adjusts for windspeed

                %Correct radiance data for sky radiance reflectance
                disp('Applying sky correction...');
                if isempty(li_dat)~=1
                    %Interpolation to match wavelengths to Lt
                    li_dat_linterp=interp1(lambda_li,li_dat_dcorr',lambda_lt,'linear','extrap')'; 
                    
                    %Eliminate duplicate values
                    [C2,ic,id]=unique(li_timer);

                    %Interpolate to match lt_timer values and correct lt data
                    li_dat_interp=interp1(li_timer(ic),li_dat_linterp(ic,:),lt_timer,'linear'); %Interpolated sky reflected radiance
                    rho_interp=interp1(datenum(gps_time(gps_index,:)),rho,datenum(lt_time),'linear','extrap');
                    li_dat_rcorr=rho_interp.*li_dat_interp;   %Compute reflected sky radiance 
                    lt_dat_corr=lt_dat_dcorr-li_dat_rcorr;  %Reflectance water leaving radiance

                    %For alternative reflectance/NIR correction after Lange
                    % et al. (2020) - See below
                    li_dat_dcorr_linterp=interp1(lambda_li,li_dat_dcorr',lambda_lt,'linear','extrap')'; %Sky radiance interpolated to Lt wavelengths
                    li_dat_dcorr_interp=interp1(li_timer(ic),li_dat_dcorr_linterp(ic,:),lt_timer,'linear'); %Interpolated sky radiance
                end

                %*********************DISCRETE SENSOR PROCESSING**************************
                if discrete_process==1
                    %Load discrete sensor data to determine valid time intervals within 60 sample bins
                        hyper_index= strcmp(all_hyper.hyper_sensor(run_index),'DR7068e')==1;
                        dlt_run_index=run_index(hyper_index);
                        disp('Getting discrete Lt indices...');
                        dlt_index=find(strncmp(all_hyper.hyper_data{dlt_run_index}.discrete_hdr,'Lt',2)==1);
                        lam_hdr=all_hyper.hyper_data{dlt_run_index}.discrete_hdr(dlt_index);
                        [m,~]=size(lam_hdr);

                        lambda_dlt=zeros(m,1);
                        for j=1:m
                            str=char(lam_hdr(j));
                            lambda_dlt(j)=str2double(str(4:8));
                        end

                        dlt_dat=all_hyper.hyper_data{dlt_run_index}.discrete_rad_dat;
                        dlt_time=all_hyper.hyper_data{dlt_run_index}.discrete_time;
                        dlt_timer=all_hyper.hyper_data{dlt_run_index}.discrete_timer;
                        dlt_samp_delay=all_hyper.hyper_data{dlt_run_index}.discrete_samp_delay; %Sample delay (s)

                        %dlt_dat_interp=interp1(dlt_timer,dlt_dat,es_timer);
                        % Commented out because of too few data to interpolate

                        %Find dlt_timer indices within es_timer range
                        min_dlt_index=find(dlt_timer>=es_timer(min_index));
                        min_dlt_index=min_dlt_index(1);
                        max_dlt_index=find(dlt_timer<=es_timer(max_index));
                        max_dlt_index=max_dlt_index(end);

                        %Load discrete irradiance data (des)
                        hyper_index= strcmp(all_hyper.hyper_sensor(run_index),'DI7130D')==1;
                        des_run_index=run_index(hyper_index);
                        disp('Getting discrete Es indices...');
                        des_index=find(strncmp(all_hyper.hyper_data{des_run_index}.discrete_hdr,'ES',2)==1);
                        lam_hdr=all_hyper.hyper_data{des_run_index}.discrete_hdr(des_index);
                        [m,~]=size(lam_hdr);

                        lambda_des=zeros(m,1);
                        for j=1:m
                            str=char(lam_hdr(j));
                            lambda_des(j)=str2double(str(4:8));
                        end

                        des_dat=all_hyper.hyper_data{des_run_index}.discrete_rad_dat(:,1:4);  %Only UV wavelengths measured for discrete Es
                        des_time=all_hyper.hyper_data{des_run_index}.discrete_time;
                        des_timer=all_hyper.hyper_data{des_run_index}.discrete_timer;
                        des_samp_delay=all_hyper.hyper_data{des_run_index}.discrete_samp_delay; %Sample delay (s)

                        %Load discrete sky radiance data
                        hyper_index=find(strcmp(all_hyper.hyper_sensor(run_index),'DR7069d')==1);
                        dli_run_index=run_index(hyper_index);
                        disp('Getting discrete Li indices...');
                        dli_index=find(strncmp(all_hyper.hyper_data{dli_run_index}.discrete_hdr,'Li',2)==1);
                        lam_hdr=all_hyper.hyper_data{dli_run_index}.discrete_hdr(dli_index);
                        [m, n]=size(lam_hdr);

                        lambda_dli=zeros(m,1);
                        for j=1:m
                            str=char(lam_hdr(j));
                            lambda_dli(j)=str2double(str(4:8));
                        end

                        dli_dat=all_hyper.hyper_data{dli_run_index}.discrete_rad_dat;
                        dli_time=all_hyper.hyper_data{dli_run_index}.discrete_time;
                        dli_timer=all_hyper.hyper_data{dli_run_index}.discrete_timer;
                        dli_samp_delay=all_hyper.hyper_data{dli_run_index}.discrete_samp_delay; %Sample delay (s)

                       %Correct discrete radiance data
                        disp('Applying sky correction to discrete radiance...');

                        if isempty(dli_dat)
                            disp('Exiting due to no good data...');
                            continue
                        end

                        %Compute reflected sky discrete radiance using Mobley (2015)
                        dli_dat_linterp=interp1(lambda_dli,dli_dat',lambda_dlt,'linear','extrap')'; %Match wavelengths

                        %Eliminate duplicate values
                        [C3,ie,im]=unique(dli_timer);

                        %Interpolate to match dlt_timer values and correct dlt data
                        dli_dat_interp=interp1(dli_timer(ie),dli_dat_linterp(ie,:),dlt_timer,'linear');
                        rhod_interp=interp1(datenum(gps_time(gps_index,:)),rho,datenum(dlt_time));

                        %Interpolated sky reflected radiance
                        dli_dat_rcorr=rhod_interp.*dli_dat_interp;   %Compute reflected sky radiance 

                        dlt_dat_corr=dlt_dat-dli_dat_rcorr;   %Correct dlt data for reflected sky radiance

                        %Limit data to lowest 5% of values to reduce glint contamination effects 
                        dlt_good_index=[];
                        for j=min_dlt_index(1):60:max_dlt_index(1)-60
                            index_array=(j:j+60)';
                            dlt_include=index_array;
                            [dlt_sort,ix]=sort(dlt_dat_corr(dlt_include,5),1);
                            dlt_include=dlt_include(ix(1:3));  %Limit to lowest 5% of nir (680 nm) data
                            index_array=sort(dlt_include);
                            dlt_good_index=[dlt_good_index;index_array];
                        end

                        %Remove duplicate entries
                        dlt_good_index=sort(unique(dlt_good_index));            
                end
                %**********END DISCRETE PROCESSING*************************

                %*****BEGIN RRS CALCUATION******

                %Eliminate lt data with negative slope exceeding threshold in 750-800 nm range
                slope_test_wl=lambda_lt(lambda_lt>770 & lambda_lt<800);
                lt_test=lt_dat_corr(:,lambda_lt>770 & lambda_lt<800)./lt_dat_corr(:,lambda_lt>769 & lambda_lt<771);
                %lt_test=lt_dat_corr(:,lambda_lt>770 & lambda_lt<800)./abs(lt_dat_corr(:,lambda_lt>769 & lambda_lt<771));
                lt_test_indx=[];
                test_m_all=zeros(length(min_lt_index(1):max_lt_index(1)),1);
                for it=min_lt_index(1):max_lt_index(1)
                    [test_m,test_b,test_r,test_smx,test_sbx]=lsqfitx(slope_test_wl,lt_test(it,:)');
                    test_m_all(it)=test_m;
                    if test_m>-0.007 && test_m<0.01
                    %if test_m>-0.006 && test_m<0.006
                        lt_test_indx=[lt_test_indx;it];
                    else
                        continue
                    end
                end
                
                %This statement will bypass the slope test above
                %lt_test_indx=find(lt_test(:,1));
                
                %Limit data to lowest 10% of spectra for selected near infrared range 
                % to reduce glint contamination effects 
                num_int=60;  %Interval over which lowest values selected
                lt_good_index=[];
                for jt=1:num_int:length(lt_test_indx)-num_int %min_lt_index(1):60:max_lt_index(1)-60  %Evaluate every 60 observations
                    index_array=(jt:jt+num_int)';
                    lt_include=lt_test_indx(index_array);      %index_array;
                    [lt_sort,ix]=sort(median(lt_dat_corr(lt_include,nir_index),2),1);  %Sort data by median value within nir specified range
                    lt_include=lt_include(ix(1:6));  %Limit to lowest 10% of data
                    index_array=sort(lt_include);
                    lt_good_index=[lt_good_index;index_array];
                end
                
                %Exit loop if no good data
                if isempty(lt_good_index)
                    disp('Exiting due to no good data...');
                    continue
                end

                %Remove duplicate entries
                lt_good_index=sort(unique(lt_good_index));            

                %For alternative reflectance calculation - cost minimization function for 750-800 nm (Lange et al., 2020)
                nir_lt=lt_dat_dcorr(lt_good_index,alt_nir_index);
                nir_li=li_dat_dcorr_interp(lt_good_index,alt_nir_index);
                %Omit values with nan
                nir_lt=nir_lt(~isnan(nir_li));
                nir_li=nir_li(~isnan(nir_li));
                %Anonymous function for minization 
                fun = @(x)abs(nir_lt-x(1)*nir_li-x(2));  %x(1) is nir reflectance (nir_rho) and x(2) is residual
                x0 = [0.028,0.01];  %Initial values for rho and LNIR
                lb = [0.02,-0.5]; 
                ub = [0.1,5];
                options = optimoptions(@lsqnonlin); %Optimization options for lsqnonlin
                x = lsqnonlin(fun,x0,lb,ub); %Solve for values
                lt_dat_alt_corr=lt_dat_dcorr-x(1).*li_dat_dcorr_interp-x(2)+x(2);  %Compute corrected radiance
                dlt_dat_alt_corr=dlt_dat_corr-x(1).*dli_dat_interp-x(2)+x(2);  %Compute corrected radiance
                
                %Find closest matching indices in Es data corresponding to selected Lt times
                m=size(lt_good_index,1);
                %Find closest es_timer value less than given lt_timer value
                es_good_index=[];
                for j=1:m
                    t=es_timer<lt_timer(lt_good_index(j));
                    [mx,ix] = max(es_timer(t));    %mx is maximum value in es_timer less than lt_timer(lt_good_index(j)) and ix is its index with respect to es_timer.
                    f = find(t);
                    good_es_add = f(ix);
                    %good_es_add=find(abs(es_timer-lt_timer(lt_good_index(j)))<=.1);
                    es_good_index=[es_good_index;good_es_add];
                end

                if discrete_process==1
                    m=size(dlt_good_index,1);

                    %Find closest matching indices in des data corresponding to selected dlt times
                    des_good_index=[];
                    for j=1:m
                        t=des_timer<dlt_timer(dlt_good_index(j));  
                        [mx,ix] = max(des_timer(t));
                        f = find(t);
                        good_des_add = f(ix);
                        %good_des_add=find(abs(des_timer-dlt_timer(dlt_good_index(j)))<=.1);
                        des_good_index=[des_good_index;good_des_add];
                    end
                end

                %Find matching Es and LT (and Li) data for estimating rrs

                %Remove duplicate entries
                [c,ia,ib]=unique(es_timer(es_good_index));

                %Interpolate irradiance es_timer to radiance lt_timer
                es_dat_interp=interp1(es_timer(es_good_index(ia)),es_dat_dcorr(es_good_index(ia),:),...
                    lt_timer(lt_good_index),'linear','extrap');
                %es_datenum_interp=interp1(es_timer(es_good_index(ia)),datenum(es_time(es_good_index(ia),:)),...
                es_datenum_interp=interp1(es_timer(es_good_index(ia)),datenum(es_time(es_good_index(ia))),...
                    lt_timer(lt_good_index),'linear','extrap');
                
                %******Compute Rrs********
                rrs=lt_dat_corr(lt_good_index,:)./es_dat_interp;
                %For alternative NIR correction (Lange et al., 2020) - See above
                alt_rrs=lt_dat_alt_corr(lt_good_index,:)./es_dat_interp;
                
                if discrete_process==1
                    %For discrete Rrs calculations

                    %Remove duplicate entries (NOT NEEDED?)
                    [c,ia,ib]=unique(des_timer(des_good_index));

                    des_dat_interp=interp1(des_timer(des_good_index(ia)),des_dat(des_good_index(ia),:),...
                        dlt_timer(dlt_good_index),'linear','extrap');
                    drrs=dlt_dat_corr(dlt_good_index,1:4)./des_dat_interp;
                    alt_drrs=dlt_dat_alt_corr(dlt_good_index,1:4)./des_dat_interp;

                    des_datenum_interp=interp1(des_timer(des_good_index(ia)),datenum(des_time(des_good_index(ia))),...
                        dlt_timer(dlt_good_index),'linear','extrap');
                end

                %Plot irradiance and radiance at selected wavelengths for data quality checking
                clf('reset');
                disp('Plotting data...');
                hes=plot(es_timer(es_good_index),es_dat_dcorr(es_good_index,wave_index));
                hold on
                
                %Get time range for axis parameters 
                time_range=es_timer(max_index)-es_timer(min_index);

                %Set up axes and labelling
                incr=time_range.*.1;
                xlimit=[es_timer(min_index) es_timer(max_index)+incr];
                ylabel('E_s (\muW/cm^{2}/nm)','FontSize',10,'FontWeight','bold','FontName','Arial');
                set(gca,'FontSize',10,'FontName','Arial','YLimMode','Auto','Xlim', xlimit,'YAxisLocation','left',...
                 'XtickLabel',[],'Position',[0.1000    0.5838    0.3347    0.3412]);
                xtick_dat=get(gca,'XTick');

                xtick_es_timer=es_timer(min_index:max_index);
                xtick_es_datenum=es_datenum(min_index:max_index);

                htit=title([num2str(Y),' Day ',num2str(floor(DOY(1)),'%03.0f'),' Run ',num2str(lrun)],'FontSize',12);
                %tit_pos=get(htit,'Position');  %To adjust title position
                %tit_pos(1)=tit_pos(1)+500;
                %tit_pos(2)=tit_pos(2)+6;                
                %set(htit,'Position',tit_pos);

                %Eliminate duplicate values
                [C3,ie,~]=unique(xtick_es_timer);
                
                %Get date and time information for plotting by
                % interpolating datenum information for corresponding xtick es_timer values  
                xtick_datenum=interp1(xtick_es_timer(ie),xtick_es_datenum(ie),xtick_dat,'linear','extrap');
                hold on

                %xtick_lab=get(gca,'XtickLabel');
                xtick_lab=datestr(xtick_datenum,15);  %Convert timer times to clock strings; formatting for time only

                %Plot discrete sensor irradiance
                hdes=plot(des_timer(des_good_index),des_dat(des_good_index,4),':g');
                hold on

                ylimit=get(gca,'YLim');

                %Print wavelength on figure
                htext=text(xlimit(1)+0.1.*(xlimit(2)-xlimit(1)),ylimit(2)-0.1.*(ylimit(2)-ylimit(1)),num2str(lambda_es(wave_index)),...
                    'FontSize',10,'FontName','Arial');
                %xlimit(1)+0.1.*xlimit(2),ylimit(2)-0.1.*ylimit(2)
                hold on

                %Set up axes for radiance (right axis)
                ax_pos=get(gca,'Position');

                hrad=axes('YLimMode','auto','XTick',xtick_dat,'Xlim',xlimit,'XtickLabel',xtick_lab);
                set(hrad,'Visible','on','Color','none','YAxisLocation','right','FontSize',10,...
                    'FontWeight','normal','FontName','Arial','Position',ax_pos); %,'Ylim',[0 10]);
                hold on

                xlabel('Time (UTC)','FontSize',10,'FontWeight','bold','FontName','Arial');
                hold on

                %Plot discrete radiance
                hdlt=plot(dlt_timer(dlt_good_index),dlt_dat_corr(dlt_good_index,4),'m>-','MarkerSize',2); %
                hold on

                %Plot radiance
                hlt=plot(lt_timer(lt_good_index),lt_dat_corr(lt_good_index,wave_index),'r.:');
                hold on

                %Axis labels
                hylab=ylabel('L_T (\muW/cm^2/nm/sr)','FontSize',10,'FontWeight','bold','FontName','Arial',...
                    'Rotation',270);
                ylab_pos=get(hylab,'Position');
                ylab_pos=[ylab_pos(1)+0.05*time_range ylab_pos(2) ylab_pos(3)];
                set(hylab,'Position',ylab_pos);
                hold on

                hleg=legend([hes,hdes,hlt,hdlt],'Es(382)','Es-dis(380)','Lt(382)','Lt-dis(380)');
                set(hleg,'Fontsize',8);

                %Plot Rrs - second figure
                disp('Plotting Rrs...');

                hrrs_time=axes('FontSize',10,'FontName','Arial','YLimMode','Auto','Xlim', xlimit,'YAxisLocation','left',...
                 'XtickLabel',xtick_lab,'Position',[0.1000    0.13    0.3347    0.3412]);
                hold on

                hrrs=plot(lt_timer(lt_good_index),rrs(:,refl_index),'+b');
                hold on
                halt_rrs=plot(lt_timer(lt_good_index),alt_rrs(:,refl_index),'xm');
                %Set x limit
                xlimit=[es_timer(min_index) es_timer(max_index)+incr];

                ylabel('R_{RS} (sr^{-1})','FontSize',10,'FontWeight','bold','FontName','Arial');
                xlabel('Time (UTC)','FontSize',10,'FontWeight','bold','FontName','Arial');
                %xtick_dat=get(gca,'XTick');
                %xlimit=get(gca,'Xlim');
                %xtick_lab=get(gca,'XtickLabel');

                box on
                hold on

                hleg=legend([hrrs,halt_rrs],'Rrs(550)','Alt-RRS(550)');
                set(hleg,'Fontsize',8);

                %Set up axes for radiance comparison - third figure
                hrad=axes('YLimMode','auto','XTick',400:200:800,'Xlim',[290 810],'Position',[0.60    0.5838    0.3347    0.3412]);
                box on
                hold on

                hlt=plot(lambda_lt,lt_dat_dcorr,'b.');
                hold on

                hli=plot(lambda_lt,li_dat_rcorr,'r-');
                hold on

                hleg=legend([hlt(1),hli(1)],'Lt','Li correction');
                set(hleg,'Fontsize',8);

                hylab=ylabel('L_x (\muW/cm^2/nm/sr)','FontSize',10,'FontWeight','bold','FontName','Arial');
                xlabel('Wavelength(nm)','FontSize',10,'FontWeight','bold','FontName','Arial');

                %Plot Rrs Spectrum - fourth figure

                %Set up axes
                hrrs=axes('YLimMode','auto','XTick',400:200:800,'Xlim',[290 810],'Position',[0.60    0.13    0.3347    0.3412]);
                box on
                hold on

                set(hrrs,'Visible','on','Color','none','YAxisLocation','left','FontSize',10,...
                    'FontWeight','normal','FontName','Arial');%,'Ylim',[0 .03]); %,'XtickLabel',xtick_lab);
                hold on

                if alt_rrs_plot==1
                    plot(lambda_lt,alt_rrs);  %For alternative NIR reflectance correction
                else
                    plot(lambda_lt,rrs);
                end
                
                %Axis labels
                hylab=ylabel('R_{RS} (sr^{-1})','FontSize',10,'FontWeight','bold','FontName','Arial',...
                    'Rotation',90);
                xlabel('Wavelength(nm)','FontSize',10,'FontWeight','bold','FontName','Arial');
                hold on

                lt_datenum=interp1(xtick_es_timer(ie),xtick_es_datenum(ie),lt_timer);
                %lt_time=cellstr(datestr(lt_datenum,13));

                %Plot discrete Rrs

                if alt_rrs_plot==1
                    plot(lambda_dlt(1:4),alt_drrs,'o');  %For alternative NIR reflectance correction
                else
                    plot(lambda_dlt(1:4),drrs,'o');
                end
                hold on

                m=size(lt_good_index,1);

                %Find indices for closest gps field corresponding to rrs times
                gps_good_index=zeros(m,1);
                for j=1:m
                    [diff]=min(abs(datenum(gps_time)-es_datenum_interp(j)));
                    good_gps_add=find(abs(datenum(gps_time)-es_datenum_interp(j))==diff);
                    gps_good_index(j)=good_gps_add(1);
                end

                rrs_lat=gps_lat(gps_good_index);
                rrs_lon=gps_lon(gps_good_index);

                %Create structured output variable
                if ~exist('rad_rrs_output','var')
                    disp('Creating output variable');
                    rad_rrs_output=struct('run_no',lrun,'es_time',{{es_time(min_index:max_index)}},...
                        'es_timer',{es_timer(min_index:max_index)},'es_datenum',{es_datenum(min_index:max_index)},...
                        'es_data_interp',{es_dat_interp},'es_data',{es_dat_dcorr(min_index:max_index,:)},...
                        'es_int_time',{es_int(min_index:max_index)},'es_samp_delay',{es_samp_delay(min_index:max_index)},...
                        'lambda_es',{lambda_es},'start_datenum',es_datenum(min_index),...
                        ...
                        'lt_timer_all',{lt_timer},'lt_data_all',{lt_dat_dcorr},'lt_datenum_all',{lt_datenum},...
                        'lt_timer',{lt_timer(lt_good_index)},'lt_time',{{lt_time(lt_good_index)}},'lt_datenum',{lt_datenum(lt_good_index)},...
                        'lt_int_time',{lt_int(lt_good_index)},'lt_samp_delay',{lt_samp_delay(lt_good_index)},...
                        'lambda_lt',{lambda_lt},'lt_data',{lt_dat_dcorr(lt_good_index,:)},'lt_data_ref_corr',{lt_dat_corr(lt_good_index,:)},...
                        'alt_lt_data',{lt_dat_alt_corr(lt_good_index,:)},...
                        ...
                        'lambda_li',{lambda_li},'li_data',{li_dat_dcorr(min_li_index:max_li_index,:)},...
                        'li_dat_interp',{li_dat_interp(lt_good_index,:)},'li_time',{{li_time(min_li_index:max_li_index)}},...
                        'li_timer',{li_timer(min_li_index:max_li_index)},...
                        'li_int_time',{li_int(min_li_index:max_li_index)},'li_samp_delay',{li_samp_delay(min_li_index:max_li_index)},...
                        'rrs',{rrs},'alt_rrs',{alt_rrs},...
                        ...
                        'dlt_timer',{dlt_timer(dlt_good_index)},'dlt_time',{{dlt_time(dlt_good_index)}},...
                        'lambda_dlt',{lambda_dlt},'dlt_data',{dlt_dat_corr(dlt_good_index,:)},...
                        'alt_dlt_data',{dlt_dat_alt_corr(dlt_good_index,:)},...
                        'dlt_samp_delay',{dlt_samp_delay(dlt_good_index)},'drrs',{drrs},'alt_drrs',{alt_drrs},...
                        ...
                        'des_time',{{des_time(des_good_index)}},'des_timer',{des_timer(des_good_index)},'des_data',{des_dat(des_good_index,:)},...
                        'des_samp_delay',{des_samp_delay(des_good_index)},'des_datenum_interp',{des_datenum_interp},'des_data_interp',{des_dat_interp},...
                        ...
                        'ths_pitch',{ths_pitch(ths_index)},'ths_roll',{ths_roll(ths_index)},...
                        'ths_timer',{ths_timer(ths_index)},'ths_time',{{ths_time(ths_index)}},'ths_comp',{ths_comp(ths_index)},...
                        ...
                        'solazi',{AZ_deg},'solelev',{El},'azimuth_viewing_angle',azimuth_viewing_angle,...
                        'wind_speed',{buoy_wspd},'rho',{rho},'nir_rho',{x(1)},'lnir',{x(2)},'ship_ws',ship_ws,...
                        'ship_cloud_cover',ship_cloud_cover,'ship_waveh',ship_waveh,...
                        ...
                        'gps_course',{gps_course(gps_index)},'gps_time',{gps_time(gps_index,:)},'gps_datenum',{gps_datenum(gps_index)},...
                        'gps_speed',{gps_speed(gps_index)},'gps_lat',{gps_lat(gps_index)},'gps_lon',{gps_lon(gps_index)},...
                        'rrs_lat',{rrs_lat},'rrs_lon',{rrs_lon});
                else
                    disp('Adding data to output file');
                    rad_rrs_output.run_no=[rad_rrs_output.run_no;lrun];
                    rad_rrs_output.es_time=[rad_rrs_output.es_time;{es_time(min_index:max_index)}];
                    rad_rrs_output.es_timer=[rad_rrs_output.es_timer;{es_timer(min_index:max_index)}];
                    rad_rrs_output.es_int_time=[rad_rrs_output.es_int_time;{es_int(min_index:max_index)}];
                    rad_rrs_output.es_samp_delay=[rad_rrs_output.es_samp_delay;{es_samp_delay(min_index:max_index)}];
                    rad_rrs_output.es_datenum=[rad_rrs_output.es_datenum;{es_datenum(min_index:max_index)}];
                    rad_rrs_output.es_data_interp=[rad_rrs_output.es_data_interp;{es_dat_interp}];
                    rad_rrs_output.es_data=[rad_rrs_output.es_data;{es_dat_dcorr(min_index:max_index,:)}];
                    
                    rad_rrs_output.lt_timer_all=[rad_rrs_output.lt_timer_all;{lt_timer}];
                    rad_rrs_output.lt_data_all=[rad_rrs_output.lt_data_all;{lt_dat_dcorr}];
                    rad_rrs_output.lt_datenum_all=[rad_rrs_output.lt_datenum_all;{lt_datenum}];
                    rad_rrs_output.lt_timer=[rad_rrs_output.lt_timer;{lt_timer(lt_good_index)}];
                    rad_rrs_output.lt_data=[rad_rrs_output.lt_data;{lt_dat_dcorr(lt_good_index,:)}];
                    rad_rrs_output.lt_data_ref_corr=[rad_rrs_output.lt_data_ref_corr;{lt_dat_corr(lt_good_index,:)}];
                    rad_rrs_output.alt_lt_data=[rad_rrs_output.alt_lt_data;{lt_dat_alt_corr(lt_good_index,:)}];
                    rad_rrs_output.lt_datenum=[rad_rrs_output.lt_datenum;{lt_datenum(lt_good_index)}];
                    rad_rrs_output.lt_time=[rad_rrs_output.lt_time;{lt_time(lt_good_index)}];
                    rad_rrs_output.lt_int_time=[rad_rrs_output.lt_int_time;{lt_int(lt_good_index)}];
                    rad_rrs_output.lt_samp_delay=[rad_rrs_output.lt_samp_delay;{lt_samp_delay(lt_good_index)}];
                    rad_rrs_output.start_datenum=[rad_rrs_output.start_datenum;{es_datenum(min_index)}];

                    rad_rrs_output.li_timer=[rad_rrs_output.li_timer;{li_timer(min_li_index:max_li_index)}];
                    rad_rrs_output.li_data=[rad_rrs_output.li_data;{li_dat_dcorr(min_li_index:max_li_index,:)}];
                    rad_rrs_output.li_dat_interp=[rad_rrs_output.li_dat_interp;{li_dat_interp(lt_good_index,:)}];
                    rad_rrs_output.li_time=[rad_rrs_output.li_time;{li_time(min_li_index:max_li_index)}];
                    rad_rrs_output.li_int_time=[rad_rrs_output.li_int_time;{li_int(min_li_index:max_li_index)}];
                    rad_rrs_output.li_samp_delay=[rad_rrs_output.li_samp_delay;{li_samp_delay(min_li_index:max_li_index)}];

                    rad_rrs_output.rrs=[rad_rrs_output.rrs;{rrs}];
                    rad_rrs_output.alt_rrs=[rad_rrs_output.alt_rrs;{alt_rrs}];

                    rad_rrs_output.dlt_timer=[rad_rrs_output.dlt_timer;{dlt_timer(dlt_good_index)}];
                    rad_rrs_output.dlt_time=[rad_rrs_output.dlt_time;{dlt_time(dlt_good_index)}];
                    rad_rrs_output.dlt_samp_delay=[rad_rrs_output.dlt_samp_delay;{dlt_samp_delay(dlt_good_index)}];
                    rad_rrs_output.dlt_data=[rad_rrs_output.dlt_data;{dlt_dat_corr(dlt_good_index,:)}];
                    rad_rrs_output.alt_dlt_data=[rad_rrs_output.alt_dlt_data;{dlt_dat_alt_corr(dlt_good_index,:)}];
                    
                    rad_rrs_output.des_time=[rad_rrs_output.des_time;{des_time(des_good_index)}];
                    rad_rrs_output.des_timer=[rad_rrs_output.des_timer;{des_timer(des_good_index)}];
                    rad_rrs_output.des_data=[rad_rrs_output.des_data;{des_dat(des_good_index,:)}];
                    rad_rrs_output.des_samp_delay=[rad_rrs_output.des_samp_delay;{des_samp_delay(des_good_index)}];
                    rad_rrs_output.des_data_interp=[rad_rrs_output.des_data_interp;{des_dat_interp}];
                    rad_rrs_output.des_datenum_interp=[rad_rrs_output.des_datenum_interp;{des_datenum_interp}];
                   
                    rad_rrs_output.drrs=[rad_rrs_output.drrs;{drrs}];
                    rad_rrs_output.alt_drrs=[rad_rrs_output.alt_drrs;{alt_drrs}];
                    
                    rad_rrs_output.ths_pitch=[rad_rrs_output.ths_pitch;{ths_pitch(ths_index)}];
                    rad_rrs_output.ths_timer=[rad_rrs_output.ths_timer;{ths_timer(ths_index)}];
                    rad_rrs_output.ths_time=[rad_rrs_output.ths_time;{ths_time(ths_index)}];
                    rad_rrs_output.ths_roll=[rad_rrs_output.ths_roll;{ths_roll(ths_index)}];
                    
                    rad_rrs_output.solazi=[rad_rrs_output.solazi;{AZ_deg}];
                    rad_rrs_output.solelev=[rad_rrs_output.solelev;{El}];
                    rad_rrs_output.azimuth_viewing_angle=[rad_rrs_output.azimuth_viewing_angle;{azimuth_viewing_angle}];
                    rad_rrs_output.wind_speed=[rad_rrs_output.wind_speed;{buoy_wspd}];
                    rad_rrs_output.rho=[rad_rrs_output.rho;{rho}];
                    rad_rrs_output.nir_rho=[rad_rrs_output.nir_rho;{x(1)}];
                    rad_rrs_output.lnir=[rad_rrs_output.lnir;{x(2)}];
                    
                    rad_rrs_output.ship_ws=[rad_rrs_output.ship_ws;ship_ws];
                    rad_rrs_output.ship_waveh=[rad_rrs_output.ship_waveh;ship_waveh];
                    rad_rrs_output.ship_cloud_cover=[rad_rrs_output.ship_cloud_cover;ship_cloud_cover];

                    rad_rrs_output.gps_course=[rad_rrs_output.gps_course;{gps_course(gps_index)}];
                    rad_rrs_output.gps_time=[rad_rrs_output.gps_time;{gps_time(gps_index,:)}];
                    rad_rrs_output.gps_datenum=[rad_rrs_output.gps_datenum;{gps_datenum(gps_index)}];
                    rad_rrs_output.gps_speed=[rad_rrs_output.gps_speed;{gps_speed(gps_index)}];
                    rad_rrs_output.gps_lat=[rad_rrs_output.gps_lat;{gps_lat(gps_index)}];
                    rad_rrs_output.gps_lon=[rad_rrs_output.gps_lon;{gps_lon(gps_index)}];
                    rad_rrs_output.rrs_lat=[rad_rrs_output.rrs_lat;{rrs_lat}];
                    rad_rrs_output.rrs_lon=[rad_rrs_output.rrs_lon;{rrs_lon}];
                end

                %if strcmp(utccorr,'True')==1
                    time_lab=strrep(es_time{1}(1:8),':','');
                %else
                %    time_lab=strrep(es_time{1}(1,1:8),':','');
                %end

                %Statements to print graphs
                %eval(['print ''',inpath,'rad_rrs_plot_',num2str(Y),'_',num2str(floor(DOY(1)),'%03.0f'),'_run',num2str(lrun),'_',time_lab,'.jpg'' -djpeg -r300']);
                %eval(['print ''',inpath,'rad_rrs_plot_Es253_',num2str(Y),'_',num2str(floor(DOY(1)),'%03.0f'),'_run',num2str(lrun),'_',time_lab,'.jpg'' -djpeg -r300']); 

            end

                %Statements to print output file
                if exist('rad_rrs_output','var')
                    %%eval(['save ',inpath,'rad_rrs_',outfile_label,'_run',num2str(lrun),'.mat rad_rrs_output']); %For individual run files
                    %%eval(['save ',inpath,'rad_rrs_',outfile_label,'.mat rad_rrs_output']);
                    %eval(['save ',inpath,'rad_rrs_Es253_',outfile_label,'.mat rad_rrs_output']);
                end

        end 
    end
    
    disp({'Execution Complete'});

% hypersas_rrs_Apr2009_output_to_seabass_v4.m - Program to read HyperSAS processed data
% and output to NASA SeaBASS format
% 
% Syntax:  hypersas_rrs_Apr2009_output_to_seabass_v4.m
%
% Inputs:
%    1) Folder location with compiled data files ('*.mat')
%      produced by 'generate_hypersas_reflectance_radiance_files_all_cruise_rev6.m'
%    2) Spreadsheet with station metadata information 
%      ('apr2009_hypersas_plot_qaa_parameters_v2a.xls')
%
% Outputs:
%    output - '*.sb' files in NASA SeaBASS format
%   
% Other m-files required: 
%  1) The script 'hypersas_data_read.m' is required to read the
%      HyperSAS data files and load variables for processing
%  2) The script 'hypersas_data_filter.m' is required to apply filter and 
%      baseline corrections to HyperSAS reflectance
%
% MAT-files required: Uses '*.mat' files produced from 
%      'generate_hypersas_reflectance_radiance_files_all_cruise_rev6.m'
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 7 Sep 2021

%% ------------- BEGIN CODE --------------%

clear
clc
clf reset
close all

%Input folder name and file information
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
infile_name=[inpath,'apr2009_hypersas_plot_qaa_parameters_v2.xlsx'];

single_sta='Run10';  %Single station for processing just one

%Input parameters
[~,~,all_par_3]=xlsread(infile_name,1,'C3:R77');
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

run_indx=[1:4,6:14,16,18:23,25:28,30:31,33:46,48:49,52:56,58:61,63];
run_array=cell(1,length(run_indx));
for runn=run_indx(1):run_indx(length(run_indx))
    run_array(runn)={['Run',num2str(runn)]};
end
station_names=cat(2,{'A1','A1b','A2','A3','B4','C1','C5','E0','E1','E3','D5','F4','F5','MR2','MR3',...
'G2','G3','G4','G5','H1','H2','H3','UWS1','UWS2'},run_array);

[nsta,~]=size(hyper_sta);

%hfig=figure(1);

%% Begin loop to process each station and run 

for im=1:nsta %1:nsta 

    sta_txt=hyper_sta{im};  %Select station to plot

    %To do all stations, comment out the continue statement
    if strcmp(single_sta,sta_txt)~=1
        %continue
    end

    %HyperSAS parameters
    file_no=hyper_doy(im);
    run_no=hyper_run(im);
    file_date=datestr(datenum(hyper_date(im,4:7),'yyyy')+str2num(hyper_date(im,8:10))-1,'yyyymmdd');
    file_label=hyper_date(im,:);
    time_str1=datestr(hyper_start_time(im),'HH:MM:SS');
    time_str2=datestr(hyper_end_time(im),'HH:MM:SS');

    station_txt=sta_lab{im};
    title_txt=title_lab{im};

    %Baseline correction flags
    base_corr=sas_basecorr_flag(im);    %HyperSAS base correction (1=on; 0=off)

    %file_name=[inpath,'rad_rrs_',num2str(file_no),'_',num2str(run_no),'.mat'];
    file_name=[inpath,'rad_rrs_Es253_',num2str(file_no),'.mat'];
    text_msg=['loading ',file_name,'...'];

    %% Read and filter data
    
    %Run matlab script to load file and populate variables
    hypersas_data_read;

    %Run matlab script to apply data filters and baseline corrections
    hypersas_data_filter;

    %Move on to next station if no acceptable data
    if isempty(good_lt_index)==1 || isempty(good_dlt_index)==1
        continue
    end
    
    %% Output data to seabass format 

    %Create field headings for radiometric wavelengths
    wavel_str=num2str(lt_lambda,'%4.1f');
    dwavel_str=num2str(drrs_lambda,'%4.1f');
    Lt_fields=[repmat(',Lt',137,1),wavel_str];
    Lsky_fields=[repmat(',Lsky',137,1),wavel_str];
    Lw_fields=[repmat(',Lw',137,1),wavel_str];
    Es_fields=[repmat(',Es',137,1),wavel_str];
    Rrs_fields=[repmat(',Rrs',137,1),wavel_str];

    %Discrete sensors
    dLt_fields=[repmat(',Lt',5,1),dwavel_str,repmat('_1id',5,1)];
    dLsky_fields=[repmat(',Lsky',5,1),dwavel_str,repmat('_1id',5,1)];
    dLw_fields=[repmat(',Lw',5,1),dwavel_str,repmat('_1id',5,1)];
    dEs_fields=[repmat(',Es',4,1),dwavel_str(1:4,:),repmat('_1id',4,1)];
    dRrs_fields=[repmat(',Rrs',4,1),dwavel_str(1:4,:),repmat('_1id',4,1)];

    %dLsky_fields
    %dLw_fields
    %dEs_fields
    
    fields=['/fields=date,time,lat,lon,heading,speed_f_w,SZA,SAZ,RelAz,wind,depth,pitch,roll',...
       reshape(Lt_fields',1,[]),reshape(Lsky_fields',1,[]),reshape(Lw_fields',1,[]),...
       reshape(Es_fields',1,[]),reshape(Rrs_fields',1,[]),reshape(dLt_fields',1,[]),...
       reshape(dLsky_fields',1,[]),reshape(dLw_fields',1,[]),...
       reshape(dEs_fields',1,[]),reshape(dRrs_fields',1,[])];
    units=['/units=yyyymmdd,hh:mm:ss,degrees,degrees,degrees,m/s,degrees,degrees,degrees,m/s,m,degrees,degrees'...
        repmat(',uW/cm^2/nm/sr',1,137.*3),repmat(',uW/cm^2/nm',1,137),repmat(',1/sr',1,137),...
        repmat(',uW/cm^2/nm/sr',1,5.*3),repmat(',uW/cm^2/nm',1,4),repmat(',1/sr',1,4)];

    %rad_data_file_name=GulfCarbon_gulfcarbon2_hypersas_St#_2009MMDD_R1.sb
    if contains(sta_txt,'UWS') || contains(sta_txt,'Run')
        rad_data_file_name=['GulfCarbon_gulfcarbon2_hypersas_',sta_txt,'_',file_label,'_R1.sb'];
    else
        rad_data_file_name=['GulfCarbon_gulfcarbon2_hypersas_St',sta_txt,'_',file_label,'_R1.sb'];
    end
        
    textstring=metafile_txt;
    [hdr_m]=size(char(textstring{1}));
    file_creation_date=datestr(datenum(date),'yyyymmdd');
    
    %Change longitude to negative to conform to SeaBASS convention
    rrs_lon_filt=-rrs_lon_filt;

    %% Input metadata into file header
    for i=1:hdr_m
        if ~contains(char(textstring{1}(i)),'/north_latitude')~=1
            textstring{1}(i)={['/north_latitude=',num2str(max(rrs_lat_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/south_latitude'))~=1 %#ok<*STREMP>
            textstring{1}(i)={['/south_latitude=',num2str(min(rrs_lat_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/west_longitude'))~=1
            textstring{1}(i)={['/west_longitude=',num2str(min(rrs_lon_filt)),'[DEG]']};
        elseif isempty(strfind(char(textstring{1}(i)),'/east_longitude'))~=1
            textstring{1}(i)={['/east_longitude=',num2str(max(rrs_lon_filt)),'[DEG]']};
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

    
    %% Format output 
    
    %Create dummy depth field to satisfy SeaBASS submission requirements
    rad_depth=zeros(size(rrs_lat_filt,1),1);

    anc_data_array=[rrs_lat_filt,rrs_lon_filt,gps_course_intrp_filt,gps_speed_intrp_filt,zenith_intrp_filt,...
        solazi_intrp_filt,az_view_intrp_filt,wind_speed_intrp_filt,rad_depth,pitch_intrp_filt,roll_intrp_filt];
    rad_data_array=[lt_filt,li_filt,lw_filt,es_filt,rrs_corr,dlt_filt,dli_filt,dlw_filt,des_filt,drrs_corr];
    
    %Drop decimal seconds from time
    lt_time_char=char(lt_time_filt);
    lt_time_char=lt_time_char(:,1:end-4);
    
    txt_data_array=[char(datestr(date_intrp_filt,'yyyymmdd')),repmat(',',size(date_intrp_filt,1),1),lt_time_char];
    
    [rad_data_m,rad_data_n]=size(rad_data_array);
    [anc_data_m,anc_data_n]=size(anc_data_array);
    format_str=['%s,%6.4f,%6.4f,%.1f,%.1g,%3.1f,%4.1f,%4.1f,%3.1f,%.1g,%4.2f,%4.2f,',repmat('%.3E,',1,137.*5),...
        repmat('%.3E,',1,5.*3),repmat('%.3E,',1,4),repmat('%.3E,',1,3),'%.3E\n'];

    %% Write SeaBASS file
    
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

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu



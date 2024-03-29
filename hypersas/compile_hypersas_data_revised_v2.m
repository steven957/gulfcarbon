% compile_hypersas_data_revised_v2.m - Program to read hypersas data from
%  sensor dat files and compile 
% 
% Syntax:  compile_hypersas_data_revised_v2.m
%
% Inputs:
%    1) Folder location with Level 2 sensor data files ('*.dat') 
%    2) File parameter spreadsheet file 'MMM2009_hypersas_file_list.xls'
%
% Outputs:
%    output - Mat file for each cruise day with structured variable 'all_hyper"
%    containing various fields 
%   
% Other m-files required: 
%  1) The function 'hyper_readdat' is required to read the hypersensor 
%     data files 
%  1) The function 'discrete_readdat' is required to read the discrete rsensor 
%     data files 
%  1) The function 'gps_readdat' is required to read the gps data files 
%  1) The function 'ths_readdat' is required to read the pitch and roll sensor 
%     data files 
%
% MAT-files required: none
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 13 Aug 2021

%------------- BEGIN CODE --------------%

clear all
close all
clc

%Input data file info
homepath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\';
%infile_name=[homepath,inpath,'jan2009_hypersas_file_list.xls'];
%infile_name=[homepath,inpath,'jul2009_hypersas_file_list.xls'];
%infile_name=[homepath,inpath,'nov2009_hypersas_file_list.xls'];
%infile_name=[homepath,inpath,'mar2010_hypersas_file_list.xls'];

%inpath='Jan2009\HyperSAS\HyperSAS_data_GC1\Output\';
%inpath='GC_hypersensor_data\GC3_hypersensors\HyperSAS_GC3\Output\';
%%inpath='GC_hypersensor_data\GC5_hyperSAS\Output\';
inpath='GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\';
infile_name=[homepath,inpath,'apr2009_hypersas_file_list.xls'];

%Read metadata from parameter spreadsheet
[~,hyper_type] = xlsread(infile_name,1,'E3:E821');
[hyper_run, ~] = xlsread(infile_name,1,'F3:F821');
[~,FName] = xlsread(infile_name,1,'C3:C821');
[~,hyper_sensor] = xlsread(infile_name,1,'D3:D821');

%Loop through different DOY for cruise and create compiled file for each
%  day
[maxrow,n]=size(hyper_sensor);
for file_day=110:120 %12:15  %70:79   %302:310  %201:210    %110:120
    clear all_hyper
    file_index=[];
    %Identify files corresponding to given DOY
    for file_nn=1:maxrow
            infile=char(FName(file_nn));
            if str2double(infile(6:8))==file_day
                file_index=[file_index;file_nn];
            end
    end

    %Initialize structured variable with file info
    all_hyper.hyper_type=hyper_type(file_index);
    all_hyper.hyper_sensor=hyper_sensor(file_index);
    all_hyper.hyper_run=hyper_run(file_index);
    all_hyper.hyper_name=FName(file_index);

    file_count=size(file_index,1);
    
    %Loop through each file, determine sensor, and extract data
    for k=1:file_count    
        infile=[homepath,inpath,char(FName(file_index(k)))];
            % Hyper sensors  
            if strncmp(all_hyper.hyper_type(k),'Hyper',5)==1
                [hyper_index,hyper_rad_dat,hyper_hdr,hyper_date,hyper_timer,hyper_frame,hyper_time,...
                    hyper_int_time,hyper_samp_delay]=hyper_readdat(infile);

                all_hyper.hyper_data{k}=struct('hyper_index',hyper_index,'hyper_rad_dat',hyper_rad_dat,...
                    'hyper_hdr',{hyper_hdr},'hyper_date',{hyper_date},'hyper_timer',hyper_timer,...
                        'hyper_frame',hyper_frame,'hyper_time',{hyper_time},'hyper_int_time',{hyper_int_time},...
                        'hyper_samp_delay',{hyper_samp_delay});
            end

            % Discrete irradiance sensor DI7130D
            if strncmp(all_hyper.hyper_sensor(k),'DI7130D',5)==1
                [discrete_index,discrete_rad_dat,discrete_hdr,discrete_date,discrete_timer,discrete_frame,...
                    discrete_time,discrete_samp_delay]=discrete_readdat(infile);

                    all_hyper.hyper_data{k}=struct('discrete_index',discrete_index,'discrete_rad_dat',discrete_rad_dat,...
                        'discrete_hdr',{discrete_hdr},'discrete_date',{discrete_date},...
                        'discrete_timer',discrete_timer,'discrete_frame',discrete_frame,...
                        'discrete_time',{discrete_time},'discrete_samp_delay',{discrete_samp_delay});
            end

            % Discrete radiance sensors DR7068e and DR7069d
            if ((strncmp(all_hyper.hyper_type(k),'Discrete',5)==1)&&(strncmp(hyper_sensor(k),'DI7130D',5)==0))
                [discrete_index,discrete_rad_dat,discrete_hdr,discrete_date,discrete_timer,discrete_frame,...
                    discrete_time,discrete_samp_delay]=discrete_readdat(infile);

                    all_hyper.hyper_data{k}=struct('discrete_index',discrete_index,'discrete_rad_dat',discrete_rad_dat,...
                        'discrete_hdr',{discrete_hdr},'discrete_date',{discrete_date},...
                        'discrete_timer',discrete_timer,'discrete_frame',discrete_frame,...
                        'discrete_time',{discrete_time},'discrete_samp_delay',{discrete_samp_delay});
            end

            % GPS data
            if strcmp(all_hyper.hyper_type(k),'GPS')==1
                [gps_dat, gps_hdr, gps_datenum, gps_hour, gps_min, gps_sec, gps_latdeg, gps_latmin, ...
                        gps_londeg, gps_lonmin, gps_course]=gps_readdat(infile);        

                all_hyper.hyper_data{k}=struct('gps_dat',{gps_dat},'gps_hdr',{gps_hdr},'gps_datenum',{gps_datenum},...
                        'gps_hour',{gps_hour},'gps_min',{gps_min},'gps_sec',{gps_sec},...
                        'gps_latdeg',{gps_latdeg},'gps_latmin',{gps_latmin},'gps_londeg',{gps_londeg},...
                        'gps_lonmin',{gps_lonmin},'gps_course',{gps_course});

            end

            if strcmp(all_hyper.hyper_type(k),'THS')==1
                [ths_index,ths_dat,ths_hdr,ths_roll,ths_pitch,ths_date,ths_timer,ths_frame,ths_time,ths_comp]=...
                    ths_readdat(infile);            

                all_hyper.hyper_data{k}=struct('ths_index',ths_index,'ths_dat',{ths_dat},'ths_hdr',{ths_hdr},...
                    'ths_roll',ths_roll,'ths_pitch',ths_pitch,'ths_date',{ths_date},'ths_timer',ths_timer,...
                        'ths_frame',ths_frame,'ths_time',{ths_time},'ths_comp',{ths_comp});
             end

        end

        disp('Completed, saving...');

        %save all_hyper_Apr2009_110.mat all_hyper;
        %eval(['save ',homepath,inpath,'all_hyper_Jan2009_',num2str(file_day),'.mat all_hyper']);
        %eval(['save ',homepath,inpath,'all_hyper_Apr2009_',num2str(file_day),'.mat all_hyper']);
        %eval(['save ',homepath,inpath,'all_hyper_Nov2009_',num2str(file_day),'.mat all_hyper']);
        %eval(['save ',homepath,inpath,'all_hyper_Mar2010_',num2str(file_day),'.mat all_hyper']);
    end

    disp({'Execution Complete'});

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu



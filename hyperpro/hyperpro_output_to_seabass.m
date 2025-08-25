% hyperpro_output_to_seabass.m
% 
% Program to filter hyperpro data output to seabass format files (*.sb)
%
% Syntax: hyperpro_output_to_seabass.m
%
% Inputs:
%    1) Processed and compiled HyperPro data ('*.mat') produced by
%     'hyperpro_matfile_compile.m' and 'hyperpro_process.m'
%
% Outputs:
%    1) '*.sb' seabass formatted files 
%   
% Other m-files required:  hyperpro_data_filter.m  - This filters hyperpro rrs data and calculates near surface kd and ed_0
%
% MAT-files required: 
%  1) Processed and compiled HyperPro data ('*.mat') produced by
%     'hyperpro_matfile_compile.m' and 'hyperpro_process.m'
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
%
% Last revision: 25 Aug 2025

%% ------------- BEGIN CODE --------------%

clc

% Load file locations and input parameters


icr = 3;  % For batch processing, commment this statement out

crz = ['GC',num2str(icr)];

% for icr = 2:5  % To loop through all cruises
    
    switch crz
        case 'GC2'
            inpath = 'C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_hyperPro_data\';
            cruise_txt = 'GC2';
    
            input_meta_file={'hyperpro_seabass_header_template_file.csv'};
            %Load metadata file
            disp(['Reading ',char(input_meta_file)]);
            
            fid1=fopen([inpath,'Spreadsheet\',char(input_meta_file)]);
                metafile_txt=textscan(fid1,...
                    '%s','delimiter','\n','TreatasEmpty','EOF');
            fclose(fid1);
            
            %Spreadsheet with input parameters
            infile_name=[inpath,'Spreadsheet\GC2_hyperpro_list.xlsx'];
        case 'GC3'
            inpath = 'C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC3_hypersensors\HyperPro_GC3\';
            cruise_txt = crz;
    
            input_meta_file={'hyperpro_seabass_header_template_file.csv'};
            %Load metadata file
            disp(['Reading ',char(input_meta_file)]);
            
            fid1=fopen([inpath,'Spreadsheet\',char(input_meta_file)]);
                metafile_txt=textscan(fid1,...
                    '%s','delimiter','\n','TreatasEmpty','EOF');
            fclose(fid1);
            
            %Spreadsheet with input parameters
            infile_name=[inpath,'Spreadsheet\GC3_hyperpro_list.xlsx'];
        case 'GC4'
            inpath ='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC4_hyper\HyperPro_data_GC4\';
            cruise_txt = 'GC4';  % Specify cruise for file loading and naming
    
            input_meta_file={'hyperpro_seabass_header_template_file.csv'};
            %Load metadata file
            disp(['Reading ',char(input_meta_file)]);
            
            fid1=fopen([inpath,'Spreadsheet\',char(input_meta_file)]);
                metafile_txt=textscan(fid1,...
                    '%s','delimiter','\n','TreatasEmpty','EOF');
            fclose(fid1);
            
            %Spreadsheet with input parameters
            infile_name=[inpath,'Spreadsheet\GC4_hyperpro_list.xlsx'];
       case 'GC5'
            inpath ='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC5_hyperPro\';
            cruise_txt = 'GC5';
    
            input_meta_file={'hyperpro_seabass_header_template_file.csv'};
            %Load metadata file
            disp(['Reading ',char(input_meta_file)]);
            
            fid1=fopen([inpath,'Spreadsheet\',char(input_meta_file)]);
                metafile_txt=textscan(fid1,...
                    '%s','delimiter','\n','TreatasEmpty','EOF');
            fclose(fid1);
            
            %Spreadsheet with input parameters
            infile_name=[inpath,'Spreadsheet\GC5_hyperpro_list.xlsx'];
    end
    
    % To do single station - need to modify code below as well
    % sta_txt = 'g4';
    
    % Flag to save output if desired
    save_output = false;  %true; %

    % Baseline correction flag (1 - yes, 0 - no)
    base_corr2 = 1;
    
    procolor1='#F5973A';
    procolor2='#0071bc';
    procolor3='#cd2026';
    procolor4='#4c2c92';
    procolor5='#2e8540';
    
    disp('Loading data...');
    % Load HyperPro data
    load([inpath,cruise_txt,'_hyperpro_processed.mat']);

    % Load metadata
    opts = detectImportOptions(infile_name);
    opts.VariableNamingRule='preserve'; 
    hyper_meta=readtable(infile_name,opts);

    disp('Loaded.')
    
    [~,ngrp]=size(pro_group_struct);
    
    % Initialize variables
    all_stations={};
    
    %% Begin loop to select station(s) and plot 
    
    for im=1:ngrp 
    
        % If group does not exist, advance
        if isempty(pro_group_struct{im})
            continue
        end

        profile_data=pro_group_struct{im};
        [n,pro_num]=size(profile_data.pro_date);
                  
        pro_indices = 1:pro_num;

        % Loop through profiles and output as SeaBASS files

        for ipro = 1:pro_num
            single_sta=profile_data.pro_station{ipro};  %Select station to plot
            
            % Count number of repeated profiles and add to station number
            target_string = single_sta;
            count = sum(strcmp(all_stations, target_string));
            
            if count>0
                sta_label=strrep(single_sta,single_sta,[single_sta,'_',num2str(count+1)]);
            else
                sta_label = single_sta;
            end
            all_stations=cat(2,all_stations,single_sta);
            
            %To do all stations, comment out the if continue statement; otherwise only
            % station given in single_sta will be processed
            % if ~contains(single_sta,sta_txt(1,:),'IgnoreCase',true)
            %     continue
            % else
                 sta_txt=char(single_sta);
            % end
            
            %% Get HyperPro station data and plot if desired

            sta_indx = find(strcmp(hyper_meta.Station,sta_txt));
            
            % Station metadata
            pro_bz=hyper_meta.bz(sta_indx(1));
            pro_ws=hyper_meta.Windspd(sta_indx(1));
            pro_cld=hyper_meta.CldCvr(sta_indx(1));
            pro_ssh=hyper_meta.SSH(sta_indx(1));

            pro_date=char(profile_data.pro_date(ipro),'yyyyMMdd');
            pro_time=char(profile_data.pro_date(ipro),'HH:mm:ss');
            pro_lat=profile_data.pro_lat(ipro);
            pro_lon=profile_data.pro_lon(ipro);
            pro_raw_file=profile_data.raw_file{ipro};

            if contains(sta_txt,'UWS')
                file_label=[char(profile_data.pro_date(ipro),'MMMuuuu'),'_',upper(sta_label)];
            else
                file_label=[char(profile_data.pro_date(ipro),'MMMuuuu'),'_St',upper(sta_label)];
            end

            % Profile data
            pres_interp=profile_data.pres_interp{ipro};
            pro_tilt=profile_data.pro_tilt{ipro};
            pro_vel=profile_data.pro_vel{ipro};
            kdz=profile_data.k_ed{ipro};
            edz=profile_data.ed_corr_es{ipro};

            lu=profile_data.lu_corr_es{ipro};
            es=profile_data.es_corr_interp{ipro};
            rrs=profile_data.rrs{ipro};
        
            ed_lambda=profile_data.ed_lambda{ipro};
        
            % Apply filter to hyperpro rrs if desired
    
            hyperpro_data_filter
            
        %{
            %% Plot spectra for QA/QC if desired
            for ipro=pro_indices  
                lambda_min = 375;  % 399.8
                lambda_max = 750;   % 700.2
                plt_lambda=ed_lambda{ipro}(ed_lambda{ipro}>lambda_min & ed_lambda{ipro}<lambda_max);
                plt_kd=k_ed{ipro}(ed_lambda{ipro}>lambda_min & ed_lambda{ipro}<lambda_max);
                plt_edz=edz{ipro}(:,ed_lambda{ipro}>lambda_min & ed_lambda{ipro}<lambda_max);
                plt_es=es{ipro}(:,ed_lambda{ipro}>lambda_min & ed_lambda{ipro}<lambda_max);
                plt_edz(plt_edz<0.1)=NaN;
            end
        %} 
    
            %% Output data to seabass format 
        
            %Create field headings for radiometric wavelengths
            wavel_str=num2str(ed_lambda,'%4.1f');
            Ed_fields=[repmat(',Ed',137,1),wavel_str];
            Lu_fields=[repmat(',Lu',137,1),wavel_str];
            Es_fields=[repmat(',Es',137,1),wavel_str];
            Rrs_fields=[repmat(',Rrs',137,1),wavel_str];
        
            % The statements below include the velocity for the profiler,
            %      but SeaBASS does not have a field for that
            % fields=['/fields=depth,tilt,vel',reshape(Ed_fields',1,[]),reshape(Lu_fields',1,[]),...
            %    reshape(Es_fields',1,[]),reshape(Rrs_fields',1,[])];
            % units=['/units=m,degrees,m/s',repmat(',uW/cm^2/nm',1,137),repmat(',uW/cm^2/nm/sr',1,137),...
            %     repmat(',uW/cm^2/nm',1,137),repmat(',1/sr',1,137)];
        
            fields=['/fields=depth,tilt',reshape(Ed_fields',1,[]),reshape(Lu_fields',1,[]),...
               reshape(Es_fields',1,[]),reshape(Rrs_fields',1,[])];
            units=['/units=m,degrees',repmat(',uW/cm^2/nm',1,137),repmat(',uW/cm^2/nm/sr',1,137),...
                repmat(',uW/cm^2/nm',1,137),repmat(',1/sr',1,137)];

            %rad_data_file_name=GulfCarbon_gulfcarbon2_hypersas_St#_2009MMDD_R1.sb
            rad_data_file_name=['GulfCarbon',num2str(icr),'_hyperpro_',file_label,'_R1.sb'];
                
            textstring=metafile_txt;
            [hdr_m]=size(char(textstring{1}),1);
            file_creation_date=char(datetime,'yyyyMMdd');
            
       
            %% Input metadata into file header
            for i=1:hdr_m
                if ~contains(char(textstring{1}(i)),'/north_latitude')~=1
                    textstring{1}(i)={['/north_latitude=',num2str(pro_lat),'[DEG]']};
                elseif isempty(strfind(char(textstring{1}(i)),'/south_latitude'))~=1 %#ok<*STREMP>
                    textstring{1}(i)={['/south_latitude=',num2str(pro_lat),'[DEG]']};
                elseif isempty(strfind(char(textstring{1}(i)),'/west_longitude'))~=1
                    textstring{1}(i)={['/west_longitude=',num2str(-pro_lon),'[DEG]']};
                elseif isempty(strfind(char(textstring{1}(i)),'/east_longitude'))~=1
                    textstring{1}(i)={['/east_longitude=',num2str(-pro_lon),'[DEG]']};
                elseif isempty(strfind(char(textstring{1}(i)),'/start_time'))~=1
                    textstring{1}(i)={['/start_time=',pro_time,'[GMT]']};
                elseif isempty(strfind(char(textstring{1}(i)),'/end_time'))~=1
                    textstring{1}(i)={['/end_time=',pro_time,'[GMT]']};
                elseif isempty(strfind(char(textstring{1}(i)),'/start_date'))~=1
                    textstring{1}(i)={['/start_date=',pro_date]};
                elseif isempty(strfind(char(textstring{1}(i)),'/end_date'))~=1
                    textstring{1}(i)={['/end_date=',pro_date]};
                elseif isempty(strfind(char(textstring{1}(i)),'/cruise'))~=1
                    textstring{1}(i)={['/cruise=gulfcarbon',num2str(icr)]};
                elseif isempty(strfind(char(textstring{1}(i)),'/station'))~=1
                    textstring{1}(i)={['/station=',sta_txt]};
                elseif isempty(strfind(char(textstring{1}(i)),'/cloud_percent'))~=1
                    textstring{1}(i)={['/cloud_percent=',num2str(pro_cld,'%i')]};
                elseif isempty(strfind(char(textstring{1}(i)),'/waveht'))~=1
                    textstring{1}(i)={['/wave_height=',num2str(pro_ssh,'%3.1f')]};
                elseif isempty(strfind(char(textstring{1}(i)),'/wind_speed'))~=1
                    textstring{1}(i)={['/wind_speed=',num2str(pro_ws,'%3.1f')]};
                elseif isempty(strfind(char(textstring{1}(i)),'/water_depth'))~=1
                    textstring{1}(i)={['/water_depth=',num2str(pro_bz,'%4.1f')]};
                elseif isempty(strfind(char(textstring{1}(i)),'/measurement_depth'))~=1
                    textstring{1}(i)={'/measurement_depth=na'};
                elseif isempty(strfind(char(textstring{1}(i)),'/fields'))~=1
                    textstring{1}(i)={fields};
                elseif isempty(strfind(char(textstring{1}(i)),'/units'))~=1
                    textstring{1}(i)={units};
                elseif isempty(strfind(char(textstring{1}(i)),'/original_file_name='))~=1
                    textstring{1}(i)={['/original_file_name=',char(pro_raw_file)]};
                elseif isempty(strfind(char(textstring{1}(i)),'/data_file_name='))~=1
                    textstring{1}(i)={['/data_file_name=',rad_data_file_name]};
                elseif isempty(strfind(char(textstring{1}(i)),'!file_creation_date='))~=1
                    textstring{1}(i)={['!file_creation_date=',file_creation_date]};
                end
            end
            
            %% Format output 
    
            % fields=['/fields=depth,tilt,',...
            %    reshape(Ed_fields',1,[]),reshape(Lu_fields',1,[]),...
            %    reshape(Es_fields',1,[]),reshape(Rrs_fields',1,[])];
            % units=['/units=m,degrees,...
            %     repmat(',uW/cm^2/nm/sr',1,137.*3),repmat(',uW/cm^2/nm',1,137),...
            %     repmat(',uW/cm^2/nm/sr',1,5.*3),repmat(',1/sr',1,4)];
            
            rad_data_array=[pres_interp,pro_tilt,edz,lu,es];
            
            % Filter data based on quality criteria
            rad_data_array=rad_data_array(~isnan(pres_interp) & pro_tilt<=5.5,:);
            
            % Eliminate rows containing NaN if desired
            % Create mask
            % contains_number = (isnan(rad_data_array));  % This creates a logical
            %    matrix of the same size as rad_data_array, where true indicates the presence
            %    of number_to_find and false indicates its absence. Identify rows
            %    containing the number.
            % rows_to_delete = any(contains_number, 2); % The any(..., 2) function checks each
            %    row (2 indicates dimension 2, i.e., all elements in a given row) contains_number. If any element
            %    in a row is true, the corresponding element in rows_to_delete will be true.
            % rad_data_array(rows_to_delete,:) = [];
            
            rrs_output = nan(size(rad_data_array,1),size(edz,2));
            rrs_output(1,:) = rrs2;
            rad_data_array = [rad_data_array,rrs_output];
            rad_data_array(isnan(rad_data_array))=-9999;  % Replace NaNs with -9999 for SeaBASS missing data convention
            [rad_data_m,rad_data_n]=size(rad_data_array);
            % format_str=['%.1f,%4.2f,%4.2f,',repmat('%.3E,',1,137.*3),...
            %     repmat('%.3E,',1,136),'%.3E\n'];
            format_str=['%.1f,%4.2f,',repmat('%.3E,',1,137.*3),...
                repmat('%.3E,',1,136),'%.3E\n'];
        
            %% Write SeaBASS file
            
            disp(['Writing ',rad_data_file_name]);
            fid2 = fopen([inpath,'To_SeaBASS\',rad_data_file_name], 'w');
                for ih=1:hdr_m
                    fprintf(fid2,'%s\n',char(textstring{1}(ih)));
                end
                for id=1:rad_data_m
                    fprintf(fid2,format_str,rad_data_array(id,:));
                end
                
            fclose(fid2);
        end    
    end
% end

fclose('all');

disp('Hyperpro output to seabass completed')

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu

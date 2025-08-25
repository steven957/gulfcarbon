% hyperpro_hdf_extract_attributes.m
% 
% Syntax: hyperpro_hdf_extract_attributes.m
%
% This program extracts profile attributes from *L2.hdf files generated using
% Prosoft and writes to a text file
%
% Inputs:
%    1) Folder location with HyperPro *L2.hdf files processed using Prosoft
%
% Outputs:
%    output - 'GC#_hyperpro_list.txt' file with file metadata for cruise
%   
% MAT-files required: 'None'
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
%
% Last revision: 25 Aug 2025 

%% ------------- BEGIN CODE --------------%% 

clearvars
clc

for ic = 1:5
    cruise_txt = ['GC',num2str(ic)];

    switch ic
        case 1
            inputfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\Gulf Carbon 1\Hyperpro_Data\';
        case 2
            inputfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_hyperPro_data\';
        case 3
            inputfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC3_hypersensors\HyperPro_GC3\';
        case 4
            inputfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC4_hyper\HyperPro_data_GC4\';
        case 5
            inputfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC5_hyperPro\';
    end
    
    file_lst=dir([inputfolder,'*L2.hdf']);
    [filen,m]=size(file_lst);
    
    lat_dec=cell(filen,1);
    lon_dec=cell(filen,1);
    date_info=cell(filen,1);
    date_val=cell(filen,1);
    pres_tare=cell(filen,1);
    station_id=cell(filen,1);
    raw_file=cell(filen,1);
    
    
    for nfile=1:filen % Processing files
        % Load data
        input_file=file_lst(nfile).name;
        input_file=char(input_file);
        labeltext=cellstr(input_file(1:6));
        filepath=[inputfolder,input_file];
    
        disp(['Reading ',filepath]);
    
        hdffile_info=hdfinfo(filepath);
    
        [u,sgroup_n]=size(hdffile_info.Vgroup);
    
        [n,m]=size(hdffile_info.Attributes);
    
        hdfmetadata_name=cell(m,1);
        hdfmetadata_value=cell(m,1);
        for i=1:m
            hdfmetadata_name{i}=hdffile_info.Attributes(i).Name;
            hdfmetadata_value{i}=hdffile_info.Attributes(i).Value;
        end
        
        lat_indx=find(strcmp('LATITUDE',hdfmetadata_name)); 
        lon_indx=find(strcmp('LONGITUDE',hdfmetadata_name));
        master_date_indx=find(strcmp('TIME-STAMP',hdfmetadata_name));
        tare_indx=find(strcmp('PRESSURE-TARE',hdfmetadata_name));
        stationid_indx=find(strcmp('STATION-ID',hdfmetadata_name));
        cal_indx=find(strcmp('CAL_FILE_NAMES',hdfmetadata_name)); 
        raw_file_indx=find(strcmp('RAW_FILE_NAME',hdfmetadata_name)); 
    
        lat=hdfmetadata_value{lat_indx};
        lon=hdfmetadata_value{lon_indx};
        if ~contains(lat,'"')
            lat_dec{nfile} = str2num(lat(1:2)) + str2num(lat(4:9))./60;
            lon_dec{nfile} = str2num(lon(1:2)) + str2num(lon(4:9))./60;
        else
            lat_dec{nfile} = str2num(lat(1:2)) + str2num(lat(4:5))./60  + str2num(lat(8:9))./3600;
            lon_dec{nfile} = str2num(lon(1:2)) + str2num(lon(4:5))./60  + str2num(lon(8:9))./3600;
        end
        date_info{nfile}=hdfmetadata_value{master_date_indx};
        date_val{nfile}=datetime([hdfmetadata_value{master_date_indx}(5:10),', 2009 ',hdfmetadata_value{master_date_indx}(12:19)]); 
        pres_tare{nfile}=hdfmetadata_value{tare_indx};
        station_id{nfile}=hdfmetadata_value{stationid_indx};
        raw_file{nfile}=hdfmetadata_value{raw_file_indx};
    end
    
    outtextfile=[cruise_txt,'_hyperpro_list.txt'];
    file_hdr_format='Raw File \t UTC Date \t Lat \t Lon \t PresTare \t Station \n';
    file_format_spec='%s \t %s \t %s \t %s \t %s \t %s \n';
    fileID=fopen([inputfolder,outtextfile],'w');
        fprintf(fileID,file_hdr_format);
        for dm=1:filen
            fprintf(fileID,file_format_spec,raw_file{dm},char(date_val{dm}),lat_dec{dm},lon_dec{dm},pres_tare{dm},station_id{dm});
        end
    fclose(fileID);
end

disp('Completed');


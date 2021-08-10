%read_hypersas_seabass.m - Function called by 'extract_hypersas_seabass.m' to 
% read SeaBASS data files
% 
% Syntax: [hyper_meta,hyper_run_date,hyper_run_time,hyper_dat]=read_hypersas_seabass(input_file)
%
% Inputs:
%    input - Folder location with SeaBASS format '*.sb' files
%
% Outputs:
%    output - Variables as defined below:
%     hyper_meta: metadata for HyperSAS observations including lat, lon,
%      date, UTC, station, data fields, and data units
%     hyper_run_date:  data of acquisition
%     hyper_run_time:  time of acquisition (UTC)
%     hyper_run_date:  HyperSAS data
%   
% Other m-files required: 
%  1) The SeaBASS data files were generated with 
%    'hypersas_rrs_Apr2009_output_to_seabass.m' 

% MAT-files required: none
%
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 26 July 2021

%------------- BEGIN CODE --------------%

function [hyper_meta,hyper_run_date,hyper_run_time,hyper_dat]=read_hypersas_seabass(input_file)

%Load hyper data file

disp(['Reading ',input_file]);

fid=fopen(input_file);

    if fid==0
        hyper_sensor=0;    
        hyper_hdr=0;
        return
    end

    textstring1='';

     while ~contains(textstring1,'end_header')==1
        textstring1=textscan(fid,'%s',1,'delimiter','\n','headerLines',0);
        textstring1=char(textstring1{1});
        index=strfind(textstring1,'=');

        if ~contains(textstring1,'/north_latitude')~=1
            hyper_lat=textstring1(index+1:end-5);
        end
        if ~contains(textstring1,'/east_longitude')~=1
            hyper_lon=textstring1(index+1:end-5);
        end
        if ~contains(textstring1,'/start_time')~=1
            hyper_time=textstring1(index+1:end-5);
        end
        %Reformat date from weird SeaBASS format
        if ~contains(textstring1,'/start_date')~=1
            Y=textstring1(index+1:index+4);
            M=textstring1(index+5:index+6);
            D=textstring1(index+7:index+8);
            hyper_date=[M,'/',D,'/',Y];
        end
        if ~contains(textstring1,'/measurement_depth')~=1
            hyper_dep=textstring1(index+1:end);
        end
        if ~contains(textstring1,'/station')~=1
            hyper_sta=textstring1(index+1:end);
        end
        if ~contains(textstring1,'/fields')~=1
            hyper_hdr=textscan(textstring1(index+1:end),'%s','delimiter',',');
        end
        if ~contains(textstring1,'/units')~=1
            hyper_units=textscan(textstring1(index+1:end),'%s','delimiter',',');
        end
    end

    textstring2=textscan(fid,...
        '%s','delimiter','\n','TreatasEmpty','EOF');

fclose(fid);

%Convert to string array
hyper_txt=char(textstring2{1});

hyper_run_date=hyper_txt(:,1:8);
hyper_run_time=hyper_txt(:,10:17);
hyper_dat=str2num(hyper_txt(:,19:end));

hyper_meta=struct('lat',hyper_lat,'lon',hyper_lon,'date',hyper_date,'UTC',hyper_time,...
    'depth',hyper_dep,'station',hyper_sta,'fields',hyper_hdr,'units',hyper_units);

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu
% read_ndbc_buoy_met_data_v2.m - Function to read in ndbc buoy met data file
% 
% Syntax: [met_hdr,met_units,met_datenum,met_wspd,met_wdir]=read_ndbc_buoy_met_data_v2(file_name)
%
% Inputs:
%    1) File name of NDBC buoy data file
%
% Outputs:
%    output - date, windspeed (m/s), wind direction (deg from N)
%   
% Other m-files required: 
%    None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 14 Aug 2021

%------------- BEGIN CODE --------------%

function [met_hdr,met_units,met_datenum,met_wspd,met_wdir]=read_ndbc_buoy_met_data_v2(file_name)

%Load header file
%hdr_file_name='C:\Steve\DATA\NASA\HyperSAS GEO-CAPE\Data\met_data_format_MetDat_s.C99.txt';

%Load data file
%file_name='C:\Steve\DATA\NASA\HyperSAS GEO-CAPE\Data\met_2015_MetDat_.C99.txt';

input_file=file_name;      

if(exist(input_file,'file')==0)
    disp(['File ',input_file,' does not exist!!!']);
    disp(['Continuing to next file...']);
    metdat=0;
    return
end

disp(['Reading ',input_file]);

fid1=fopen(input_file);
    met_hdr=textscan(fid1,'%s',1,'delimiter','\r');
    met_units=textscan(fid1,'%s',1,'delimiter','\r');
    format_str=['%u',repmat(' %u',1,5),repmat(' %f',1,5),'%u',repmat(' %f',1,6)];
    textdat=textscan(fid1,format_str,'delimiter','\r');  %read data 
fclose(fid1);

[m,~]=size(textdat{1});

if m<2
    met_hdr={};
    met_units={};
    met_rad_dat={};
    met_date={};
    met_time={};
    disp('No data found - exiting file read');
end

met_hdr=met_hdr{1};
met_units=met_units{1};

%met_dat=cell2mat(textdat);

met_year=double(textdat{1});
met_month=double(textdat{2});
met_day=double(textdat{3});
met_hh=double(textdat{4});
met_mm=double(textdat{5});
met_wspd=textdat{7};
met_wdir=textdat{6};

met_datenum=datenum(met_year,met_month,met_day,met_hh,met_mm,zeros(m,1));

disp('Completed');

return

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu

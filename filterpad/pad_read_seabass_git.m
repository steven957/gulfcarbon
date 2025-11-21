%pad_read_seabass_git.m
%
% Syntax: pad_extract_seabass_v2.m
%
% Inputs:
%    1) SeaBASS formatted filterpad files generated with 'pad_output_to_seabass_v2.m' 
%    2) Spreadsheet with station cluster and phytoplankton size fraction information
%    3) Spreadsheet with photosynthetron light spectrum 'ptron_light_spectrum.xlsx'
%
% Outputs:
%    1) Compiled SeaBASS filterpad data in 'all_pad.mat' including data and metadata
%   
% Other m-files required: None 
%
% MAT-files required: None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 21 November 2025

%% ------------- BEGIN CODE --------------%

function [pad_meta,pad_dat]=pad_read_seabass_git(input_file)

%Load pad data file

disp(['Reading ',input_file]);

fid=fopen(input_file);

if fid==0
    pad_sensor=0;    
    pad_hdr=0;
    return
end

textstring1='';

 while ~contains(textstring1,'end_header')==1
    textstring1=textscan(fid,'%s',1,'delimiter','\r','headerLines',0);
    textstring1=char(textstring1{1});
    index=strfind(textstring1,'=');

    if ~contains(textstring1,'/north_latitude')~=1
        pad_lat=textstring1(index+1:end-5);
    end
    if ~contains(textstring1,'/east_longitude')~=1
        pad_lon=textstring1(index+1:end-5);
    end
    if ~contains(textstring1,'/start_time')~=1
        pad_time=textstring1(index+1:end-5);
    end
    %Reformat date from weird SeaBASS format
    if ~contains(textstring1,'/start_date')~=1
        Y=textstring1(index+1:index+4);
        M=textstring1(index+5:index+6);
        D=textstring1(index+7:index+8);
        pad_date=[M,'/',D,'/',Y];
    end
    if ~contains(textstring1,'/measurement_depth')~=1
        pad_dep=textstring1(index+1:end);
    end
    if ~contains(textstring1,'/station')~=1
        pad_sta=textstring1(index+1:end);
    end
    if ~contains(textstring1,'/fields')~=1
        pad_hdr=textscan(textstring1(index+1:end),'%s','delimiter',',');
    end
    if ~contains(textstring1,'/units')~=1
        pad_units=textscan(textstring1(index+1:end),'%s','delimiter',',');
    end

end

textstring=textscan(fid,...
    '%s','delimiter','\r','TreatasEmpty','EOF');

fclose(fid);

%Convert to string array
pad_dat=char(textstring{1});
pad_dat=str2num(pad_dat);

pad_meta=struct('lat',pad_lat,'lon',pad_lon,'date',pad_date,'utc',pad_time,...
    'depth',pad_dep,'station',pad_sta,'fields',pad_hdr,'units',pad_units);

%End of pad_readdat.m


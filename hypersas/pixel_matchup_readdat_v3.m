% pixel_matchup_readdat_v3.m - Function to read satellite pixel matchup files
%
% Syntax: [sat_hdr, sat_dat, hyper_rrs_hdr, rrs_time, rrs_lat, rrs_lon, hyper_rrs_dat]=pixel_matchup_readdat_v3(file_name)
%
% Inputs:
%    1) Input file name for satellite matchup file ('*.txt') produced by 
%      ESA SNAP software
%
% Outputs:
%    output - Variables from satellite retrievals 
%   
% MAT-files required: none
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 6 Sep 2021

%% ------------- BEGIN CODE --------------%% 

function [sat_hdr, sat_dat, hyper_rrs_hdr, rrs_time, rrs_lat, rrs_lon, hyper_rrs_dat]=pixel_matchup_readdat_v3(file_name)

    %Load data file

    input_file=file_name;      

    if exist(input_file,'file')==0
        disp(['File ',input_file,' does not exist!!!']);
        disp(['Continuing to next file...']);
        hyperdat=0;
        return
    end

    disp(['Reading ',input_file]);
    %hdr_str=['%s',repmat(' %s',1,183),' %*[^\n]'];  %For April 2009
    hdr_str=['%s',repmat(' %s',1,534),' %*[^\n]'];
    format_str=['%u %s',repmat(' %f',1,139),' %s %u %u %s %f %f %f %f %s %s',repmat(' %f',1,384),'%*[^\n]']; 
    %format_str=['%u %s %s',repmat(' %f',1,8),' %s %d %u %s %f %f %f %f %s %s ',repmat(' %f',1,273),' %*[^\n]']; 
    %[repmat('%n ',1,146),'%8s %12s'];

    fid=fopen(input_file);
        null_str=textscan(fid,'%s',5,'delimiter','\r');  %eliminate header lines
        hdrdat=textscan(fid,hdr_str,1,'delimiter','\t','MultipleDelimsAsOne',0);   %read header line
        %null_str=textscan(fid,'%s',2,'delimiter','\r');  %eliminate headers

        %read data 
        textdat=textscan(fid,format_str,'delimiter','\t','MultipleDelimsAsOne',0);
    fclose(fid);
    
    [~,hdrn]=size(hdrdat);
    sat_hdr={};
    for hdri=1:hdrn
        sat_hdr=[sat_hdr,hdrdat{hdri}];
    end

    sat_dat=textdat;
    
    hyper_rrs_hdr=sat_hdr(5:141);
    hyper_rrs_dat=textdat(5:141);
    hyper_rrs_dat=cell2mat(hyper_rrs_dat);
     
    time_index=2; %strncmp('Time',sat_hdr,4)==1;
    lat_index=3; %strncmp('Lat',sat_hdr,8)==1;
    lon_index=4; %strncmp('Lon',sat_hdr,8)==1;
    
    rrs_time=textdat{time_index};
    rrs_lat=textdat{lat_index};
    rrs_lon=textdat{lon_index};

    disp('Completed Readat Function');
    
return

% ------------- END OF CODE --------------
% Please send suggestions for improvement
%  to Steven Lohrenz at this email address: slohrenz@umassd.edu   

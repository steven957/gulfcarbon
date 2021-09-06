% pixel_matchup_master_v2.m - Script to input satellite pixel matchup data
%
% Syntax: pixel_matchup_master_v2.m
%
% Inputs:
%    1) Folder with matchup files for satellite matchup file ('*.txt') produced by 
%      ESA SNAP software
%
% Outputs:
%    output - Variables from satellite retrievals 
%   
% Other m-files required: 
%  1) pixel_matchup_readdat_v3.m - function to read data files
% 
% MAT-files required: none
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 6 Sep 2021

%% ------------- BEGIN CODE --------------%% 
clc
clear all

%Input data file info
inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\Satellite\';
%inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC3_hypersensors\HyperPro_GC3\Output\Satellite\';

FName=dir([inpath,'rrs_matchup_Apr2009*_measurements.txt']);
%FName=dir([inpath,'jul2009_gc_product_matchup.txt']);

[maxrow,n]=size(FName);

%Setup structured variable for input
all_sat=struct('lat',{},'lon',{},'time',{},'run_date',{},'stn',{},'sat_fields',{},'sat_data',{},...
    'sat_rrs',{},'hyper_rrs',{},'fname',{},'sat_lambda',{},'hyper_lambda',{});

%Loop through each file in folder and input data
for k=1:maxrow
    filename=char(FName(k).name);
	infile=[inpath,filename];
	
    %Get run data and station from filenamee
    strtok=textscan(filename,'%s','Delimiter','_');
    strtok=strtok{1};
    
    run_date=strtok{3};
    stn=strtok{5};
    
    %------------------sat readdat function------------------------

    [sat_hdr, sat_dat, hyper_rrs_hdr, rrs_time, rrs_lat, rrs_lon, hyper_rrs_dat]=pixel_matchup_readdat_v3(infile);  %Call subroutine to read in data
    
    %-------------End sat readdat function------------------
    
    if isempty(sat_dat)
        continue
    end

    all_sat(k).lat=rrs_lat;
    all_sat(k).lon=rrs_lon;
    all_sat(k).time=rrs_time;
    all_sat(k).run_date=run_date;
    all_sat(k).stn=stn;
    all_sat(k).sat_fields=sat_hdr;
    
    sat_indx1=contains(sat_hdr,'Rrs_');
    sat_indx2=contains(sat_hdr,'_mean');
    sat_rrs_indx=find(contains(sat_hdr,'Rrs_') & contains(sat_hdr,'_mean'));
        
    sat_lambda=str2double(strrep(strrep(sat_hdr(sat_rrs_indx),'Rrs_',''),'_mean',''));    %Remove prefix and suffix and convert to double
    all_sat(k).sat_lambda=sat_lambda;

    hyper_lambda=str2double(hyper_rrs_hdr);
    all_sat(k).hyper_lambda=hyper_lambda;
    
    all_sat(k).fname=filename;

    %Data
    sat_rrs_cell=sat_dat(:,sat_rrs_indx);
    sat_rrs=cell2mat(sat_rrs_cell);
    all_sat(k).sat_rrs=sat_rrs;
    all_sat(k).sat_data=sat_dat;
    all_sat(k).hyper_rrs=hyper_rrs_dat;
    
end

all_sat_GC2=all_sat;

%save([inpath,'all_sat_GC2_v2.mat'],'all_sat_GC2');
%save([inpath,'all_sat_GC2_v2.mat'],'all_sat_GC2');

msgbox({'Execution Complete'});

% ------------- END OF CODE --------------
% Please send suggestions for improvement
%  to Steven Lohrenz at this email address: slohrenz@umassd.edu   




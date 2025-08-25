% hyperpro_extract_seabass - Program to read hyperpro data from SeaBASS
% formatted data files and compiles into a single .mat file 
% 
% Syntax:  hyperpro_extract_seabass.m
%
% Inputs:
%    input - Folder location with '*.sb' files 
%
% Outputs:
%    output - Single mat file (GC#_all_hyper_seabass.mat) with structured variable 'all_hyper"
%    containing various fields
%   
% Other m-files required: 
%  1) The function 'hyperpro_read_seabass.m' is required
%    to read the SeaBASS formatted data
%  2) The SeaBASS data files were generated with 
%    'hyperpro_output_to_seabass.m' 
%
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 24 Jul 2024

%% ------------- BEGIN CODE --------------%

clc

% For batch processing, commment statement out below
% icr = 5;

cruise = ['GC',num2str(icr)];

%Input data file info

switch cruise
    case 'GC2'
        inpath='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_hyperPro_data\To_SeaBASS\';
        cruise_name = 'Apr2009';
    case 'GC3'
        inpath='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC3_hypersensors\Hyperpro_GC3\To_SeaBASS\';
        cruise_name = 'Jul2009';    
    case 'GC4'
        inpath='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC4_hyper\HyperPro_data_GC4\To_SeaBASS\';
        cruise_name = 'Nov2009';
    case 'GC5'
        inpath='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC5_hyperPro\To_SeaBASS\';
        cruise_name = 'Mar2010';
end

FName=dir([inpath,'*.sb']);

[maxrow,~]=size(FName);

all_hyper=struct('lat',{},'lon',{},'utc_datetime',{},'station',{},...
    'name',{},'fields',{},'units',{},'data',{});
for k=1:maxrow
	infile=[inpath,char(FName(k).name)];
	
        [hyper_meta,hyper_utc_datetime,hyper_dat]=hyperpro_read_seabass(infile);  %Call subroutine to read in data
        if hyper_dat==0
            continue
        end
        
        %Metadata
        all_hyper(1).lat(k,1)={hyper_meta.lat};
        all_hyper(1).lon(k,1)={hyper_meta.lon};
        all_hyper(1).utc_datetime(k,1)={hyper_utc_datetime};
        all_hyper(1).station(k,1)={hyper_meta.station};
        all_hyper(1).fields(k,1)={hyper_meta.fields};
        all_hyper(1).units(k,1)={hyper_meta.units};
        all_hyper(1).name(k,1)={FName(k).name};
        
        %Data
        all_hyper(1).data(k,1)={hyper_dat};
end
        
	save([inpath,'all_hyper_seabass_',cruise_name,'.mat'],'all_hyper');

    disp('Seabass extraction complete');

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu



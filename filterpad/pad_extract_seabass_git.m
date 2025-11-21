% pad_extract_seabass_git.m - Program to read pad from SeaBASS formatted data files
%
% Syntax: pad_extract_seabass_git.m
%
% Inputs:
%    1) SeaBASS formatted filterpad files generated with 'pad_output_to_seabass_v2.m' 
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

clearvars
clc
close all

%Input data file info

inpath='YOUR FOLDER NAME WITH SEABASS FILES';

FName=dir([inpath,'*.sb']);

[maxrow,~]=size(FName);

all_pad=struct('lat',{},'lon',{},'date',{},'utc',{},'depth',{},'station',{},...
    'name',{},'fields',{},'data',{},'units',{});
for k=1:maxrow
	infile=[inpath,char(FName(k).name)];
	
        [pad_meta,pad_dat]=pad_read_seabass_git(infile);  %Call subroutine to read in data
        if pad_dat==0
            continue
        end
        
        %Metadata
        all_pad(1).lat(k,1)={pad_meta.lat};
        all_pad(1).lon(k,1)={pad_meta.lon};
        all_pad(1).date(k,1)={pad_meta.date};
        all_pad(1).utc(k,1)={pad_meta.utc};
        all_pad(1).station(k,1)={pad_meta.station};
        all_pad(1).depth(k,1)={pad_meta.depth};
        all_pad(1).fields(k,1)={pad_meta.fields};
        all_pad(1).units(k,1)={pad_meta.units};
        all_pad(1).name(k,1)={FName(k).name};
        
        %Data
        all_pad(1).data(k,1)={pad_dat};
end
        
   
	save([inpath,'all_pad.mat'],'all_pad');

    disp('Execution Complete');





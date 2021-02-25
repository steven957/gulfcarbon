%Program to read hyper from SeaBASS formatted data files
%This program uses function 'read_hypersas_seabass.m' to read SeaBASS data
% files generated with 'hypersas_rrs_Apr2009_output_to_seabass.m' and 
% compiles results in a single .mat file

clc
clear all

%Input data file info

inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\To_SeaBASS\';

FName=dir([inpath,'*.sb']);

[maxrow,~]=size(FName);

all_hyper=struct('lat',{},'lon',{},'start_date',{},'start_utc',{},'depth',{},'station',{},...
    'name',{},'fields',{},'units',{},'run_date',{},'run_time',{},'data',{});
for k=1:maxrow
	infile=[inpath,char(FName(k).name)];
	
        [hyper_meta,hyper_run_date,hyper_run_time,hyper_dat]=read_hypersas_seabass(infile);  %Call subroutine to read in data
        if hyper_dat==0
            continue
        end
        
        %Metadata
        all_hyper(1).lat(k,1)={hyper_meta.lat};
        all_hyper(1).lon(k,1)={hyper_meta.lon};
        all_hyper(1).start_date(k,1)={hyper_meta.date};
        all_hyper(1).start_utc(k,1)={hyper_meta.UTC};
        all_hyper(1).station(k,1)={hyper_meta.station};
        all_hyper(1).depth(k,1)={hyper_meta.depth};
        all_hyper(1).fields(k,1)={hyper_meta.fields};
        all_hyper(1).units(k,1)={hyper_meta.units};
        all_hyper(1).name(k,1)={FName(k).name};
        
        %Run date and run time
        all_hyper(1).run_date={hyper_run_date};
        all_hyper(1).run_time={hyper_run_time};
        
        %Data
        all_hyper(1).data(k,1)={hyper_dat};
end
        
   
	save([inpath,'all_hyper_April2009.mat'],'all_hyper');

    msgbox({'Execution Complete'});





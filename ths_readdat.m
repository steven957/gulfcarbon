
function [ths_index,ths_dat,ths_hdr,ths_roll,ths_pitch,ths_date,ths_timer,ths_frame,ths_time,ths_comp]=...
    ths_readdat(file_name)

%Load data file

input_file=file_name;      

if(exist(input_file)==0)
    disp(['File ',input_file,' does not exist!!!']);
    disp(['Continuing to next file...']);
    thsdat=0;
    return
end

disp(['Reading ',input_file]);
hdr_str='%s';
format_str=[repmat('%n ',1,10),'%8s %12s'];

fid=fopen(input_file);
    null_str=textscan(fid,'%s',3,'delimiter','\r');  %eliminate headers
    hdrdat=textscan(fid,'%s',1,'delimiter','\r');
    null_str=textscan(fid,'%s',1,'delimiter','\r');  %eliminate headers
    textdat=textscan(fid,format_str,'delimiter','\r');  %read data 
fclose(fid);

ths_hdr=parse_string(char(hdrdat{1}));
ths_dat=textdat;

%Input data into matrix format
%[m,n]=size(textdat{1});

    ths_index=textdat{1};             %Convert to numeric
    ths_roll=textdat{4};    
    ths_pitch=textdat{5};    
    ths_date=textdat{11};    
    ths_time=textdat{12};
    ths_frame=textdat{2};
    ths_timer=textdat{3};
    ths_comp=textdat{10};

disp('Completed');

return

%End of ths_readdat_Oct2005.m


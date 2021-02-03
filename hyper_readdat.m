
function [hyper_index,hyper_rad_dat,hyper_hdr,hyper_date,hyper_timer,hyper_frame,hyper_time,hyper_int_time,hyper_samp_delay]=...
    hyper_readdat(file_name)

%Load data file

input_file=file_name;      

if exist(input_file,'file')~=2
    disp(['File ',input_file,' does not exist!!!']);
    disp(['Continuing to next file...']);
    hyperdat=0;
    return
end

disp(['Reading ',input_file]);
format_str=['%d %d ',repmat('%n ',1,139),'%d %d %n %d %n %d %8s %12s'];

fid=fopen(input_file);
    null_str=textscan(fid,'%s',3,'delimiter','\r');  %eliminate headers
    hdrdat=textscan(fid,'%s',1,'delimiter','\r');
    null_str=textscan(fid,'%s',1,'delimiter','\r');  %eliminate headers
    textdat=textscan(fid,format_str,'delimiter','\r');  %read header data 
fclose(fid);

hyper_hdr=parse_string(char(hdrdat{1}));

%Input data into matrix format
[m,n]=size(textdat{1});

    hyper_index=textdat{1};             %Convert to numeric
  for j=5:141
      hyper_rad_dat(:,j-4)=textdat{j};
  end
    hyper_int_time=textdat{3};
    hyper_samp_delay=textdat{4};
    hyper_date=textdat{148};    
    hyper_time=textdat{149};
    hyper_frame=textdat{145};
    hyper_timer=textdat{146};

disp('Completed');

return

%End of hyper_readdat.m


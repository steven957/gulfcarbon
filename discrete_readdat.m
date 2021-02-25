
function [discrete_index,discrete_rad_dat,discrete_hdr,discrete_date,discrete_timer,...
    discrete_frame,discrete_time,discrete_samp_delay]=discrete_readdat(file_name)

%Load data file

input_file=file_name;      

if exist(input_file,'file')~=2
    disp(['File ',input_file,' does not exist!!!']);
    disp(['Continuing to next file...']);
    discretedat=0;
    return
end

disp(['Reading ',input_file]);
%hdr_str='%s';

%Need to account for different file format for discrete ES sensor
if contains(input_file,'DI7130')==1
    format_str=[repmat('%n ',1,13),'%8s %12s'];
else
    format_str=[repmat('%n ',1,14),'%8s %12s'];
end

fid=fopen(input_file);
    null_str=textscan(fid,'%s',3,'delimiter','\r');  %eliminate headers
    hdrdat=textscan(fid,'%s',1,'delimiter','\r');
    null_str=textscan(fid,'%s',1,'delimiter','\r');  %eliminate headers
    textdat=textscan(fid,format_str,'delimiter','\r');  %read data 
fclose(fid);

discrete_hdr=parse_string(char(hdrdat{1}));

%Input data into matrix format
%[~,n]=size(textdat{1});

    discrete_index=textdat{1};             %Convert to numeric
  for j=5:9
      discrete_rad_dat(:,j-4)=textdat{j};
  end
 
if contains(input_file,'DI7130')==1
    discrete_date=textdat{14};
    discrete_time=textdat{15};
    discrete_frame=textdat{12};
    discrete_timer=textdat{3};
    discrete_samp_delay=textdat{4};
else
    discrete_date=textdat{15};
    discrete_time=textdat{16};
    discrete_frame=textdat{13};
    discrete_timer=textdat{3};
    discrete_samp_delay=textdat{4};
end

disp('Completed');

return

%End of discrete_readdat_Oct2005.m


%gps_readdat.m - Function called by 'compile_hypersas_data_revised_v2.m' to 
% read gps data files
% 
% Syntax: [gps_dat, gps_hdr, gps_datenum, gps_hour, gps_min, gps_sec, gps_latdeg, gps_latmin,...
%     gps_londeg, gps_lonmin, gps_course]=gps_readdat(file_name)
%
% Inputs:
%    1) File name of Level 2 sensor data files ('*.dat') 
%
% Outputs:
%    output - Variables used in processing gps data
%   
% Other m-files required: 
%    None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 13 Aug 2021

%------------- BEGIN CODE --------------%

function [gps_dat, gps_hdr, gps_datenum, gps_hour, gps_min, gps_sec, gps_latdeg, gps_latmin, ...
    gps_londeg, gps_lonmin, gps_course]=gps_readdat(file_name)

%Load data file

input_file=file_name;      

if(exist(input_file)==0)
    disp(['File ',input_file,' does not exist!!!']);
    disp(['Continuing to next file...']);
    hyperdat=0;
    return
end

disp(['Reading ',input_file]);
hdr_str='%s';
format_str='%s %s %s %s %s %s %s %n %n %s %n %s %n %s %s' ;
%[repmat('%n ',1,146),'%8s %12s'];

fid=fopen(input_file);
    null_str=textscan(fid,'%s',3,'delimiter','\r');  %eliminate headers
    hdrdat=textscan(fid,'%s',1,'delimiter','\r');
    null_str=textscan(fid,'%s',2,'delimiter','\r');  %eliminate headers
    %textdat=textscan(fid,format_str,'delimiter','\r');  %read data 
    textdat=textscan(fid,format_str,'delimiter','\t','MultipleDelimsAsOne',1);
      %read data 
fclose(fid);

gps_hdr=parse_string(char(hdrdat{1}));
gps_dat=textdat;

%Input data into matrix format
    strmat=find(strcmp('DATE',gps_hdr)==1);
    gps_date=textdat{strmat(1)};
    gps_datenum=datenum(gps_date,'dd/mm/yy');
    
    gps_time=textdat{strcmp('UTCPOS',gps_hdr)==1};
    indx=strfind(gps_time,'.');
    [m, n]=size(indx);
    
   gps_hour=zeros(m,n);
   gps_min=zeros(m,n);
   gps_sec=zeros(m,n);
   
   for j=1:m
        if strcmp(char(gps_time{j}(1,3)),':')==1
            gps_hour(j)=str2double(char(gps_time{j}(:,1:2)));
            gps_min(j)=str2double(char(gps_time{j}(:,4:5)));
            gps_sec(j)=str2double(char(gps_time{j}(:,7:end)));
        else
            gps_hour(j)=str2double(char(gps_time{j}(:,1)));
            gps_min(j)=str2double(char(gps_time{j}(:,3:4)));
            gps_sec(j)=str2double(char(gps_time{j}(:,5:end)));
        end
   end
         
     gps_latdeg=zeros(m,n);
     gps_latmin=zeros(m,n);
     gps_londeg=zeros(m,n);
     gps_lonmin=zeros(m,n);
     
     gps_lat=textdat{strcmp('LATPOS',gps_hdr)==1};
     gps_lon=textdat{strcmp('LONPOS',gps_hdr)==1};
     
    for j=1:m
        gps_latdeg(j)=str2double(char(gps_lat{j}(:,1:2)));
        gps_latmin(j)=str2double(char(gps_lat{j}(:,3:end-1)));
        gps_londeg(j)=str2double(char(gps_lon{j}(:,1:2)));
        gps_lonmin(j)=str2double(char(gps_lon{j}(:,3:end-1)));
    end
  
    gps_course=textdat{strcmp('COURSE(TRUE)',gps_hdr)==1};

disp('Completed');

return

%End of gps_readdat.m

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu


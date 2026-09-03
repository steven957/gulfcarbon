% read_pace_kd490_PAR_aph442_netcdf_single_file.m
% Syntax:  read_pace_kd490_PAR_aph442_netcdf_single_file
%
% Script reads PACE netcdf file and extracts kd490, aph442, and par
%
% Inputs:
%    1) Directory locations for imagery and coast files
%    2) PACE level 2 IOP nc file
%
% Outputs:
%    1) Variables with kd490, aph442, and par
%   
% Other m-files required: None 
%
% MAT-files required: 
%    1) None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 2 Mar 2026

%% ------------- BEGIN CODE --------------%

% clearvars

% Read variables from L3 PACE or MODIS netcdf files

switch sensor_type
    case 'PACE'
        aph_target_wv = 440;  %m-1
    case 'MODIS'
        aph_target_wv = 443;  %m-1
end

varnames = {finfo.Variables.Name};
vardata = cell(length(varnames),1);

for ivar = 1:length(varnames)
    vardata{ivar} = ncread(filepath,['/',varnames{ivar}]);   
end

% Extract variables for wavelength of interest
switch iop
    case 'kd'
        var_indx = find(contains(varnames,'Kd_'));
        wv_kd = str2num(char(strrep(varnames(var_indx),'Kd_','')));
        kd = single(cat(3, vardata{var_indx}));
        kd490 = single(vardata{wv_kd==490});
        kd490(kd490<0 | kd490>20)=nan;
    case 'aph'
        var_indx = find(contains(varnames,'aph_'));
        scrap = strrep(varnames(var_indx),'aph_','');
        wv_aph = str2num(char(strrep(scrap,'_qaa','')));
        aph = single(cat(3, vardata{var_indx}));
        aph440 = single(vardata{wv_aph==aph_target_wv});
        aph440(aph440<0 | aph440>5)=0;
    case 'chlor_par'
        var_indx1 = find(contains(varnames,'ipar'));  % micromol m-2 s-1
        par = single(vardata{var_indx1});
        par(par<0 | par>3000)=0;
        % Handle case for chlorophyll concentration if needed
        var_indx2 = find(contains(varnames,'chlor_a'));
        chl = single(vardata{var_indx2});
        chl(chl<0 | chl>200)=nan;  % Adjust threshold as necessary
    case 'micro'
        var_indx = find(contains(varnames,'microplankton_uitz'));  % micromol m-2 s-1
        fmicro = single(vardata{var_indx});
        fmicro(fmicro<0 | fmicro>1)=0;
end

% Get image datetime and spatial location information
dt_start=ncreadatt(filepath,"/","time_coverage_start");
dt_end=ncreadatt(filepath,"/","time_coverage_end");
lat_south = ncreadatt(filepath,"/","southernmost_latitude");
lat_north = ncreadatt(filepath,"/","northernmost_latitude");
lon_west = ncreadatt(filepath,"/","westernmost_longitude");
lon_east = ncreadatt(filepath,"/","easternmost_longitude");

% Convert to datetime format
utc_start = char(datetime(dt_start,'Format','HH:mm:ss.SSS','InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS''Z'));
utc_end = char(datetime(dt_end,'Format','HH:mm:ss.SSS','InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS''Z'));
dt_start = datetime(dt_start,'Format','dd-MMM-uuuu HH:mm:ss.SSS','InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS''Z');
dt_end = datetime(dt_end,'Format','dd-MMM-uuuu HH:mm:ss.SSS','InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS''Z');

dt = char(dt_start + (dt_end - dt_start)./2,'dd-MMM-uuuu');
utc = char(datetime(utc_start) + (datetime(utc_end) - datetime(utc_start))./2,'HH:mm:ss.SSS');
dtm = datetime([dt,' ',utc],'InputFormat','dd-MMM-uuuu HH:mm:ss.SSS');

dt = char(dtm,'dd-MMM-uuuu HH:mm:ss.SSS');

% Get image lat and lon
lat=ncread(filepath,'lat');
lon=ncread(filepath,'lon');

disp('Completed reading satellite data...');




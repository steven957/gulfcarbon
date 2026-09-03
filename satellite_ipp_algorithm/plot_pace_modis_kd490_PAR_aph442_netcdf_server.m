% plot_pace_kd490_PAR_aph440_netcdf_server.m
% Syntax:  plot_pace_kd490_PAR_aph440_netcdf_server
%
% Script reads PACE netcdf file and extracts and plots Kd490 and other
% products: this version reads in multiple wavelengths for aph and Kd
%
% Inputs:
%    1) Directory locations for imagery and coast files
%    2) PACE level 2 or 3 IOP nc file (current version works with Level 3,
%    would require editing for Level 2)
%    3) High resolution coastal shape file ('GSHHS_f_L1.shp')
%
% Outputs:
%    1) Color map of optical variables (aph, kd)
%    2) Spectral plot of optical variables (aph, kd)
%    3) Spreadsheets with compiled mean kd spectra for each cluster and cruise, 
%       'kd_cluster_vector_means_PACE_monthly_composites.xlsx' or aph shape 
%        vectors for each size class and cruise, 'aph_size_class_means_PACE_monthly_composites.xlsx'
%   
% Other m-files required: None calhttps://umassd.zoom.us/j/97551651457?pwd=AF0Gp6oWyKUGVopirM8buiakyKbIJq.1c
%
% MAT-files required: 
%    1) None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 10 Aug 2026

%% ------------- BEGIN CODE --------------%

clc
clearvars
close all

% Processing level flag
level_flag = 3;
sensor_type = 'PACE';
hbar_flag = false;

%Input image data and plot

switch sensor_type
    case 'PACE'
        % Laptop folders
        % inputfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Satellite Imagery\PACE\';
        % ancfolder='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Satellite Imagery\PACE\Ancillary\';  %Enter the folder where your coastal shape file can be found
        % codepath='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\Steve\DATA\NSF\GulfCarbon\Matlab\Server Scripts\';
    
        % Server folders
        inputfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/pace/';
        ancfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/pace/Ancillary/';  %Enter the folder where your coastal shape file can be found
        codepath='~/ocean_color/matlab/';
    case 'MODIS'
        % Server folders
        inputfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/modis/';
        ancfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/pace/Ancillary/';  %Enter the folder where your coastal shape file can be found
        codepath='~/ocean_color/matlab/';
end

var_area_mean_all_cruise = table();  % Compiled mean spectra for all cruises (for averaging purposes)
monthrange = [4,7,11,3]; % 1:12;

for mnthn = monthrange

    % List of input files for IOP and PAR
    switch sensor_type
        case 'PACE'
            % prod_flag = {'kd','aph','chlor_par'};
            % iop_list = {'kd','aph','chlor_par'};
            iop_list = {'kd'};
            prod_flag = {'kd'};
        case 'MODIS'
            prod_flag = {'OC'}; % {'OC','aph'};calhttps://umassd.zoom.us/j/97551651457?pwd=AF0Gp6oWyKUGVopirM8buiakyKbIJq.1c
            endall_cruise_aph_size_class_tab


    nprod = length(prod_flag);

    for iprod = 1:nprod
        prod_type = prod_flag{iprod};
        sensor_type = sensor_type; % 'MODIS'; % 'PACE'

        % List of input files for IOP and PAR
        switch sensor_type
            case 'PACE'
                if strcmp(iop_list,'micro')
                    prod_type = 'chlor_par';
                end
                datafolder=[inputfolder,'L3_',prod_type,'/'];
                file_lst=dir([datafolder,'PACE_OCI.2025',sprintf('%02d', mnthn),'*_',prod_type,'.L3.nc']);
                iop_list = {prod_type};
                % Load data
                fname=char(file_lst.name);

                labeltext=fname(1:12);
                filepath=[datafolder,fname];
                finfo = ncinfo(filepath);

                disp(['Reading image data from ',fname]);

                for ivar = 1:length(iop_list)
                    iop = iop_list{ivar};
                    run([codepath,'read_l3_kd490_PAR_aph440_netcdf_single_file']);
                end

            case 'MODIS'
                if strcmp(prod_type,'OC')
                    iop_list = {'kd','chlor_par','aph'};
                elseif strcmp(prod_type,'aph')
                    iop_list = {'aph'};
                end
                datafolder=inputfolder;
                file_lst=dir([datafolder,'A2025',sprintf('%02d', mnthn),'*_',prod_type,'.L3.nc']);

                % Load data
                fname=char(file_lst.name);

                labeltext=fname(1:9);
                filepath=[datafolder,fname];
                finfo = ncinfo(filepath);

                disp(['Reading image data from ',fname]);
                all_cruise_aph_size_class_tab

                for ivar = 1:length(iop_list)
                    iop = iop_list{ivar};
                    run([codepath,'read_l3_kd490_PAR_aph440_netcdf_single_file']);
                end
        end
    end

    mnth = month(dtm,'shortname');
    switch mnth{1}
        case 'Jan'
            cruise_name = 'GC1'; % For station cluster averaging
        case 'Apr'
            cruise_name = 'GC2'; % For station cluster averaging
        case 'Jul'
            cruise_name = 'GC3'; % For station cluster averaging
        case 'Nov'
            cruise_name = 'GC4'; % For station cluster averaging
        case 'Mar'
            cruise_name = 'GC5'; % For station cluster averaging
    end

    yr = year(dtm); 
    dt_label = [char(mnth),' ',num2str(yr)];
    
    lat=ncread(filepath,'lat');
    lon=ncread(filepath,'lon');

    % Get image dimensions
    [imagez,imagem,imagen] = size(vardata{1});

    % Reshape lat and lon to allow plotting
    lat_rshp = repmat(lat,imagez,1);   % lat replicated to match number of elements
    lon_rshp = repmat(lon,imagem,1);   % lon replicated to match number of elements
    lat_rshp = reshape(lat_rshp,imagem,imagez)';
    lon_rshp = reshape(lon_rshp,imagez,imagem);
    % lat = lat_rshp;
    % lon = lon_rshp;
    
   
    %% Plot contour map 
    % Override lat/lon limit settings
    % Study Area
    latmin = 27;
    latmax = 31;
    lonmin = -94;
    lonmax = -87.5;

    % Entire Gulf
    % latmin = 18;
    % latmax = 32;
    % lonmin = -98;
    % lonmax = -80;
    
    alat = 0.5*(latmin+latmax);
    f = cos(pi*alat/180);

    % Vertical to horizontal ratio
    % [left bottom width height]
    v2h = (latmax-latmin)/((lonmax-lonmin).*f) ;
    if (v2h <= 1.1) 
	    fig_pos = [0.1 0.15 0.65./v2h .65];
    else
	    fig_pos = [0.1 .15 0.65 .65./v2h];
    end

    %% Plot figures
    
    % Loop through variables
    for  ivar = 1:length(iop_list)
        iop = iop_list{ivar};
        hf1 = figure(1);
        clf
        % scrsz = get(groot,'ScreenSize');
        % set(hf1,'Position',[scrsz(4).*.1 scrsz(3).*.1 scrsz(3).*.65 scrsz(4).*.65])
        % set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired
    
        latlim = [latmin latmax];
        lonlim = [lonmin lonmax];
        
        switch iop
            case 'kd'
                wv = wv_kd;
                log_kd = log10(kd490);
                log_kd(isnan(log_kd)) = -5;
                G = geoshow(lat_rshp,lon_rshp,log_kd,'DisplayType','surface');
                clim_min = 0.01;
                clim_max = 2;
                textlabel = '{\it{K_d}}(490) (m^{-1})';
            case 'aph'
                wv = wv_aph;
                log_aph = log10(aph440);
                log_aph(isnan(log_aph)) = -5;
                G = geoshow(lat_rshp,lon_rshp,log_aph,'DisplayType','surface');
                clim_min = 0.005;
                clim_max = 1.5;
                textlabel = '{\it{a_{ph}}}(440) (m^{-1})';
            case 'chlor_par'
                wv = 1;  % 
                log_chl = log10(chl);
                log_chl(isnan(log_chl)) = -5;
                G = geoshow(lat_rshp,lon_rshp,log_chl,'DisplayType','surface');
                clim_min = 0.01;
                clim_max = 30;
                textlabel = 'Chlor\_a (mg m^{-3})';
            case 'micro'
                wv = 1;  % 
                G = geoshow(lat_rshp,lon_rshp,fmicro,'DisplayType','surface');
                clim_min = 0.0;
                clim_max = 1;
                textlabel = 'Microplankton(Uitz)';
            % case 'rrs'
            %     G = geoshow(lat_rshp,lon_rshp,log10(rrs532),'DisplayType','surface');
            %     clim_min = 0.001;
            %     clim_max = 0.05;
            %     textlabel = '{\it{R_{rs}}}(532) (sr^{-1})';
        end
    
        ax = gca;
        set(ax,'Xlim',lonlim,'YLim',latlim,'FontSize',14,'LineWidth',1.25,'Position',fig_pos,'PositionConstraint','innerposition');
        drawnow
        hold on
        % 
        xlabel('Longitude (^{o}W)','Fontsize',16,'Fontname','Arial','Fontweight','Bold');
        ylabel('Latitude (^{o}N)','Fontsize',16,'Fontname','Arial','Fontweight','Bold');
        % what do i do if my archive file from Blackboard is too large to import to Canvas
        xtick_dat=get(gca,'Xtick');
        xtick_label=-xtick_dat;
        set(ax,'XtickLabel',xtick_label);
        % 
        if strcmp(iop,'micro')
            set(ax,'clim',[clim_min,clim_max]);
        else
            set(ax,'clim',log10([clim_min,clim_max]));
        end
    
        if hbar_flag
            hbar=colorbar;
        end

        hmap = colormap("parula"); %jet
        hmap(1,:) = [0,0,0];
        colormap(hmap);
    
        if strcmp(iop,'chlor_par')
            ctickvals = [0.01;0.1;0.5;1;5;10;30];  % clim_min+round([0.0;0.04;0.09;0.4;0.9].*(clim_max - clim_min),1);
        elseif strcmp(iop,'micro')
            ctickvals = 0.2:0.2:1;  % clim_min+round([0.0;0.04;0.09;0.4;0.9].*(clim_max - clim_min),1);
        else
            ctickvals = [0.01;0.03;0.1;0.2;0.5;1;2];  % clim_min+round([0.0;0.04;0.09;0.4;0.9].*(clim_max - clim_min),1);
        end
        
        % hbar_pos=get(hbar,'Position');
        % hbar_pos(1)=hbar_pos(1) +.02; %move off page if desired (set to 1.0)
        if strcmp(iop,'micro')
            if hbar_flag
                set(hbar,'YLim',[clim_min,clim_max],'Ytick',ctickvals,...
                    'fontname','arial','fontsize',14); %'YtickLabel',num2str(ctickvals),'Position',hbar_pos);
            end
        else
            if hbar_flag
                set(hbar,'YLim',[log10(clim_min) log10(clim_max)],'Ytick',[log10(ctickvals)],...
                 'YtickLabel',num2str(ctickvals),'fontname','arial','fontsize',14); %,'Position',hbar_pos);
            end
        end
        % hold on
    
        %Coastline data from http://www.soest.hawaii.edu/pwessel/gshhg/ (using high resolution)
        coastfile=[ancfolder,'GSHHS_f_L1.shp'];all_cruise_aph_size_class_tab

        map_info=shapeinfo(coastfile);
        coast_dat=shaperead(coastfile,'UseGeoCoords',true,'BoundingBox',[lonmin,latmin;lonmax,latmax]);  
        coast_n=size(coast_dat);
        
        % Load coast file
        % load([inputfolder,'GSHHS_f_L1_coast_dat.mat']);
        % coast_n=size(coast_dat);
        
        %Loop to plot Lat/Lon to current map axes
        for mapn=1:coast_n(1)
            coast_dat(mapn).Lon=coast_dat(mapn).Lon; 
            hcst = fill(coast_dat(mapn).Lon(1:end-1),coast_dat(mapn).Lat(1:end-1),...
                [0.7,0.6,0.3]); %Use fill command to plot filled polygons
            hcst.EdgeColor = [0,0,0];
            hcst.Marker = '.';
            hcst.MarkerSize = 0.1;
            hold on
        end
        
        box on;
    
        htxt1=text(ax,lonmin+0.5,latmax-0.3,textlabel,'fontname','arial','fontsize',16,'Color','k','FontWeight','bold');
        htxt2=text(ax,lonmin+0.5,latmax-0.6,dt_label,'fontname','arial','fontsize',15,'Color','k','FontWeight','normal');
    
        % Consolidate all wavelengths into a single 3D array
        img_dat = cat(3, vardata{var_indx});
    
        lat_vect = lat_rshp(:);
        lon_vect = lon_rshp(:);
    
        % Read in cluster stations for a given cruise and determine means (SEE NEXT SECTION)
    
        % Read in spreadsheet with Cluster station locations
        clstr_sht = readtable([ancfolder,'all_cruise_stations_clusters_fmicro.xlsx']);
        
        switch iop
            case 'kd'
                % Select indices for stations associated with each cluster
                cluster_name = {'Estuary','Inner','Mid','Outer'};
                sta_symbol = {'+','o','^','*'};
                
                sta_indx = [];  % Clear prior indices
                plt_indx = [];  % For indexing legend entries
                for iclst = 1:4
                    sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & contains(clstr_sht.Cluster,cluster_name{iclst}) & clstr_sht.Hpro_flag==1);
                    if isempty(sta_indx)
                        continue
                    end
                    plt_indx = [plt_indx,iclst];  % This limits indices to valid graphics objects
                    sta_loc = [clstr_sht.Lon(sta_indx),clstr_sht.Lat(sta_indx)];
                    htarg1(iclst) = plot3(sta_loc(:,1),sta_loc(:,2),repmat(10,size(sta_loc,1),1),sta_symbol{iclst},'MarkerSize',12,'LineWidth',2);
                    
                    % Add labels to station locations
                    sta_txt = clstr_sht.Station(sta_indx);
                    htxt1 = text(sta_loc(:,1)+0.1,sta_loc(:,2)-0.1,repmat(10,size(sta_loc,1),1),sta_txt,'Fontsize',14,'Color',[0.3 0.3 0.3]);
                    
                    sta_indx = [];
                end
                
                hleg_kd = legend(htarg1(plt_indx),cluster_name(plt_indx),'Location','southeast');
                
                % print([inputfolder,labeltext,'.tif'],'-dtiff','-r300');
                % print([inputfolder,labeltext,'_',iop,'.eps'],'-depsc','-r300');  % Encapsulated postscript color
            case 'aph'
                % Reshape matrix to vector format
                img_dat_vect = reshape(img_dat,imagem.*imagez,size(wv,1));
    
                % % Predictors (e.g., lambda)
                % wvfit_indx = find(wv(wv<=510 & wv>=440));
                % wvl_pred = wv(wvfit_indx);
                % X_reg = [ones(size(wvl_pred,1),1) wvl_pred];
                % aphfit_indx = find(~isnan(img_dat_vect(:,1)));
                % 
                % % Perform regression
                % coeffs(:,aphfit_indx) = X_reg \ img_dat_vect(aphfit_indx,wvfit_indx)'; % Linear fit for every pixel
                % fmicro_pace = 0.853./(1+exp(-1.93.*(log(abs(coeffs(2,:)))+8.32)));
                % micro_indx = find(log(abs(coeffs(2,:)))>-8.3);
                % pico_indx = find(log(abs(coeffs(2,:)))<=-8.3);
    
                % Select indices for stations associated with each cluster
                class_name = {'Micro','Pico_/Nano'};
                sta_symbol = {'+','o'};
    
                sta_indx = [];  % Clear prior indices
                plt_indx = [];  % For indexing legend entries
                for iclss = 1:2
                    if iclss == 1
                        sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & clstr_sht.fmicro>0.5 & clstr_sht.Hpro_flag==1);
                    else
                        sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & clstr_sht.fmicro<=0.5 & clstr_sht.Hpro_flag==1);
                    end
                    
                    if isempty(sta_indx)
                        continue
                    end
    
                    plt_indx = [plt_indx,iclss];  % This limits indices to valid graphics objects
                    sta_loc = [clstr_sht.Lon(sta_indx),clstr_sht.Lat(sta_indx)];
                    htarg2(iclss) = plot3(sta_loc(:,1),sta_loc(:,2),repmat(10,size(sta_loc,1),1),sta_symbol{iclss},'MarkerSize',12,'LineWidth',2);
        
                    % Add labels to station locations
                    sta_txt = clstr_sht.Station(sta_indx);
                    htxt1 = text(sta_loc(:,1)+0.1,sta_loc(:,2)-0.1,repmat(10,size(sta_loc,1),1),sta_txt,'Fontsize',14,'Color',[0.3 0.3 0.3]);
        
                    sta_indx = [];
    
                end
    
                hleg_aph = legend(htarg2(plt_indx),class_name(plt_indx),'Location','southeast');
    
            case 'micro'
                % Select indices for stations associated with each cluster
                class_name = {'Micro','Pico_/Nano'};
                sta_symbol = {'+','o'};

                sta_indx = [];  % Clear prior indices
                plt_indx = [];  % For indexing legend entries
                for iclss = 1:2
                    if iclss == 1
                        sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & clstr_sht.fmicro>0.5 & clstr_sht.Hpro_flag==1);
                    else
                        sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & clstr_sht.fmicro<=0.5 & clstr_sht.Hpro_flag==1);
                    end
    
                    if isempty(sta_indx)
                        continue
                    end
    
                    plt_indx = [plt_indx,iclss];  % This limits indices to valid graphics objects
                    sta_loc = [clstr_sht.Lon(sta_indx),clstr_sht.Lat(sta_indx)];
                    htarg2(iclss) = plot3(sta_loc(:,1),sta_loc(:,2),repmat(10,size(sta_loc,1),1),sta_symbol{iclss},'MarkerSize',12,'LineWidth',2);
    
                    % Add labels to station locations
                    sta_txt = clstr_sht.Station(sta_indx);
                    htxt1 = text(sta_loc(:,1)+0.1,sta_loc(:,2)-0.1,repmat(10,size(sta_loc,1),1),sta_txt,'Fontsize',14,'Color',[0.3 0.3 0.3]);
    
                    sta_indx = [];
                end
    
                hleg_micro = legend(htarg2(plt_indx),class_name(plt_indx),'Location','southeast');
                file_suffix = '_fmicro.tif';
                plt_ttl = 'fmicro';
        end
    
        % End program if plotting chlor_par
        if strcmp(iop,'chlor_par')
            disp('Completed');
            continue
        end
        
        % Calculate means for target areas (individual stations)
        
        switch iop
            case 'aph'
            
            var_area_mean = zeros(length(wv),2);
            var_area_std = zeros(length(wv),2);
            var_area_se = zeros(length(wv),2);
            micro_sta_indx = [];
            piconano_sta_indx = [];

            for iclss = 1:2
                micro_sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & clstr_sht.fmicro>0.5 & clstr_sht.Hpro_flag==1);
                piconano_sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & clstr_sht.fmicro<=0.5 & clstr_sht.Hpro_flag==1);
    
                if iclss == 1 & isempty(micro_sta_indx)
                    continue
                elseif iclss == 2 & isempty(piconano_sta_indx)
                    continue
                end
    
                if iclss == 1
                    sta_loc = [clstr_sht.Lon(micro_sta_indx),clstr_sht.Lat(micro_sta_indx)];
                elseif iclss == 2
                    sta_loc = [clstr_sht.Lon(piconano_sta_indx),clstr_sht.Lat(piconano_sta_indx)];
                end
    
                % Loop through station locations to select pixels for averaging
                nsta = size(sta_loc,1);
                % area_indx = [];
                % for stan = 1:nsta
                %     area_indx = [area_indx;find(lat_vect>sta_loc(stan,2)-0.1 & lat_vect<sta_loc(stan,2)+0.1 & ...
                %         lon_vect>sta_loc(stan,1)-0.1 & lon_vect<sta_loc(stan,1)+0.1)];
                % end
                % var_area_mean(:,iclss) = mean(img_dat_vect(area_indx,:),1,'omitmissing');
                % var_area_std(:,iclss) = std(img_dat_vect(area_indx,:),[],1,'omitmissing');
                % sta_indx = [];
    
                switch sensor_type
                    % Compute means for each station, normalizing to aph(440)
                    case 'PACE'
                        matchup_interval = 0.025;
                        station_means = [];
                        station_std = [];
                        station_se = [];
                        for stan = 1:nsta
                            px_indx = find(lat_vect>sta_loc(stan,2)-matchup_interval & lat_vect<sta_loc(stan,2)+matchup_interval & ...
                                lon_vect>sta_loc(stan,1)-matchup_interval & lon_vect<sta_loc(stan,1)+matchup_interval);
                            station_means = [station_means;mean(img_dat_vect(px_indx,:)./img_dat_vect(px_indx,(wv==440 | wv==443)),1,'omitmissing')]; % 
                            station_std = [station_std;std(img_dat_vect(px_indx,:)./img_dat_vect(px_indx,(wv==440 | wv==443)),[],1,'omitmissing')];
                            station_se = [station_se;station_std./sqrt(length(px_indx))];
                        end
                        var_area_mean(:,iclss) = mean(station_means,1,'omitmissing');
                        var_area_std(:,iclss) = std(station_means,[],1,'omitmissing');
                        var_area_se(:,iclss) = sqrt(sum(station_se.^2)./nsta);
                        sta_indx = [];
                end
            end
            
            class_name = {'Micro','Pico/Nano'};
            file_suffix = '_aph.tif';
            plt_ttl = ' \it{a_{ph}}';
        
            % Normalize aph spectra to aph440 to get shape vectors\
            % var_area_mean = var_area_mean./var_area_mean(wv==440,:);
        case 'kd'
            % Cluster means for kd
        
            % Normalize to kd490 and reshape matrix to vector format
            % img_dat_vect = reshape(img_dat./kd490,imagem.*imagez,size(wv,1));
            
            % Unnormalized 
            img_dat_vect = reshape(img_dat,imagem.*imagez,size(wv,1));
    
            var_area_mean = zeros(length(wv),4);
            var_area_std = zeros(length(wv),4);
            var_area_se = zeros(length(wv),4);

            switch sensor_type
                % Compute means for each station
                case 'PACE'
                    % sta_indx = [];
                    % for iclst = 1:4
                    %     sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & contains(clstr_sht.Cluster,cluster_name{iclst}) & clstr_sht.Hpro_flag==1);
                    %     if isempty(sta_indx)
                    %         continue
                    %     end
                    %     sta_loc = [clstr_sht.Lon(sta_indx),clstr_sht.Lat(sta_indx)];
                    % 
                    %     % Loop through station locations to select pixels for averaging
                    %     nsta = size(sta_loc,1);
                    %     area_indx = [];
                    %     for stan = 1:nsta
                    %         area_indx = [area_indx;find(lat_vect>sta_loc(stan,2)-0.1 & lat_vect<sta_loc(stan,2)+0.1 & ...
                    %          lon_vect>sta_loc(stan,1)-0.1 & lon_vect<sta_loc(stan,1)+0.1)];
                    %     end
                    %     var_area_mean(:,iclst) = mean(img_dat_vect(area_indx,:),1,'omitmissing');
                    %     var_area_std(:,iclst) = std(img_dat_vect(area_indx,:),[],1,'omitmissing');
                    %     sta_indx = [];
                    % end

                % Compute means for each station, normalizing to kd(490)
                matchup_interval = 0.025;
                sta_indx = [];
                for iclst = 1:4
                    sta_indx = find(contains(clstr_sht.Cruise,cruise_name) & contains(clstr_sht.Cluster,cluster_name{iclst}) & clstr_sht.Hpro_flag==1);
                    if isempty(sta_indx)
                        continue
                    end
                    sta_loc = [clstr_sht.Lon(sta_indx),clstr_sht.Lat(sta_indx)];
        
                    % Loop through station locations to select pixels for averaging
                    nsta = size(sta_loc,1);
        
                    station_means = [];
                    station_std = [];
                    station_se = [];
                    for stan = 1:nsta
                        px_indx = find(lat_vect>sta_loc(stan,2)-matchup_interval & lat_vect<sta_loc(stan,2)+matchup_interval & ...
                            lon_vect>sta_loc(stan,1)-matchup_interval & lon_vect<sta_loc(stan,1)+matchup_interval);
                        station_means = [station_means;mean(img_dat_vect(px_indx,:),1,'omitmissing')]; % ./img_dat_vect(px_indx,wv==490)
                        station_std = [station_std;std(img_dat_vect(px_indx,:),[],1,'omitmissing')];
                        station_se = [station_se;station_std./sqrt(length(px_indx))];
                    end
                    var_area_mean(:,iclst) = mean(station_means,1,'omitmissing');
                    var_area_std(:,iclst) = std(station_means,[],1,'omitmissing');
                    var_area_se(:,iclst) = sqrt(sum(station_se.^2)./nsta);
                    sta_indx = [];
                end
            end

            file_suffix = '_kd.tif';
            plt_ttl = ' Cluster Averages';
        end
        
        if strcmp(sensor_type,'PACE')
            switch iop
                case 'kd'
                    var_area_mean_all_cruise = [var_area_mean_all_cruise,...
                    table(wv,var_area_mean(:,1),var_area_std(:,1),var_area_se(:,1),var_area_mean(:,2),var_area_std(:,2),...
                    var_area_se(:,2),var_area_mean(:,3),var_area_std(:,3),var_area_se(:,3),var_area_mean(:,4),var_area_std(:,4),...
                    var_area_se(:,4),'VariableNames',{[mnth{1},'_Wavelength'],[mnth{1},'_',cluster_name{1}],...
                    [mnth{1},'_',cluster_name{1},'_std'],[mnth{1},'_',cluster_name{1},'_se'],[mnth{1},'_',cluster_name{2}],[mnth{1},'_',cluster_name{2},'_std'],...
                    [mnth{1},'_',cluster_name{2},'_se'],[mnth{1},'_',cluster_name{3}],[mnth{1},'_',cluster_name{3},'_std'],[mnth{1},'_',cluster_name{3},'_se'],...
                    [mnth{1},'_',cluster_name{4}],[mnth{1},'_',cluster_name{4},'_std'],[mnth{1},'_',cluster_name{4},'_se']})];   
                    [max_var,max_indx] = max(var_area_mean,[],1);
                    [min_var,min_indx] = min(var_area_mean,[],1);
                case 'aph'
                    varnames = {[mnth{1},'_Lambda'],[mnth{1},'_Micro'],[mnth{1},'_Micro_std'],[mnth{1},'_Micro_se'],[mnth{1},'_Pico-Nano'],...
                        [mnth{1},'_Pico-Nano_std'],[mnth{1},'_Pico-Nano_se']};
                    var_area_mean_all_cruise = [var_area_mean_all_cruise,...
                        table(wv,var_area_mean(:,1),var_area_std(:,1),var_area_se(:,1),var_area_mean(:,2),var_area_std(:,2),var_area_se(:,2),...
                        'VariableNames',varnames)];
                    [max_var,max_indx] = max(var_area_mean,[],1);
                    [min_var,min_indx] = min(var_area_mean,[],1);
            end
        end

        print([datafolder,labeltext,file_suffix],'-dtiff','-r300');
    
        % End program if plotting micro
        if strcmp(iop,'micro') | strcmp(sensor_type,'MODIS')
            disp('Completed');
            continue
        end
    
        %% Plot spectral Kd or aph
    
        figure(2);
        clf
        
        h1 = plot(wv,var_area_mean,'-','LineWidth',2);
        hold on
        box on
        
        xlabel('Wavelength (nm)','FontSize',24,'FontWeight','bold');
        ax = gca;
        
        switch prod_type
            case 'kd'
                ylabel('{\it{K_d(\lambda)}} (m^{-1})','FontSize',14,'FontWeight','bold');  % /{\itK_d(490)}}
                ylimits = [0.02 6]; %[0.3 50]; 
                lab_xy = [wv(max_indx(1))-10,1.1.*max_var(1,1);wv(max_indx(2))-10,1.1.*max_var(1,2);...
                    wv(max_indx(3))-10,1.1.*max_var(1,3);wv(1)-10,1.1.*var_area_mean(1,4)];
                set(ax,'YScale','log','Fontsize',22,'YTick',[0.01,0.03,0.1,0.3,1,3,6,10,20,50], ...
                            'XLim',[390 710],'YLim',ylimits,'XTick',400:100:800,'Linewidth',1,'FontWeight','bold');
                pbaspect([4 3 1]);
                labeltext = [labeltext,'_cluster_averages'];
            case 'aph'
                ylabel('{\it{a_{ph}(\lambda)}} Shape Vector','FontSize',14,'FontWeight','bold');
                ylimits = [0.0 1.5];
                lab_xy = [wv(max_indx(1))-10,1.2.*max_var(1,1);wv(max_indx(2))-10,1.2.*max_var(1,2)];
                set(ax,'YScale','linear','Fontsize',16,'YTick',0.2:0.2:1, ...
                    'XLim',[390 550],'YLim',ylimits);
                labeltext = [labeltext,'_size_class_averages'];
            case 'rrs'
                ylabel('{\it{R_{rs}(\lambda)}} (m^{-1})','FontSize',14,'FontWeight','bold');
                ylimits = [0.005 0.1];
                lab_xy = [wv(max_indx(1))-10,1.2.*max_var(1,1);wv(max_indx(2))-10,1.2.*max_var(1,2);...
                    wv(max_indx(3))-10,1.2.*max_var(1,3);wv(max_indx(4))-10,1.2.*max_var(1,4)];
                set(ax,'YScale','log','Fontsize',16,'YTick',[0.01,0.03,0.1,0.3,1,3,6], ...
                            'XLim',[350 750],'YLim',ylimits);
        end % End switch prod_type       
        
        % Add labels to curves
        switch iop
            case 'kd'
                plt_labels = cluster_name;
            case 'aph'
                plt_labels = class_name;
        end
        % text(gca,lab_xy(:,1),lab_xy(:,2),plt_labels,'Fontsize',14)    
        hleg = legend(plt_labels,'FontSize',16,'Location','Southeast');
        
        htl = title(ax,[char(mnth),plt_ttl],'FontSize',18); 
        
        print([datafolder,labeltext,file_suffix],'-dtiff','-r300');
        
    end % End month range loop
end

if strcmp(sensor_type,'PACE')
    switch iop
        case 'kd'
            writetable(var_area_mean_all_cruise,[datafolder,'kd_cluster_vector_means_PACE_monthly_composites.xlsx']);
        case 'aph'
            writetable(var_area_mean_all_cruise,[datafolder,'aph_size_class_means_PACE_monthly_composites.xlsx']);
    end
end

disp('Completed');




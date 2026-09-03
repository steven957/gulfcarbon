% plot_pace_integrated_PP.m
% Syntax:  plot_pace_pace_integrated_PP
%
% Script reads PACE .mat file and generated with pp_calculation_image_master
% plots water column-integrated primary production (gC m-2 h-1)
% products
%
% Inputs:
%    1) Directory locations for .mat and coast files
%    2) *.mat file generated with pp_calculation_image_master
%    3) High resolution coastal shape file ('GSHHS_f_L1.shp')
%
% Outputs:
%    1) map plot of water-column-integrated primary production
%    2) *.tif file of above plot (if desired)
%   
% Other m-files required: None 
%
% MAT-files required: 
%    1) *.mat file generated with pp_calculation_image_master
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 7 Aug 2026

%% ------------- BEGIN CODE --------------%

clc
clearvars

sensor_type = 'PACE'; % 'MODIS';

monthrange = 1; %1:12;

all_pp_mean = table();

for mnthn = monthrange

    %Input image data and plot
    switch sensor_type
        case 'PACE'
            inputfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/pace/';
            pp_file_lst=dir([inputfolder,'PACE_OCI.2025',sprintf('%02d', mnthn),'*.L3.PP.mat']);
        case 'MODIS'
            inputfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/modis/';
            pp_file_lst=dir([inputfolder,'A2025',sprintf('%02d', mnthn),'*.L3.PP.mat']);
    end
    
    coastfolder='/run/media/slohrenz/PrimaryDrive/ocean_color/pace/Ancillary/';  %Enter the folder where your coastal shape file can be found
    
    latmin = 25;
    latmax = 31;
    lonmin = -96;
    lonmax = -86;
    
    [filen,m]=size(pp_file_lst);
    
    mean_pp = table('Size',[filen,3],'VariableNames',{'Date','Plume','Outer Shelf'},'VariableTypes',["datetime","double","double"]);
    
    for ifile=1:filen
        
        % Load data
        pp_file=char(pp_file_lst(ifile).name);
        switch sensor_type
            case 'PACE'
                date_text = [pp_file(10:17)]; %,' ',pp_file(14:17),'-',pp_file(23:26)];
            case 'MODIS'
                date_text = [pp_file(2:9)]; %,' ',pp_file(14:17),'-',pp_file(23:26)];
        end
        dt = datetime(date_text,'InputFormat','uuuuMMdd');
        mnth = char(month(dt,'name'));
        yr = num2str(year(dt));
    
        disp(['Plotting ',pp_file]);
        
        pp_path=[inputfolder,pp_file];
        load(pp_path);
    
        %% Determine average production for target areas
    
        bxwdth = 0.2;  % Width of box in degrees
    
        % River plume
        lat_min1 = 28.8;
        lat_max1 = lat_min1 + bxwdth;
        lon_min1 = -89.85;
        lon_max1 = lon_min1 + bxwdth;
    
        pp_indx1 = find(lat_rshp<=lat_max1 & lat_rshp>=lat_min1 & lon_rshp<=lon_max1 & lon_rshp>=lon_min1);
        pp_mean1 = mean(PP_est_int_array(pp_indx1),"omitmissing").*12; % Units converted to gC m-2 d-1
    
        % Outer shelf
        lat_min2 = 28.0;
        lat_max2 = lat_min1 + bxwdth;
        lon_min2 = -91.5;
        lon_max2 = lon_min2 + bxwdth;
    
        pp_indx2 = find(lat_rshp<=lat_max2 & lat_rshp>=lat_min2 & lon_rshp<=lon_max2 & lon_rshp>=lon_min2);
        pp_mean2 = mean(PP_est_int_array(pp_indx2),'omitmissing').*12; % Units converted to gC m-2 d-1
        
        mean_pp(ifile,:)={dt,pp_mean1,pp_mean2};

        %% Plot data
        % close all
        figure(1);
        clf
        set(gcf, 'renderer', 'zbuffer')
        
        latlim = [latmin latmax];
        lonlim = [lonmin lonmax];
    
        % G = geoshow(lat,lon,log10(PP_est_int_array.*12),'DisplayType','surface'); % Factor of 12 converts from molC to gC
        G = geoshow(lat_rshp,lon_rshp,log10(PP_est_int_array.*12),'DisplayType','surface'); % Factor of 12 converts from molC to gC
        % G = imagesc(lon,lat,log10(rot90(PP_est_int_array).*12)); % Factor of 12 converts from molC to gC
        % set(gca,'YDir','reverse'); 
        clim_min = 0.1;
        clim_max = 10;
        textlabel = 'Integrated PP (gC m^{-2} d^{-1})';
    
        ax = gca;
        set(ax,'Xlim',[lonmin,lonmax],'YLim',[latmin,latmax],'FontSize',18);
        aspect_rat = cosd((latmin+latmax)./2);
        ax.DataAspectRatio = [1 aspect_rat 1];
        box on 
        hold on
    
        % 
        xlabel('Longitude (^{o}W)','Fontsize',20,'Fontname','Arial','Fontweight','Bold');
        ylabel('Latitude (^{o}N)','Fontsize',20,'Fontname','Arial','Fontweight','Bold');
        % 
        xtick_dat=get(gca,'Xtick');
        xtick_label=-xtick_dat;
        set(ax,'XtickLabel',xtick_label);
        % 
        set(ax,'clim',log10([clim_min,clim_max]));
        hbar=colorbar;
        hmap = colormap(jet);
        colormap(hmap);
    
        % ctickvals = clim_min+round([0.0;0.05;0.10;0.5;1].*(clim_max - clim_min),-2);
        ctickvals = [0.010;0.02;0.05;0.1;0.2;0.5;1;2;5;10];
        % hbar_pos=get(hbar,'Position');
        % hbar_pos(1)=hbar_pos(1) +.02; %move off page if desired (set to 1.0)
        set(hbar,'YLim',[log10(clim_min) log10(clim_max)],'Ytick',[log10(ctickvals)],...
             'YtickLabel',num2str(ctickvals),'fontname','arial','fontsize',14); %,'Position',hbar_pos);
        
        % Adjust axes position
        axpos = get(ax,'Position');
        axpos(2) = axpos(2)+0.05;
        set(ax,'Position',axpos);
    
        % Plot coastline (NOTE: THIS SECTION OF COMMANDS ONLY NECESSARY FOR INTIAL RUN TO SAVE COAST FILE) 
        % coastfile=[coastfolder,'GSHHS_f_L1.shp'];
        % %Coastline data from http://www.soest.hawaii.edu/pwessel/gshhg/ (using low resolution)
        % map_info=shapeinfo(coastfile);
        % coast_dat=shaperead(coastfile,'UseGeoCoords',true,'BoundingBox',[lonmin,latmin;lonmax,latmax]);  
    
        % Save coast file for later use
        % save([inputfolder,'GSHHS_f_L1_coast_dat.mat'],'coast_dat');
        
        % Load coast file
        load([coastfolder,'GSHHS_f_L1_coast_dat.mat']);
        coast_n=size(coast_dat);
        
        %Loop to plot Lat/Lon to current map axes
        for mapn=1:coast_n(1)
            % coast_dat(mapn).Lon=coast_dat(mapn).Lon; 
            fill(coast_dat(mapn).Lon(1:end-1),coast_dat(mapn).Lat(1:end-1),[0.8,0.6,0.4]); %Use fill command to plot filled polygons
            hold on
        end
        
        htxt1=text(ax,lonmin+0.5,latmax-0.25,textlabel,'fontname','arial','fontsize',16,'Color','k','FontWeight','bold');
        htxt2=text(ax,lonmin+0.5,latmax-0.5,[mnth,' ',yr],'fontname','arial','fontsize',16,'Color','k','FontWeight','bold');
    
        % hrect1 = rectangle(ax,'position',[lon_min1,lat_min1,0.2,0.2],'LineWidth',1.25,'EdgeColor',[0.1 0.1 0.1]);
        % hold on
        
        % Define corners of the rectangle in 3D space [x y z]
        xrct1 = [lon_min1 lon_min1+bxwdth lon_min1+bxwdth lon_min1 lon_min1];
        yrct1 = [lat_min1  lat_min1 lat_min1+bxwdth lat_min1+bxwdth lat_min1];
        zrct1 = [0.999 0.999 0.999 0.999 0.999]; % Set Z-level
        
        hrect1 = plot3(xrct1,yrct1,zrct1,'-k','Linewidth',1.5);
        
        
        xrct2 = [lon_min2 lon_min2+bxwdth lon_min2+bxwdth lon_min2 lon_min2];
        yrct2 = [lat_min2  lat_min2 lat_min2+bxwdth lat_min2+bxwdth lat_min2];
        zrct2 = [0.999 0.999 0.999 0.999 0.999]; % Set Z-level
        
        hrect2 = plot3(xrct2,yrct2,zrct2,'-k','Linewidth',1.5);
    
        % patch(xrct1,yrct1,zrct1,'LineWidth',1.75,'EdgeColor',[0 0 0],'FaceColor','none','Parent',ax);
        % hrect2 = rectangle(ax,'position',[lon_min2,lat_min2,0.2,0.2],'LineWidth',1.75);
    
        print([char(inputfolder),sensor_type,date_text,'_pp.tif'],'-dtiff','-r300');
    
    end

    all_pp_mean = [all_pp_mean;mean_pp];

end

writetable(all_pp_mean,[inputfolder,'mean_pp_table_',sensor_type,'.xls']);

disp('Completed');




% areamap_modis_true_color_selected_stations_PE.m
%
% Script plots map of northern Gulf of Mexico with
%   photosynthesis-irradiance stations and HyperPro profile locations 
%   overlaid on true color satellite image
%
% Inputs:
%    1) True color satellite image for cruise
%    2) Spreadsheet with station names and locations 
%
% Outputs:
%    1) tif image of map
%   
% Other m-files required: None
%
% MAT-files required: Gulf of Mexico coastline file
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 9 Dec 2025

%% ------------- BEGIN CODE --------------%

clc
clear variables
clf

%Load data 

plot_sta=1;  %1 for plotting stations, 0 for none

% Choose whether to plot all stations or subset
sta_plot_option='subset';   %'all';  %'subset';
    
fig=figure(1);
scrsz = get(groot,'ScreenSize');
set(fig,'Position',[scrsz(4).*.1 scrsz(3).*.1 scrsz(3).*.65 scrsz(4).*.65])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

cruise='January';  %'January';  %'March';  %'November';  %'July';  %'April';
coastpath = '';  % Location of coastline data

% For river mouth enlargement
% latmin = 28.5;     % Minimum and Maximum Latitudes and Longitudes
% latmax = 29.25;
% lonmin = -90;
% lonmax = -89;
% tickincr=0.2;

% Entire study area
latmin = 27;     % Minimum and Maximum Latitudes and Longitudes
latmax = 31;
lonmin = -94;
lonmax = -87.5;
tickincr=1;

alat = 0.5*(latmin+latmax);
f = cos(pi*alat/180);

% Vertical to horizontal ratio
%  [left bottom width height]

v2h = (latmax-latmin)/((lonmax-lonmin).*f) ;
if (v2h <= 1.1) 
	position = [.125 0.125 0.7 .6./v2h];
else
	position = [.125 .125 0.75.*v2h .75];
end

position2=position;

xticklab=-(lonmin:tickincr:lonmax);
yticklab=floor(latmin):tickincr:ceil(latmax);

% Set up map axes
h4=axes('Position',position2,'XColor','k','YColor','k','Color','none');
set(h4,'ylim',[latmin latmax],'xlim',[lonmin lonmax],...
    'ytick',floor(latmin):tickincr:ceil(latmax),...
    'xtick',lonmin:tickincr:ceil(lonmax),'Linewidth',1.2,'GridColor',[0.7,0.7,0.7],...
    'fontsize',18,'Layer','top','xticklabel',xticklab,'yticklabel',yticklab,...
    'Layer','top','DataAspectRatio',[1 f 1]);

pbaspect('manual');

hold on
hxlab=xlabel('{Longitude (^{o}W)}','Fontsize',22,'fontweight','bold','fontname','Helvetica');
hylab=ylabel('{Latitude (^{o}N)}','Fontsize',22,'fontweight','bold','fontname','Helvetica');

ylab_pos=get(hylab,'Position');
ylab_pos(1)=ylab_pos(1);
set(hylab,'Position',ylab_pos);

xlab_pos=get(hxlab,'Position');
xlab_pos(2)=xlab_pos(2);
set(hxlab,'Position',xlab_pos);

grid on
hold on

% Plot image if desired 

disp('Loading data...');
switch cruise

    case 'January'
        inpath='';  % Location of satellite image file
        [A,map]=imread([inpath,'A2009014191000_sub.RGB_QKM.png']);  % Aqua satellite image
        y=[33.75 25.22];
        x=[-98.01 -83.98];
        filenums=[12,15];  % DOY
        station_file='';  % Spreadsheet with station names and locations
        station_sheet='';
        file_txt = 'PE'; 
    case 'April'
        inpath='';  % Location of satellite image file
        [A,map]=imread([satpath,'A2009104194500_sub.RGB_QKM_2.png']);  % Aqua satellite image
        y=[32 26];
        x=[-95 -86];
        filenums=110:120;   % DOY 
        station_file='';  % Spreadsheet with station names and locations
        station_sheet='';
        file_txt ='PE';
    case 'July'
        inpath='';  % Location of satellite image file
        [A,map]=imread([inpath,'T2009202161500_sub.RGB_QKM.png']);  % Terra satellite image
        y=[32 26];
        x=[-95 -86];
        filenums=201:210;  % DOY
        station_file=''; % Spreadsheet with station names and locations
        file_txt = 'PE'; 
        station_sheet='';
    case 'November'
        inpath='';  % Location of satellite image file
        [A,map]=imread([inpath,'A2009305194000_sub.RGB_QKM.png']);  % Aqua satellite image
        y=[32 26];
        x=[-95 -86];
        filenums=302:310;  % DOY
        station_file=''; % Spreadsheet with station names and locations
        station_sheet='';
        file_txt = 'PE'; 
    case 'March'
        inpath='';  % Location of satellite image file
        [A,map]=imread([inpath,'A2010065190500_sub.RGB_QKM.png']);  % Aqua satellite image
        y=[32 26];
        x=[-95 -86];
        filenums=70:79;  % DOY
        station_file=''; % Spreadsheet with station names and locations
        station_sheet='';
        file_txt = 'PE'; 
end

% Plot image if desired 
himage=image(x,y,A);

% Plot the coastline, if available
load([coastpath,'usse_gom_coastline.mat']);  % Coastline file
nlat=coast(:,2);
nlon=coast(:,1);

plot(nlon,nlat,'MarkerSize',0.05,'Parent',h4,'Color',[.2 .2 .2],'Marker','.');
box on
hold on

% Plot bathymetry
load([coastpath,'usse_gom_bath2.mat']);

v=[20 100 500 1000 6000];
[~,hcon]=imcontour(lon, lat, F3, v);
set(hcon,'LineWidth',1.8,'LineStyle',':','EdgeColor',[.6,.6,.6]);

text(-91.4, 28.95,'20 m','Color',[0.6,0.6,0.6],'fontname','Helevetica',...
   'fontsize',8,'Rotation',0,'fontangle','italic','fontweight','normal');
text(-91.4, 28.28,'100 m','Color',[0.6,0.6,0.6],'fontname','Helevetica',...
   'fontsize',8,'Rotation',0,'fontangle','italic','fontweight','normal');
text(-91.4, 27.93,'500 m','Color',[0.6,0.6,0.6],'fontname','Helevetica',...
   'fontsize',8,'Rotation',0,'fontangle','italic','fontweight','normal');
text(-91.4, 27.67,'1000 m','Color',[0.6,0.6,0.6],'fontname','Helevetica',...
   'fontsize',8,'Rotation',5,'fontangle','italic','fontweight','normal');

% Plot additional annotation

annot_color=[1.0,1.0,0.5];

text(-89.5, 29.55,'Mississippi River','Color',annot_color,'fontname','Helevetica',...
   'fontsize',12,'Rotation',0,'fontangle','italic','fontweight','bold');
text(-89.45, 29.40,'("birdfoot" delta)','Color',annot_color,'fontname','Helevetica',...
   'fontsize',12,'Rotation',0,'fontangle','italic','fontweight','bold');

text(-92.1, 30,'Atchafalaya River','Color',annot_color,'fontname','Helevetica',...
   'fontsize',12,'Rotation',0,'fontangle','italic','fontweight','bold');
hline2=line([-91.6,-91.5],[29.9,29.45],'Color',annot_color,'Linewidth',1);

text(-90.28, 29.93,'New','Color',annot_color,'fontname','Helevetica',...
   'fontsize',12,'Rotation',0,'fontangle','italic','fontweight','bold');
text(-90.34, 29.82,'Orleans','Color',annot_color,'fontname','Helevetica',...
   'fontsize',12,'Rotation',0,'fontangle','italic','fontweight','bold');

text(-88.7, 27.4,'Clouds','Color',[0.6,0.6,0.6],'fontname','Helevetica',...
   'fontsize',10,'Rotation',0,'fontangle','italic','fontweight','bold');

hold on

% Add north compass arrow and scale
hnarrow=annotation('arrow',[0.2,0.2],[0.855,0.955],'Linewidth',2.0,'Color',annot_color);
text(-93.5,30.9,'N','Color',annot_color,'fontname','Helevetica',...
   'fontsize',10,'Rotation',0,'fontangle','italic','fontweight','bold');

axpos=get(gca,'Position');

scaledist=100; % km
londiff=lonmax-lonmin;
xwidth=axpos(3)-axpos(1);
xdist=londiff.*60*1.852*cos(30*pi()/180);
kmscale=(xwidth*scaledist/xdist);

hscale=annotation('doublearrow',[axpos(1)+0.15,axpos(1)+0.15+kmscale],[0.935,0.935],'Linewidth',2.0,'Color',annot_color);
text(-92.45,30.7,[num2str(scaledist),' km'],'Color',annot_color,'fontname','Helevetica',...
   'fontsize',10,'Rotation',0,'fontangle','italic','fontweight','bold');

% Plot stations if desired
if plot_sta==1

    % Load cruise track and station data (the next several statements are
    %   customized to the station waypoint file format - revise as needed
    latlon_tabl = readtable(station_file,'Sheet',station_sheet);
    
    latlon_char=char(latlon_tabl.EndingWaypointLat_Lon);
    sta_lab = latlon_tabl.LegNo__Name;

    latdeg=latlon_char(:,1:2);
    latmin=latlon_char(:,4:10);

    londeg=latlon_char(:,16:18);
    lonmin=latlon_char(:,21:26);

    sta_lat=str2num(latdeg(:,:))+str2num(latmin(:,:))./60;
    sta_lon=str2num(londeg)+str2num(lonmin)./60;

    % Strip unneeded characters from station label
    [txt,remain]=strtok(sta_lab);
    [sta_lab,txt2]=strtok(remain);

    switch sta_plot_option
        case 'all'    
            % Find waypoint labels for all stations
            wp_indx=find(strcmp('Point',sta_lab)~=1);  % Comment out if only specific stations desired

            sta_lab=sta_lab(wp_indx);
            
            % Plots selected stations
            plot(-sta_lon(wp_indx),sta_lat(wp_indx),'+y','Linewidth',1.5); 
            
            m=length(sta_lab);

            for i=1:m
                htext(i)=text(-sta_lon(wp_indx(i))+.075,sta_lat(wp_indx(i))-0.02,sta_lab(i),...
                    'Fontsize',13,'Linewidth',2,'Color','y');
            end
            
        case 'subset'

            % Get P-E stations
            
            switch cruise
                case 'January'
                    % GC1 P-E stations
                    subset_array1={'A1';'MR3';'C1';'D4';'D3';'E1';'F5';'H2';'G4'}; 
            
                    % GC1 HyperPro stations
                    subset_array2={'A1'};  
                case 'April'
                    % GC2 P-E stations
                    subset_array1={'B5';'B4';'C5';'D5';'F5';'H2';'H3';'G4';'G3';'A3';'A2';'A1'}; 
            
                    % GC2 HyperPro stations
                    subset_array2={'A1';'B5';'B4';'C5';'D5';'E1';'F5';'F4';'H2';'H3';'G4';'G2';'A3';'A2'};  
                case 'July'
                    % GC3 P-E stations
                    subset_array1={'A2';'B4';'B3';'C2';'D4';'D3';'F5';'F4';'G4';'G3'}; 

                    % GC3 HyperPro stations
                    subset_array2={'A1';'A2';'A3';'MR3';'B4';'B3';'C2';'C3';'D4';'D3';'D2';'E0';...
                        'E3';'F5';'F4';'H1';'H2';'H3';'G4';'G3';'NGI8'};  

                case 'November'
                    % GC4 P-E stations
                    subset_array1={'A1';'A2';'C2';'D4';'D3';'E2';'H2';'H3';'G4'}; 

                    % GC4 HyperPro stations
                    subset_array2={'A1';'A2';'C0';'C2';'D4';'D3';'D2';'E1';'E2';'F4';'F5';'G4';'G3';...
                        'H2';'H3';'MR3';'H4'};  

                case 'March'
                    % GC5 P-E stations
                    subset_array1={'C1';'C2';'D4';'E2';'E3';'F5';'F4';'H3';'G4';'G3';'B1';'B2'}; 

                    % GC5 HyperPro stations
                    subset_array2={'A2';'B2';'C1';'C2';'D3';'D4';'F4';'F5';'G2';'G3';'G4';'H2';'H3';'H4'};  

                    uwsubset_array={'UWS1','UWS2'};
            end
            %subset_array={'mr1','mr2','mr3','b2'};
            
            for isub = 1:2
                switch isub
                    case 1
                        subset_array = subset_array1;
                    case 2
                        subset_array = subset_array2;
                end

                wp_indx=[];
                wp_count=0;
                for nsub=1:length(subset_array)
                    % Find indices for specific stations if subset desired
                    if isempty(find(strcmpi(subset_array{nsub},sta_lab), 1))
                        continue
                    else
                        wp_count=wp_count+1;
                        wp_indx(wp_count)=find(strcmpi(subset_array{nsub},sta_lab)==1); %Case insensitive strcmp
                    end
                end
             
                sta_lab_sub=sta_lab(wp_indx);
                
                % Plots selected stations
                switch isub
                    case 1
                        hp1 = plot(-sta_lon(wp_indx),sta_lat(wp_indx),'+y','Linewidth',1.5,'MarkerSize',10); 
                        hold on
                    case 2
                        hp2 = plot(-sta_lon(wp_indx),sta_lat(wp_indx),'oy','Linewidth',1.5,'MarkerSize',10); 
                end
                
                m=length(sta_lab_sub);
    
                for i=1:m
                    htext(i)=text(-sta_lon(wp_indx(i))+.075,sta_lat(wp_indx(i))-0.02,upper(sta_lab_sub(i)), ...
                        'Fontsize',13,'Linewidth',2,'Color','y');
                end
    
                switch isub
                    case 2
                        hleg = legend([hp1,hp2],{'P-E stations','HyperPro stations'},'Fontsize',12,'Linewidth',1.5,...
                        'Color','#77AC30','TextColor','y');
                end                    
            end

    end

    hold off

end

exportgraphics(fig,[inpath,'rrs_true_color_modis_',cruise,'2009_',file_txt,'_stns.tif'],'Resolution',300,...
    'ColorSpace','rgb')

disp('Completed');

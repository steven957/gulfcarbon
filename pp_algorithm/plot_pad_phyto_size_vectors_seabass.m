% plot_pad_phyto_size_vectors_seabass.m - Plot selected filterpad data compiled from SeaBASS files 
% Syntax: plot_pad_phyto_size_vectors_seabass.m
%
% Inputs:
%    1) Folder location with compiled filterpad data ('all_pad.mat') produced 
%      using 'pad_extract_seabass_v2.m')
%    
% Outputs:
%    1) Figures with plotted filterpad results
%    2) '*.mat' files with fpico and fmicro shape vectors
%   
% Other m-files required: None
%  
% MAT-files required: 'all_pad.mat' (Compiled filterpad data from SeaBASS
%   files using 'pad_extract_seabass_v2.m')
%
% Author: Steven E. Lohrenz
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 30 Sep 2025

%% ------------- BEGIN CODE --------------%% 

clear all
clc
clf

% Folder path and file information for each cruise

for igc = 1:5

    cruise_name = ['gulfcarbon',num2str(igc)];

    switch cruise_name
        case 'gulfcarbon1'
            padpath='\Jan2009\Filterpad\To_SeaBASS\';
            pigsheetpath='\GC1_Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC1';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_jan2009.xlsx'],...
                'Sheet',wksheet);
            file_pref='GulfCarbon_gulfcarbon1_filterpad_';
             
        case 'gulfcarbon2'
            padpath='\Apr2009\Filterpad\To_SeaBASS\';
            pigsheetpath='\Apr2009\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC2';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_apr2009.xlsx'],...
                'Sheet',wksheet);
            file_pref='GulfCarbon_gulfcarbon2_filterpad_';
    
        case 'gulfcarbon3'
            padpath='\Jul2009\Filterpad\To_SeaBASS\';
            pigsheetpath='\Jul2009\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC3';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_jul2009.xlsx'],...
                'Sheet',wksheet);
           file_pref='GulfCarbon_gulfcarbon3_filterpad_';
            
        case 'gulfcarbon4'
            padpath='\Nov2009\Filterpad\To_SeaBASS\';
            pigsheetpath='\Nov2009\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC4';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_nov2009.xlsx'],...
                'Sheet',wksheet);
           file_pref='GulfCarbon_gulfcarbon4_filterpad_';
            
        case 'gulfcarbon5'
            padpath='\Mar2010\Filterpad\To_SeaBASS\';
            pigsheetpath='\Mar2010\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC5';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_mar2010.xlsx'],...
                'Sheet',wksheet);
           file_pref='GulfCarbon_gulfcarbon5_filterpad_';
    end
    
    % Load filterpad data
    infile=[padpath,'all_pad.mat'];
    load(infile);
    
    %Group stations into groups with fmicro > 0.5 and fpico > 0.5, and all else mixed
    fmicro_sta=stn_pig_sf.Station(stn_pig_sf.fmicro>0.5 & ~contains(stn_pig_sf.Station,'mr'));  % Omit river stations
    fpico_sta=stn_pig_sf.Station(stn_pig_sf.fpico>0.5 & ~contains(stn_pig_sf.Station,'mr'));
    fpiconano_sta=stn_pig_sf.Station(stn_pig_sf.fpico + stn_pig_sf.fnano>=0.5 & ~contains(stn_pig_sf.Station,'mr'));
    
    % Determine indices for stations for different phytoplankton size 
    %   class dominance (surface samples only)
    fmicro_indx=[];
    fpico_indx=[];
    fpiconano_indx=[];
    for ik=1:length(all_pad.station)
        if str2num(all_pad.depth{ik})<=3 && any(contains(fmicro_sta,all_pad.station{ik},'Ignorecase',true))
            fmicro_indx=[fmicro_indx;ik];
        elseif str2num(all_pad.depth{ik})<=3 && any(contains(fpiconano_sta,all_pad.station{ik},'Ignorecase',true))
            %fpiconano_test=find(contains(fpiconano_sta,all_pad.station{ik},'Ignorecase',true));
            fpiconano_indx=[fpiconano_indx;ik];
        end
        
        if str2num(all_pad.depth{ik})<=3 && any(contains(fpico_sta,all_pad.station{ik},'Ignorecase',true))
            %fpico_test=find(contains(fpico_sta,all_pad.station{ik},'Ignorecase',true));
            fpico_indx=[fpico_indx;ik];
        end
    end
    
    % % Test table to see what stations and depths are included
    % ftab = table(cell2mat(all_pad.station([fmicro_indx;fpiconano_indx])),cell2mat(all_pad.depth([fmicro_indx;fpiconano_indx])));
    
    pad_hdr=all_pad.fields{1};
    pad_units=all_pad.units{1};
    
    % Select variable to plot
    var={'aph'};  %Select either "ap", "ad", "aph", "ad_model","aph_model"
    var_str=var;
    var_char_num=length(var);
    
    %Find variable and lambda indices
    [m,n]=size(var);
    var_indx=zeros(n,1);
    for ip=1:n
        var_indx(ip)=find(strcmp(pad_hdr,var(ip))==1);
        lambda_indx(ip)=find(strcmp(pad_hdr,'wavelength'));
    end
    
    max_dat=0;
    min_dat=100;
    
    %Loop for each variable 
    for i=1:n
    
        %Loop for each size class
        for isf=1:3
            %Initialize variables
            pad_dat=cell(m,1);
            pad_depth=cell(m,1);
            pad_station=cell(m,1);
            pad_lambda=cell(m,1);
            
            switch isf
                case 1
                    stn_num=length(fmicro_indx);
                    indx1=fmicro_indx;
                case 2
                    stn_num=length(fpico_indx);
                    indx1=fpico_indx;
                case 3
                    stn_num=length(fpiconano_indx);
                    indx1=fpiconano_indx;
            end
            
            %Initialize variables
            lambdan = zeros(n,stn_num);
            pad_dat=cell(n,stn_num);
            pad_depth=cell(n,stn_num);
            pad_station=cell(n,stn_num);
            pad_lambda=cell(n,stn_num);
    
            %Loop through stations
            for stn_indx=1:stn_num
                % indx1=find(datenum(all_pad(stn).date)>datenum('12/31/08') & contains(all_pad.station,sta_txt)); 
    
                j=indx1(stn_indx); % Indices for each size class
            
                % Omit spectra with anomalous values
                if any(all_pad.data{j}(:,var_indx(i))<-0.1)
                    continue
                end
                
                pad_depth(i,stn_indx)={all_pad.depth{j}(:,1)};
                lambda_all={all_pad.data{j}(:,lambda_indx(i))};

                % Get lambda indices for 350 - 700 nm
                lambda_range=find(lambda_all{1}>=350 & lambda_all{1}<=750);
                pad_lambda(i,stn_indx)={lambda_all{1}(lambda_range)};
                pad_dat(i,stn_indx)={all_pad.data{j}(lambda_range,var_indx(i))};
                pad_station(i,stn_indx)=strrep(strrep(all_pad.name(j),file_pref,''),'_1m_R2.dat','');
                lambdan(i,stn_indx) = length(lambda_all{1}(lambda_range));
                
                %{
                % Find max for plotting
                new_max_dat=max(pad_dat{i,stn_indx},[],1);
                new_min_dat=min(pad_dat{i,stn_indx},[],1);
                if new_max_dat > max_dat
                    max_dat=new_max_dat;
                end
                if new_min_dat < min_dat
                    min_dat=new_min_dat;
                end
                %}
            end
    
            % Convert cell arrays to matrices
    
            % Determine most common number of wavelengths from
            %    spectrophotometer scans and use only those records
    
            switch isf
                case 1
                    mode_indx = find(lambdan==mode(lambdan));
                    pad_lambda_fmicro=cell2mat(pad_lambda(mode_indx));
                    pad_dat_fmicro=cell2mat(pad_dat(mode_indx));
                    pad_station_fmicro=pad_station(mode_indx);
                    pad_depth_fmicro=pad_depth(mode_indx);
                case 2
                    mode_indx = find(lambdan==mode(lambdan));
                    pad_lambda_fpico=cell2mat(pad_lambda(mode_indx));
                    pad_dat_fpico=cell2mat(pad_dat(mode_indx));
                    pad_station_fpico=pad_station(mode_indx);
                    pad_depth_fpico=pad_depth(mode_indx);
                case 3
                    mode_indx = find(lambdan==mode(lambdan));
                    pad_lambda_fpiconano=cell2mat(pad_lambda(mode_indx));
                    pad_dat_fpiconano=cell2mat(pad_dat(mode_indx));
                    pad_station_fpiconano=pad_station(mode_indx);
                    pad_depth_fpiconano=pad_depth(mode_indx);
            end
        end
    end
    
    % Override selection of min and max if desired
    min_dat=-0.02;
    max_dat=1.1;
    %
    hplot=[];
    figure(1);
    clf
    for k=1:n
        for isf=[1,3]
            switch isf
                case 1
                    pad_lambda=pad_lambda_fmicro;
                    pad_dat=pad_dat_fmicro;
                    plt_clr='#0072BD';
                case 2
                    pad_lambda=pad_lambda_fpico;
                    pad_dat=pad_dat_fpico;
                    plt_clr='#EDB120';
                case 3
                    pad_lambda=pad_lambda_fpiconano;
                    pad_dat=pad_dat_fpiconano;
                    %plt_clr='#D95319';
                    plt_clr='#EDB120';
            end
            
            normalize = 1;
            
            switch normalize
                case 1
                    % Normalized to aph440
                    lwidth = 2.5;
                    hplot(isf)=plot(pad_lambda(:,1),mean(pad_dat./pad_dat(idx_440(1),:),2,'omitnan'),'-',...
                       'Linewidth',lwidth,'Color',plt_clr);
                    hold on
                    
                    % Turn off X and Y axis labels as necessary
                    % set(gca,'YTickLabels',{}); %'XTickLabels',{},
                    
                    % Plot errorbars
                    herr(isf)=errorbar(pad_lambda(:,1),mean(pad_dat./pad_dat(idx_440(1),:),2,'omitnan'),...
                        std(pad_dat./pad_dat(idx_440,:),[],2,'omitmissing'),'-', ...
                       'Linewidth',lwidth,'Color',plt_clr,'capsize',0);
    
                    % Set transparency level (0:1)
                    alpha = 0.05;   
                    % Set transparency (undocumented)
                    set([herr(isf).Bar, herr(isf).Line], 'ColorType', 'truecoloralpha', 'ColorData', [herr(isf).Line.ColorData(1:3); 255*alpha]);               
                    
                    % Plot individual spectra if desired
                    % hpad=plot(pad_lambda(:,1),pad_dat,'-','Linewidth',0.5,'Color',[0.6,0.6,0.6]); %Show individual spectra
                    
                    % Plot error lines if desired
                    % hstd1(isf)=plot(pad_lambda(:,1),...
                    %     mean(pad_dat./pad_dat(idx_440,:),2)+std(pad_dat./pad_dat(idx_440,:),[],2,'omitmissing'),...
                    %     ':','Linewidth',lwidth./2,'Color',plt_clr);
                    % hstd2(isf)=plot(pad_lambda(:,1),...
                    %     mean(pad_dat./pad_dat(idx_440,:),2)-std(pad_dat./pad_dat(idx_440,:),[],2,'omitmissing'),...
                    %     ':','Linewidth',lwidth./2,'Color',plt_clr);
                    % 
                    %Create y label
                    ylab = '$\bf{\hat{\it{a}}_{\it{ph}}}$';
                    
                  case 0
                    % No normalization
                    % 
                    max_dat = max(mean(pad_dat,2));
            
                    hplot(isf)=plot(pad_lambda(:,1),mean(pad_dat,'omitnan'),'-','Linewidth',2.5,'Color',plt_clr);
                    hold on
                    hstd1(isf)=plot(pad_lambda(:,1),...
                        (mean(pad_dat,2)+std(pad_dat,[],2))./mean(pad_dat,2,'omitmissing'),...
                        ':','Linewidth',0.7,'Color',plt_clr);
                    hstd2(isf)=plot(pad_lambda(:,1),...
                        (mean(pad_dat,2)-std(pad_dat,[],2))./mean(pad_dat,2,'omitmissing'),...
                        ':','Linewidth',0.7,'Color',plt_clr);
                    hold on
                    % if isf==2
                        % hpad=plot(pad_lambda(:,1),pad_dat,'-','Linewidth',0.5,'Color',[0.6,0.6,0.6]); %Show individual spectra
                    % end
                    
                    %Create y label
                    ylab='\it{a_{ph}}';
            end
    
            if k==1 && isf==1
                set(gca,'Linewidth',1.3,'fontname','times','fontweight','bold','fontsize',18,...
                'xlim',[375 750],'xtick',400:100:800,'ylim',[min_dat-0.01.*abs(min_dat) max_dat+0.05.*max_dat]);
                xlabel('Wavelength (nm)','fontname','times','fontweight','bold','fontsize',18);
                ylabel(ylab,'Interpreter','latex','fontname','arial','fontweight','bold','fontsize',18);
            end
    
            if isf==1
                pad_dat_fmicro=pad_dat;
                save([pigsheetpath,'..\','fmicro_aph_dat.mat'],'pad_dat_fmicro','pad_lambda_fmicro');
            elseif isf==2
                pad_dat_fpico=pad_dat;
                save([pigsheetpath,'..\','fpico_aph_dat.mat'],'pad_dat_fpico','pad_lambda_fpico');
            else 
                pad_dat_fpiconano=pad_dat;
                save([pigsheetpath,'..\','fpiconano_aph_dat.mat'],'pad_dat_fpiconano','pad_lambda_fpiconano');
            end
        end
    end
    
    hleg=legend(hplot([1,3]),'fmicro','fpico/nano');
    set(hleg,'Fontsize',10);
    
    % Add annotation if desired
    % text(500,(max_dat-min_dat)-(max_dat-min_dat).*.1,['Station ',strrep(sta_txt,'_',' ')],'fontname','times','fontweight','bold','fontsize',[18]);
    % text(375,.2,'Shelf','fontname','times','fontweight','bold','fontsize',[18]);
    hold on
    
    % Adjust figure position
    fig_pos=get(gca,'position');
    fig_pos(1:2)=fig_pos(1:2)+0.025;
    set(gca,'Position',fig_pos);
    
    print([pigsheetpath,'..\',wksheet(end-2:end),'_pad_aph_shape_vectors.tif'],'-dtiff','-r600');
end

disp('Completed');
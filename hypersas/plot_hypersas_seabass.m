% plot_hypersas_seabass.m - Matlab script to plot HyperSAS spectra
% 
% Syntax: plot_hypersas_seabass.m
%
% Inputs:
%    1) Folder location with '*.mat' files produced by extract_hypersas_seabass.m
%    2) Spreadsheet with plot parameters
%    ('apr2009_hypersas_plot_qaa_parameters_v2.xlsx')
%
% Outputs:
%    -Two figures are produced, one with plot of HyperSAS spectra of Rrs 
%      (or other selected variable) and the other a time series of HyperSAS Rrs
%    -A print statement current produces a '*.tif' file of the Rrs time
%      series
%   
% Other m-files required: 
%  1) The '*.mat' files with the extracted SeaBASS data were generated with 
%    'extract_hypersas_seabass.m' 

% MAT-files required: see above
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 10 Aug 2021

%------------- BEGIN CODE --------------%

clearvars
clc
close all
clf('reset')

file_pref='GulfCarbon2_gulfcarbon2_hypersas_St';

%Specify information for input path and file
inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\To_SeaBASS\';
infile=[inpath,'all_hyper_seabass_April2009.mat'];

rrs_txt_output_flag=0;  %Set to 1 for txt output for pixel matchup file 

load(infile);

%Spreadsheet with input parameters
xls_path='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\';
infile_name=[xls_path,'apr2009_hypersas_plot_qaa_parameters_v2.xlsx'];

%Input parameters
[~,~,all_par_3]=xlsread(infile_name,1,'C3:R77');
hyper_sta=all_par_3(:,1);
hyper_doy=cell2mat(all_par_3(:,2));
hyper_date=cell2mat(all_par_3(:,16));
hyper_run=cell2mat(all_par_3(:,3));
hyper_start_time=cell2mat(all_par_3(:,4));
hyper_end_time=cell2mat(all_par_3(:,5));
sta_lab=all_par_3(:,12);
title_lab=all_par_3(:,13);
sas_basecorr_flag=cell2mat(all_par_3(:,14));

%hyper_date=datestr(datenum(hyper_date),'yyyymmdd');

single_sta='Run52';  %Single station for processing just one
num_char=length(file_pref)+length(single_sta);

run_indx=[1:4,6:14,16,18:23,25:28,30:31,33:46,48:49,52:56,58:61,63];
run_array=cell(1,length(run_indx));
for runn=run_indx
    run_array(runn)={['Run',num2str(runn)]};
end

%For all stations and runs
station_names=cat(2,{'A1','A1b','A2','A3','B4','C1','C5','E0','E1','E3','D5','F4','F5','MR2','MR3',...
'G2','G3','G4','G5','H1','H2','H3','UWS1','UWS2'},run_array);

[nsta,~]=size(hyper_sta);

%hfig=figure(1);

%Loop through stations and runs for plotting
for im=1:nsta 

    sta_txt=hyper_sta{im};  %Select station to plot

    %To do all stations, comment out the continue statement below
    if strcmp(single_sta,sta_txt)~=1
        continue
    end

    %HyperSAS parameters
    file_no=hyper_doy(im);
    run_no=hyper_run(im);
    file_date=hyper_date(im,:);
    time_str1=datestr(hyper_start_time(im),'HH:MM:SS');
    time_str2=datestr(hyper_end_time(im),'HH:MM:SS');
    indx1=find(datenum(all_hyper.start_date)>datenum('12/31/08') & ...
        contains(all_hyper.name,sta_txt)); 
    [m,n]=size(indx1);

    hyper_hdr=all_hyper.fields{im};
    hyper_units=all_hyper.units{im};

    %Variable options for spectral plot are 'Lt','Lsky','Lw','Es' and 'Rrs'
    var='Rrs';  
    var_str=var;
    var_char_num=length(var);

    %Corresponding ylabel options are 'L_t (\muW cm^{-2} nm^{-1} sr^{-1}','L_{sky} (\muW cm^{-2} nm^{-1} sr^{-1}',
    % 'L_w (\muW cm^{-2} nm^{-1} sr^{-1}','E_s (\muW cm^{-2} nm^{-1}', and 'R_{RS} (sr^{-1})'
    switch var
        case 'Lt'
            ylabel_txt='L_t (\muW cm^{-2} nm^{-1} sr^{-1}';
        case 'Lsky'
            ylabel_txt='L_{sky} (\muW cm^{-2} nm^{-1} sr^{-1}';
        case 'Lw'
            ylabel_txt='L_w (\muW cm^{-2} nm^{-1} sr^{-1}';
        case 'Es'
            ylabel_txt='E_s (\muW cm^{-2} nm^{-1}';
        case 'Rrs'
            ylabel_txt='R_{rs} (sr^{-1}';
    end

    var_indx=find(strncmp(hyper_hdr,var,var_char_num)==1 & ~contains(hyper_hdr,'_1id')); 
    hyper_lambda=str2double(strrep(hyper_hdr(var_indx),var,''));

    dvar_indx=find(strncmp(hyper_hdr,var,var_char_num)==1 & contains(hyper_hdr,'_1id')); 
    discrete_lambda=str2double(strrep(strrep(hyper_hdr(dvar_indx),var,''),'_1id',''));
    
    lat_indx=find(strcmp(hyper_hdr,'lat')==1);
    lon_indx=find(strcmp(hyper_hdr,'lon')==1);
    
    hyper_dat=cell(m,1);
    discrete_dat=cell(m,1);
    hyper_time=cell(m,1);
    hyper_depth=cell(m,1);
    hyper_station=cell(m,1);
    hyper_lat=cell(m,1);
    hyper_lon=cell(m,1);

    max_dat=0;
    min_dat=100;

    for ih=1:m
        j=indx1(ih);

        hyper_dat{ih}=all_hyper.data{j}(:,var_indx-2); %Adjust var_indx to account
        % for date and time fields not included in data array
        discrete_dat{ih}=all_hyper.data{j}(:,dvar_indx-2); %As above
                
        %Replace missing data tag (-9999) with NaN
        hyper_dat{ih}(hyper_dat{ih}==-9999)=NaN;
        discrete_dat{ih}(discrete_dat{ih}==-9999)=NaN;
        
        hyper_time(ih)=all_hyper.run_time(j);
        hyper_station(ih)=all_hyper.station(j);
        hyper_lat{ih}=all_hyper.data{j}(:,lat_indx-2);
        hyper_lon{ih}=all_hyper.data{j}(:,lon_indx-2);

        if rrs_txt_output_flag==1
            %Generate text file for pixel matchups
            outtextfile=[inpath,'..\','rrs_Apr2009_',num2str(file_no),'_',num2str(run_no),'_',sta_txt,'.txt'];

            lambda_str=[];
            for lm=1:136
                lambda_str=[lambda_str,num2str(hyper_lambda(lm),'%5.1f'),'\t'];
            end
            lambda_str=[lambda_str,num2str(hyper_lambda(end),'%5.1f'),'\n'];
            file_hdr_format=['Seq \t UTC Time \t Lat \t Lon \t',lambda_str];
            file_format_spec=['%u\t%s',repmat('\t%f',1,139'),'\n'];
            hdrtext='#Matchup input sequence';
            fileID=fopen(outtextfile,'w');
                fprintf(fileID,'%s\n',hdrtext);
                fprintf(fileID,file_hdr_format);
                for dm=1:size(hyper_time,1)
                    fprintf(fileID,file_format_spec,dm,hyper_time(dm,:),hyper_lat(dm),-hyper_lon(dm),hyper_dat(dm,:));
                end
            fclose(fileID);

        end
        
        %Find max for plotting
        new_max_dat=max(max(all_hyper.data{j}(:,var_indx-2),[],1));
        new_min_dat=min(min(all_hyper.data{j}(:,var_indx-2),[],1));
        if new_max_dat > max_dat
            max_dat=new_max_dat;
        end
        if new_min_dat < min_dat
            min_dat=new_min_dat;
        end
    end

    %Prepare figure for spectral plot
    hfig1=figure(1);
    %Set up axes
    hrrs=axes('YLimMode','auto','XTick',300:100:800,'Xlim',[290 810],'Position',[0.25 0.2 0.6 0.7]);

    xtick_lab=get(gca,'Xtick');
    set(hrrs,'Visible','on','Color','none','YAxisLocation','left','FontSize',14,...
        'FontWeight','normal','FontName','Arial','Linewidth',1.5,'XtickLabel',xtick_lab); 
    box on
    hold on

    %Axis labels
    hylab=ylabel(ylabel_txt,'FontSize',14,'FontWeight','bold','FontName','Arial',...
        'Rotation',90);
    xlabel('Wavelength(nm)','FontSize',14,'FontWeight','bold','FontName','Arial');
    hold on

    %Override selection of min and max if desired
    %min_dat=0;
    %max_dat=1;

    for k=1:m
        %Plot hypersas portion
        hplot=plot(hyper_lambda,hyper_dat{k},'-','Linewidth',1.2);
        set(hplot,'Linewidth',1.5);
        hold on
 
        %Plot discrete if desiged
        hplot2=plot(discrete_lambda,discrete_dat{k},'o','Linewidth',1.2);
        set(hplot2,'Linewidth',1.5);
        hold on
    end

    ylimits=get(gca,'YLim');
    if contains(sta_txt,'UWS') || contains(sta_txt,'Run')
        sta_label=strrep(sta_txt,'_',' ');
    else
        sta_label=['Station ',strrep(sta_txt,'_',' ')];
    end        

    text(500,ylimits(2)+0.05.*ylimits(2),sta_label,'fontname','times','fontweight','bold','fontsize',18);
    hold on
    fig_pos=get(gca,'position');
    fig_pos(2)=fig_pos(2)-0.0;
    set(gca,'Position',fig_pos);

    %Prepare figure for time series plot of Rrs
    hfig2=figure(2);
    scrsz = get(groot,'ScreenSize');
    set(hfig2,'Position',[scrsz(4).*.1 scrsz(3).*.1 scrsz(3).*.75 scrsz(4).*.55])
    set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

    %Set up axes
    hrrs2=axes('YLimMode','auto','Color','none','FontSize',10,...
        'FontWeight','normal','FontName','Arial','Linewidth',1.75,...
        'FontSize',12,'YAxisLocation','left'); %,'Ylim',[0 10]);
    %box on
    hold on
    
    %Axis labels
    ylabel('Rrs (sr^{-1})','FontSize',14,'FontWeight','bold','FontName','Arial');
    xlabel('Time (UTC)','FontSize',14,'FontWeight','bold','FontName','Arial');
    hold on

    for k=1:m
        hplot=plot(datenum(hyper_time{k}),hyper_dat{k}(:,floor(hyper_lambda)==553),'o','Linewidth',1.2);
        hold on
    end

    xtickdat=get(gca,'Xtick');
    xlimits=get(gca,'Xlim');
    xtick_label=datestr(xtickdat,'HH:MM');
    set(gca,'XtickLabel',xtick_label);
    hold on
    
    %Secondary axis for discrete plot
    ax_pos=get(hrrs2,'Position');
    hrrs3=axes('YLimMode','auto','XTick',xtickdat,'Xlim',xlimits,'Color','none','FontSize',10,...
        'FontWeight','normal','FontName','Arial','Linewidth',1.75,...
        'FontSize',12,'Position',ax_pos,'YAxisLocation','right',...
        'XAxisLocation','top','XTickLabel',[]); %,'Ylim',[0 10]);
    hold on
    
    for k=1:m
        hdplot=plot(datenum(hyper_time{k}),discrete_dat{k}(:,4),'m+','Linewidth',1.2,'Parent',hrrs3);
        hold on
    end
    
    %Axis labels
    hylab=ylabel('Discrete Rrs (sr^{-1})','FontSize',14,'FontWeight','bold','FontName','Arial',...
        'Rotation',270);
    ylab_pos=get(hylab,'Position');
    ylab_pos(1)=ylab_pos(1)+0.001;
    set(hylab,'Position',ylab_pos);

    hleg=legend([hplot,hdplot],'Rrs553 (sr^{-1})','Discrete Rrs380 (sr^{-1})','Location','southwest');
    set(hleg,'Fontsize',12);

    %Override selection of min and max if desired
    %min_dat=0;
    %max_dat=1;

    %hleg=legend(hplot,strrep(var_str,'_',''));
    %set(hleg,'Fontsize',10);

    ylimits=get(gca,'YLim');
    text(xlimits(1)+0.4.*(xlimits(2)-xlimits(1)),ylimits(2)+0.05.*(ylimits(2)-ylimits(1)),['Station ',strrep(sta_txt,'_',' ')],'fontname','times','fontweight','bold','fontsize',[18]);
    hold on
    fig_pos2=get(gca,'position');
    fig_pos2(2)=fig_pos2(2)-0.0;
    set(gca,'Position',fig_pos2);

    print([inpath,'..\hypersas_rrs_time_series_GC2_',sta_txt,'.tif'],'-dtiff','-r600');
end

disp('Completed');

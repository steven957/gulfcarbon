% plot_hypersas_hyperpro_rrs_Apr2009_wqaa_output_v8c.m - Program 
%   to load hypersas data from compiled '*.mat' files produced by 
%   'generate_hypersas_reflectance_radiance_files_all_cruise_rev6.m'
%
% Syntax:  plot_hypersas_hyperpro_rrs_Apr2009_wqaa_output_v8c.m
%
% Inputs:
%    1) Folder location with compiled data files ('*.mat')
%      produced by 'generate_hypersas_reflectance_radiance_files_all_cruise_rev6.m'
%    2) Spreadsheet with station metadata information 
%      ('apr2009_hypersas_plot_qaa_parameters_v2a.xls')
%
% Outputs:
%    1) '*.mat' file for quasi-analytical algorithm results for retrieved 
%      IOP's
%    2) plot of reflectance spectra comparison among HyperSAS, HyperPro,
%      and satellite observations
%   
% Other m-files required: 
%  1) The script 'hypersas_data_read.m' is required to read the
%      HyperSAS data files and load variables for processing
%  2) The script 'hypersas_data_filter.m' is required to apply filter and 
%      baseline corrections to HyperSAS reflectance
%  3) The function 'mqaa_v6.m' is required to retrieve IOPs from the
%      reflectance data
%  4) [OPTIONAL] The script 'hyperpro_data_read_and_filter.m' is used to
%      read and filter the HyperPro processed and compiled data if desired
%      (set hproplot=0 to omit this step)
%
% MAT-files required: 
%  1) Compiled data files ('*.mat') produced by 
%     'generate_hypersas_reflectance_radiance_files_all_cruise_rev6.m'
%  2) [OPTIONAL] Processed and compiled HyperPro data ('*.mat') produced by
%     'hyperpro_matfile_compile_GC2.m' and 'hyperpro_process_GC2_v4.m'
%  3) [OPTIONAL] Compiled satellite matchup data ('*.mat') from MODIS Aqua produced
%     by 'pixel_matchup_master_v2.m' and function
%     'pixel_matchup_readdat_v3.m' (set satplot=0 to omit this step)
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 6 Sep 2021

%% ------------- BEGIN CODE --------------%

clear
clc
close all

%Load file locations and input parameters
inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\';

%Spreadsheet with input parameters
infile_name=[inpath,'apr2009_hypersas_plot_qaa_parameters_v2a.xlsx'];

single_sta='G2';  %Single station for processing just one

sas_omit_index=[]; %Initialize omit_index variable
pro_omit_index=[]; %Initialize omit_index variable 
run_qaa=1;   %1 for true, 0 for false - runs QAA routine on HyperSAS and HyperPro
chl_alg='OC3M'; %'OC3M'  'Gittelson'

%Read input parameters from spreadsheet
[~,~,all_par_3]=xlsread(infile_name,1,'C3:Q77');
hyper_sta=all_par_3(:,1);
hyper_doy=cell2mat(all_par_3(:,2));
hyper_run=cell2mat(all_par_3(:,3));
hyper_start_time=cell2mat(all_par_3(:,4));
hyper_end_time=cell2mat(all_par_3(:,5));
pro_group=cell2mat(all_par_3(:,6));
sat_test_str=all_par_3(:,8);
sat_plot_flag=cell2mat(all_par_3(:,9));
pro_plot_flag=cell2mat(all_par_3(:,10));
pro_scan_omit=all_par_3(:,11);
sta_lab=all_par_3(:,12);
title_lab=all_par_3(:,13);
sas_basecorr_flag=cell2mat(all_par_3(:,14));
pro_basecorr_flag=cell2mat(all_par_3(:,15));

%Specify which runs to include (omit runs with no useful data)
run_indx=[1:4,6:14,16,18:23,25:28,30:31,33:46,48:49,52:56,58:61,63];
run_array=cell(1,length(run_indx));
for runn=run_indx
    run_array(runn)={['Run',num2str(runn)]};
end

%For all stations and runs
station_names=cat(2,{'A1','A1b','A2','A3','B4','C1','C5','E0','E1','E3','D5','F4','F5','MR2','MR3',...
'G2','G3','G4','G5','H1','H2','H3','UWS1','UWS2'},run_array);

%For stations only
%station_names={'A1','A1b','A2','A3','B4','C1','C5','E0','E1','E3','D5','F4','F5','MR2','MR3',...
%'G2','G3','G4','G5','H1','H2','H3','UWS1','UWS2'};

[nsta,~]=size(hyper_sta);

hfig=figure(1);

%% Begin loop to process each station and run 

for im=1:nsta 

    sta_txt=hyper_sta{im};  %Select station to plot

    %To do all stations, comment out the continue statement; otherwise only
    % station given in single_sta will be processed
    if strcmp(single_sta,sta_txt)~=1
        continue
    end

    %HyperSAS parameters
    file_no=hyper_doy(im);
    run_no=hyper_run(im);
    time_str1=datestr(hyper_start_time(im),'HH:MM:SS');
    time_str2=datestr(hyper_end_time(im),'HH:MM:SS');

    %HyperPro group (set of profiles at a given station)
    group_no=pro_group(im);

    %Satellite time string
    sat_test_string=sat_test_str{im};

    %Flags for plotting
    satplot=sat_plot_flag(im);
    hproplot=pro_plot_flag(im);

    station_txt=sta_lab{im};
    title_txt=title_lab{im};

    %Baseline correction flags
    base_corr=sas_basecorr_flag(im);    %HyperSAS baseline correction (1=on; 0=off)
    base_corr2=pro_basecorr_flag(im);   %HyperPro baseline correction (1=on; 0=off)

    %Omit HyperPro individual profiles if needed
    omit_str = sprintf('%s,', pro_scan_omit{im});
    pro_scan_val=sscanf(omit_str, '%g,',inf);

    if ~isnan(pro_scan_val)
        pro_omit_index=pro_scan_val;
    else
        pro_omit_index=[];
    end

    %file_name=[inpath,'rad_rrs_',num2str(file_no),'_',num2str(run_no),'.mat'];
    file_name=[inpath,'rad_rrs_Es253_',num2str(file_no),'.mat'];
    text_msg=['loading ',file_name,'...'];

    %% Call Matlab scripts to load file and variables and apply filters
    
    %Run matlab script to load file and populate variables
    hypersas_data_read;

    %Run matlab script to apply data filters and baseline corrections
    hypersas_data_filter;

    %Move on to next station if no acceptable data
    if isempty(good_lt_index)==1 || isempty(good_dlt_index)==1
        continue
    end

    mean_rrs=mean(rrs_corr,1);
    std_rrs=std(rrs_corr,0,1);

    mean_drrs=mean(drrs_corr,1,'omitnan');
    std_drrs=std(drrs_corr,0,1,'omitnan');

    %% Run QAA algorithm to retrieve IOPs
    
    %Run QAA algorithm
    if run_qaa==1
        [m,n]=size(rrs_corr);

        %Initialize output matrices
        a=zeros(m,n);
        adg=zeros(m,n);
        aph=zeros(m,n);
        bb=zeros(m,n);
        rrs_u=zeros(m,n);

        for k=1:m
            %Call qaa algorithm
            [nbands, wavel, rrs_u(k,:), aw, a(k,:), adg(k,:), aph(k,:), bbw, bb(k,:)] = ...
                mqaa_v6(rrs_corr(k,:),lt_lambda');
        end
        hyperqaa_outfile=[inpath,'qaa_v6_sas_output_',sta_txt,'.mat'];
        save(hyperqaa_outfile,'wavel','rrs_u','aw','a','adg','aph','bbw','bb',...
            'lt_time_filt','rrs_lat_filt','rrs_lon_filt');
    end

    %% Estimate chlorophyll and Kd490 using specified algorithms 

    %Find wavelength index for Rrs(443), Rrs(488), Rrs(547), and Rrs(750-800)
    index_443=find(lt_lambda>441 & lt_lambda<445);
    index_488=find(lt_lambda>487 & lt_lambda<490);
    index_547=find(lt_lambda>545 & lt_lambda<549);
    nir_index=find(lt_lambda>750 & lt_lambda<800);
    index_684=find(lt_lambda>683 & lt_lambda<685);
    index_700=find(lt_lambda>699 & lt_lambda<701);
    index_720=find(lt_lambda>719 & lt_lambda<721);

    rrs_blue=max([rrs_corr(:,index_443),rrs_corr(:,index_488)],[],2); %-rrs_nir;
    rrs_490=rrs_corr(:,index_488); 
    rrs_green=rrs_corr(:,index_547); %-rrs_nir;

    switch chl_alg
        case 'OC3M'
            %Determine chlor_a using OC3M algorithm:
            a0=0.2424;
            a_coef(1)=-2.7423;
            a_coef(2)=1.8017;
            a_coef(3)=0.0015;
            a_coef(4)=-1.2280;

            chlor_a_sum_term=0;
            for a_indx=1:4
                chlor_a_sum_term=chlor_a_sum_term+a_coef(a_indx).*log10(rrs_blue./rrs_green).^a_indx;
            end

            log10_chlor_a=a0+chlor_a_sum_term;
            chlor_a=10.^log10_chlor_a;

        case 'Gitelson'
            %Compute chla using Gitelson et al. (2011) turbid water chl algorithm
            %chl-a = 418.88{[Rrs(684)−1−Rrs(700)−1]Rrs(720)}+19.275

            chlor_a=418.88.*((1./rrs_corr(:,index_684)-1./rrs_corr(:,index_700)).*rrs_corr(:,index_720))+19.275;
    end

    %Compute Kd490 using NASA algorithm
    a0_kd(1)=-0.8813;
    a_coef_kd(1)=-2.0584;
    a_coef_kd(2)=2.5878;
    a_coef_kd(3)=-3.4885;
    a_coef_kd(4)=-1.5061;

    kd_sum_term=0;
    for a_indx=1:4
        kd_sum_term=kd_sum_term+a_coef_kd(a_indx).*log10(rrs_490./rrs_green).^a_indx;
    end

    log10_kd=a0_kd+kd_sum_term;
    kd490=10.^log10_kd + 0.0166;

    %% Plot HyperSAS results
    
    if contains(sta_txt,'Run')
        %Plot transect map of chl if desired
        %plot_HyperSAS_chl_run_transect

        %Plot transect map of Kd490 if desired
        %plot_HyperSAS_kd490_run_transect
    end

    %Return to figure(1)
    figure(hfig);
    clf reset

    %Set up axes
    hrrs=axes(hfig,'YLimMode','auto','XTick',350:100:800,'Xlim',[290,810],'Position',[0.25,0.2,0.6,0.7]);

    xtick_lab=get(gca,'Xtick');
    set(hrrs,'Visible','on','Color','none','YAxisLocation','left','FontSize',20,...
        'FontWeight','normal','FontName','Arial','Linewidth',1.5,'XtickLabel',xtick_lab); 
    box on
    hold on

    %Plot rrs
    rrscolor='#85C0F9';

    %hhyper=plot(lt_lambda,rrs_corr,'-m','Linewidth',1.2);
    hhyper=plot(lt_lambda(1:end),mean_rrs(1:end),'-','Linewidth',2,'Color',rrscolor);
    %set(gca,'YLimMode','auto');
    %set(gca,'YLim',[-0.0002,0.004]);
    hold on

    hstd1=plot(lt_lambda(1:end),mean_rrs(1:end)+std_rrs(1:end),':','Linewidth',1.5,'Color',rrscolor);
    hstd2=plot(lt_lambda(1:end),mean_rrs(1:end)-std_rrs(1:end),':','Linewidth',1.5,'Color',rrscolor);
    hold on

    %Plot individual spectra if desired
    %plot(lt_lambda,rrs_corr,'-m');

    %Plot discrete UV rrs
    discolor='#A95AA1';
    %hdisc=plot(drrs_lambda(1:4),drrs_corr,'om','Linewidth',1.2);
    hdisc=plot(drrs_lambda(1:4),mean_drrs,'o','Linewidth',2,'Color',discolor);
    hold on

    hstdd1=plot(drrs_lambda(1:4),mean_drrs+std_drrs,':','Linewidth',1.5,'Color',discolor);
    hstdd2=plot(drrs_lambda(1:4),mean_drrs-std_drrs,':','Linewidth',1.5,'Color',discolor);
    hold on

    %Plot individual spectra if desired
    %plot(drrs_lambda(1:4),drrs_corr,'+m');

    %Axis labels
    hylab=ylabel('{\itR_{rs}} (sr^{-1})','FontSize',20,'FontWeight','bold','FontName','Arial',...
        'Rotation',90);
    %ylab_pos=get(hylab,'Position');
    %ylab_pos=[ylab_pos(1)+0.1.*time_range ylab_pos(2) ylab_pos(3)];
    %set(hylab,'Position',ylab_pos);
    xlabel('Wavelength(nm)','FontSize',20,'FontWeight','bold','FontName','Arial');
    hold on

    %Format y-axis label with exponential notation
    hrrs.YAxis.Exponent = -2;

    sas_ylimits=get(gca,'YLim');

    if satplot==0 && hproplot==0

        legend_str={'HyperSAS','Discrete'};

        hleg=legend([hhyper(1),hdisc(1)],legend_str);
        set(hleg,'Fontsize',14,'FontName','Arial','Location','Best');

        set(gca,'Ylim',sas_ylimits);
        hold on
    end      

    %% Read filter and plot HyperPro results 

    if hproplot==1
        %Load data 
        hyperpro_data_read_and_filter %This calls Matlab script to read and filter HyperPro data

        %Run QAA algorithm if desired
        if run_qaa==1
        [mm,nn]=size(rrs2);

            %Initialize output matrices
            a=zeros(mm,nn);
            adg=zeros(mm,nn);
            aph=zeros(mm,nn);
            bb=zeros(mm,nn);
            rrs_u=zeros(mm,nn);

            for k=1:mm
                %Call qaa algorithm
                [nbands, wavel, rrs_u(k,:), aw, a(k,:), adg(k,:), aph(k,:), bbw, bb(k,:)] = mqaa_v6(rrs2(k,:),rrs2_lambda');
            end
            save([inpath,'qaa_pro_output_sta_',station_txt,'.mat'],'wavel','rrs_u','aw','a','adg','aph','bbw','bb');
        end

        procolor='#F5973A';
        
        %Plot spectra
        hhyper2=plot(rrs2_lambda(rrs2_lambda_index),mean_rrs2(rrs2_lambda_index),'-','Linewidth',2,'Color',procolor);
        hold on

        hstd1_2=plot(rrs2_lambda(rrs2_lambda_index),mean_rrs2(rrs2_lambda_index)+std_rrs2(rrs2_lambda_index),':','Linewidth',1.5,'Color',procolor);
        hstd2_2=plot(rrs2_lambda(rrs2_lambda_index),mean_rrs2(rrs2_lambda_index)-std_rrs2(rrs2_lambda_index),':','Linewidth',1.5,'Color',procolor);
        hold on

        %Plot individual spectra if desired
        %hhyper2=plot(rrs2_lambda,rrs2,'-b');

        if satplot==0

            legend_str={'HyperSAS','HyperPro','Discrete'};

            hleg=legend([hhyper(1),hhyper2(1),hdisc(1)],legend_str);
            set(hleg,'Fontsize',14,'FontName','Arial','Location','Best');

            %set(gca,'Ylim',ylimits);
            hold on

        end
    end

    if satplot==1

        %% Plot satellite matchup

        %Load data

        disp('Loading satellite data...');
        load('C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\all_sat_GC2_v2.mat');
        disp('Loaded.');

        [n,satfile_num]=size(all_sat_GC2);

        target_indx=[];
        for isat=1:satfile_num
            if ~isempty(strfind(all_sat_GC2(isat).fname,sat_test_str(im))) && ~isempty(strfind(all_sat_GC2(isat).stn,hyper_sta(im)))
                target_indx=isat;
            end
        end

        if ~isempty(target_indx)
            sat_rrs=all_sat_GC2(target_indx).sat_rrs;
            sat_time=all_sat_GC2(target_indx).time;
            sat_datenum=datenum(sat_time);
            %Search for 6 h window on either side of observation window
            good_sat_index=find(sat_datenum>datenum(time_str1)-0.25 & sat_datenum<datenum(time_str2)+0.25);

            sat_rrs=sat_rrs(good_sat_index,:);

            if size(sat_rrs,1)>1
                mean_sat_rrs=mean(sat_rrs,1,'omitnan');
                std_sat_rrs=std(sat_rrs,1,'omitnan');
            else
                mean_sat_rrs=sat_rrs;
                std_sat_rrs=zeros(1,length(sat_rrs));
            end                    

            sat_lambda=all_sat_GC2(target_indx).sat_lambda;

            %hsat=plot(sat_lambda,mean_sat_rrs,'k+','Linewidth',[2]);
            %hold on

            satcolor='#0F2080';
            hsat_error=errorbar(sat_lambda,mean_sat_rrs,std_sat_rrs,'+','Linewidth',1.5,'Markersize',10,'Color',satcolor);

            %hsat_std2_2=plot(sat_lambda,mean_sat_rrs-std_sat_rrs,'k:','Linewidth',[1.5]);
            hold on

        end

        if isempty(target_indx) | isnan(mean_sat_rrs)
            %Change satplot flag to 0 if no data
            satplot=0;
        end

        %% Create legend
        
        if hproplot==1 && satplot==1
            legend_str={'HyperSAS','HyperPro','Discrete','Satellite'};

            hleg=legend([hhyper(1),hhyper2(1),hdisc(1),hsat_error],legend_str);
            set(hleg,'Fontsize',14,'FontName','Arial','Location','Best');

        elseif hproplot==0 && satplot==1            
            legend_str={'HyperSAS','Discrete','Satellite'};

            hleg=legend([hhyper(1),hdisc(1),hsat_error],legend_str);
            set(hleg,'Fontsize',14,'FontName','Arial','Location','Best');

        elseif hproplot==1 && satplot==0
            legend_str={'HyperSAS','HyperPro','Discrete'};

            hleg=legend([hhyper(1),hhyper2(1),hdisc(1)],legend_str);
            set(hleg,'Fontsize',14,'FontName','Arial','Location','Best');

        elseif hproplot==0 && satplot==0
            legend_str={'HyperSAS','Discrete'};

            hleg=legend([hhyper(1),hdisc(1)],legend_str);
            set(hleg,'Fontsize',14,'FontName','Arial','Location','Best');
        end

    end

    %% Final formatting and printing 
    
    %Set axis limits to HyperSAS rrs initial settings
    set(gca,'Ylim',sas_ylimits);

    htitle=title(title_txt,'Fontsize',20,'FontName','Arial','Fontweight','normal');

    print_file=[inpath,'hypersas_hyperpro_rrs_',num2str(file_no),'_',num2str(run_no),'_',station_txt,'_',sat_test_string,'.tif'];

    print(hfig,print_file,'-dtiff','-r300');

    %print -r300 -dmeta;  %To print to Windows clipboard
end

disp('Execution complete');

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu        
  
% hyperpro_kd_averages.m
% 
% Program to load hyperpro data and calculate near surface kd and ed_0
%
% Syntax:  hyperpro_kd_averages.m
%
% Inputs:
%    1) Processed and compiled HyperPro data ('*.mat') produced by
%     'hyperpro_matfile_compile_GC#.m' and 'hyperpro_process_GC#.m'
%
% Outputs:
%    1) '*.mat' files with hyperpro kd average information
%    2) plot of HyperPro profile data
%   
% Other m-files required: None 
%
% MAT-files required: 
%  1) Processed and compiled HyperPro data ('*.mat') produced by
%     'hyperpro_matfile_compile_GC#.m' and 'hyperpro_process_GC#.m'
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 22 Sep 2025

%% ------------- BEGIN CODE --------------%

clearvars
clc
close all

% Load file locations and input parameters

for icr = 1:5
    crz = ['GC',num2str(icr)];
    
    switch crz
        case 'GC2'
            inpath = '\HyperPro_data_GC2\';
            cruise_txt = 'GC2';
    
            pigsheetpath = '\Apr2009\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC2';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_apr2009.xlsx'],...
                'Sheet',wksheet);
        case 'GC3'
            inpath = '\HyperPro_data_GC3\';
            cruise_txt = crz;
    
            pigsheetpath = '\Jul2009\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC3';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_jul2009.xlsx'],...
                'Sheet',wksheet);
        case 'GC4'
            inpath ='\HyperPro_data_GC4\';
            cruise_txt = 'GC4';  % Specify cruise for file loading and naming
    
            pigsheetpath = '\Nov2009\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC4';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_nov2009.xlsx'],...
                'Sheet',wksheet);
        case 'GC5'
            inpath ='\HyperPro_data_GC5\';
            cruise_txt = 'GC5';
    
            pigsheetpath = '\Mar2010\Filterpad\Spreadsheet\'; 
            wksheet='PigmentData_GC5';
            % Read table with station cluster and phytoplankton size fraction information
            stn_pig_sf=readtable([pigsheetpath,'pigment_size_fractions_and_station_clusters_aph_slope_mar2010.xlsx'],...
                'Sheet',wksheet);
    end
    
    sta_txt = 'g2';
    
    % Flag to save output if desired
    save_output = false;  %true; %
    
    procolor1='#F5973A';
    procolor2='#0071bc';
    procolor3='#cd2026';
    procolor4='#4c2c92';
    procolor5='#2e8540';
    
    disp('Loading data...');
    load([inpath,cruise_txt,'_hyperpro_processed.mat']);
    
    disp('Loaded.')
    
    [~,nsta]=size(pro_group_struct);
    
    % Initialize variables
    all_stations={};
    
    %% Begin loop to select station(s) and plot 
    all_kd_e0=table('Size',[nsta,11],'VariableTypes',{'string','double','double','string','cell','cell','cell','cell','cell','cell','cell'},...
        'VariableNames',{'Station','kd0_490','ed0_490','Cluster','ed0','es','kd0','ed0_sd','es_sd','kd0_sd','ed_lambda'});
    
    cluster=cell(nsta,1);
    stations=cell(nsta,1);

    lambda_min = 390;  % 399.8
    lambda_max = 710;   % 700.2

    edlambda = pro_group_struct{1}.ed_lambda{1};
    
    lamnum = length(edlambda(edlambda>=lambda_min & edlambda<=lambda_max));  % Number of wavelengths between lamda_min and lambda_max
    
    % Preallocation arrays for mean kd
    all_kd_calc_mean=nan(nsta,lamnum); 
    all_ed0_calc_mean=nan(nsta,lamnum); 
    all_es_calc_mean=nan(nsta,lamnum); 
    
    for ista=1:nsta 
    
        % Skip group if empty
        if isempty(pro_group_struct{ista})
            continue
        end
        
        single_sta=pro_group_struct{ista}.pro_station;  %Select station to plot
        % For repeated stations, add 'b' to station_id unless already
        %    labeled
        % if contains(single_sta{1},all_stations)
        %     single_sta=strrep(single_sta{1},single_sta,[single_sta{1},'b']);
        % end
        all_stations=cat(2,all_stations,single_sta(1));
        
        %To do all stations, comment out the if continue statement; otherwise only
              % station given in single_sta will be processed
        if ~contains(single_sta,sta_txt(1,:),'IgnoreCase',true)
            continue
        else
             sta_txt=char(single_sta);
        end
    
        cluster_indx = find(contains(stn_pig_sf.Station,lower(single_sta{1})));
        if isempty(cluster_indx)
            cluster{ista} = 'Unknown';
        else
            cluster{ista} = stn_pig_sf.Cluster{cluster_indx};
        end
    
        %% Get HyperPro station data and plot if desired
    
        profile_data=pro_group_struct{ista};
        [n,pro_num]=size(profile_data.pro_date);
    
        pres_interp=profile_data.pres_interp;
        k_ed=profile_data.k_ed;
        edz=profile_data.ed_corr_es;
        es=profile_data.es_corr_interp;
        pro_tilt=profile_data.pro_tilt;
    
        ed_lambda=profile_data.ed_lambda;

        kd_calc=nan(pro_num,lamnum); % 90
        ed0_calc=nan(pro_num,lamnum); % 90

        % pro_indices=1:pro_num;

        % Select good profiles  % THIS IS DONE AT THE PROCESSING LEVEL 
        % if icr == 2 && any(strcmp(cellstr(sta_txt),'h2'))
        %         pro_indices = 2;
        % elseif icr == 2 && any(strcmp(cellstr(sta_txt),'g4'))
        %         pro_indices = 1:3;
        % elseif icr == 2 && any(strcmp(cellstr(sta_txt),'e0'))
        %         continue;
        % end

        %% Plot spectra 
        for ipro=1:pro_num  
            plt_lambda=ed_lambda{ipro}(ed_lambda{ipro}>=lambda_min & ed_lambda{ipro}<=lambda_max);
            plt_kd=k_ed{ipro}(ed_lambda{ipro}>=lambda_min & ed_lambda{ipro}<=lambda_max);
            plt_edz=edz{ipro}(:,ed_lambda{ipro}>=lambda_min & ed_lambda{ipro}<=lambda_max);
            plt_es=es{ipro}(:,ed_lambda{ipro}>=lambda_min & ed_lambda{ipro}<=lambda_max);
            plt_edz(plt_edz<0.1)=NaN;  % This is to eliminate noisy data
            plt_edz(pro_tilt{ipro}>5,:)=NaN;
    
            % Adjust depth range based on kd_532 for log extrapolation to surface values
            if plt_kd(plt_lambda>532 & plt_lambda<532.6)>=2
                depth_rng=1;
            elseif plt_kd(plt_lambda>424 & plt_lambda<426)>=1.0 && plt_kd(plt_lambda>424 & plt_lambda<426)<2.0
                depth_rng=2;
            elseif plt_kd(plt_lambda>424 & plt_lambda<426)>=0.5 && plt_kd(plt_lambda>424 & plt_lambda<426)<1.0
                depth_rng=5;
            elseif plt_kd(plt_lambda>424 & plt_lambda<426)>.25 && plt_kd(plt_lambda>424 & plt_lambda<426)<0.5
                depth_rng=8;
            else
                depth_rng=10;
            end

            for lmbn=1:length(plt_lambda)

                good_indx = pres_interp{ipro}<depth_rng & ~isnan(plt_edz(:,lmbn));
                ed_mdl=polyfit(pres_interp{ipro}(good_indx),log(plt_edz(good_indx,lmbn)),1);
                kd_calc(ipro,lmbn)=-ed_mdl(1);
                ed0_calc(ipro,lmbn)=exp(ed_mdl(2));
                
                % Alternative method
                % ed_mdl=fitlm(pres_interp{ipro}(good_indx),log(plt_edz(good_indx,lmbn)),'Intercept',true);
                % if ed_mdl.Rsquared.Ordinary<0.7
                %     kd_calc(ipro,lmbn)=NaN;
                % else            
                %     kd_calc(ipro,lmbn)=-ed_mdl.Coefficients.Estimate(2);
                %     ed0_calc(ipro,lmbn)=-ed_mdl.Coefficients.Estimate(1);
                % end

            end

            kd_calc(kd_calc<=0)=NaN;
            ed0_calc(ed0_calc<=0)=NaN;

            % Eliminate rows containing NaN if desired
            % Create mask
            contains_number = (isnan(kd_calc));  % This creates a logical
            %    matrix of the same size as rad_data_array, where true indicates the presence
            %    of number_to_find and false indicates its absence. Identify rows
            %    containing the number.
            rows_to_delete = any(contains_number, 2); % The any(..., 2) function checks each
            %    row (2 indicates dimension 2, i.e., all elements in a given row) contains_number. If any element
            %    in a row is true, the corresponding element in rows_to_delete will be true.
            kd_calc(rows_to_delete,:) = NaN;  % kd_calc(rows_to_delete,:) = [];

            contains_number = (isnan(ed0_calc));  
            rows_to_delete = any(contains_number, 2); 
            ed0_calc(rows_to_delete,:) = NaN; % ed0_calc(rows_to_delete,:) = [];

            figure(1);
            clf
            hkd1=plot(plt_lambda,kd_calc(ipro,:),'-','Linewidth',2,'Color',procolor1);
            hold on
            hkd2=plot(plt_lambda,plt_kd,':','Linewidth',2,'Color',procolor2);
            xlabel('Wavelength');
            ylabel('Kd(\lambda)');
            title(sta_txt(ipro,:),'FontSize',12);
    
            figure(2);
            clf

            % 3D plot
            % hedz=plot3(plt_lambda',repmat(pres_interp{ipro},1,length(plt_lambda)),plt_edz);
            % set(gca,'YDir','reverse','ZScale','log','LineWidth',1.5,'FontSize',12);
            % xlabel('Wavelength (nm)','FontWeight','bold');
            % ylabel('Depth (m)','FontWeight','bold');
            % zlabel('{\itE_d(z,\lambda)} (\muW cm^{-2} nm^{-1})','FontWeight','bold');
            % box on

            % 2D plot
            hedz=plot(plt_lambda',plt_edz);
            set(gca,'YScale','log','LineWidth',1.5,'FontSize',12);
            xlabel('Wavelength (nm)','FontWeight','bold');
            ylabel('{\itE_d(z,\lambda)} (\muW cm^{-2} nm^{-1})','FontWeight','bold');
            box on
        end
    
        % Mean value for kd
        kd_calc_mean=mean(kd_calc,1,'omitnan');
        kd_calc_sd=std(kd_calc,1,'omitnan');
        all_kd_calc_mean(ista,:)=kd_calc_mean;
        kd_calc_mean(kd_calc_mean<0) = NaN; 
        
        %Mean values for ed0
        ed0_calc_mean=mean(ed0_calc,1,'omitnan');    
        ed0_calc_sd=std(ed0_calc,1,'omitnan');    
        all_ed0_calc_mean(ista,:)=ed0_calc_mean;
           
        %Mean values for es
        es_calc_mean = mean(plt_es,1,'omitnan');
        es_calc_sd = std(plt_es,1,'omitnan');
        all_es_calc_mean(ista,:)=es_calc_mean;
    
        % If desired, save profile data for subsequent calculations
        if save_output
            save([inpath,sta_txt(1,:),'_hyperpro_profile.mat'],'profile_data','ed0_calc_mean','es_calc_mean',...
            'kd_calc_mean','plt_lambda','kd_calc_sd','ed0_calc_sd','es_calc_sd');
        end
    
        kd490_indx = find(plt_lambda>488 & plt_lambda<491);
        all_kd_e0.Station(ista) = sta_txt(1,:);
        all_kd_e0.kd0_490(ista) = kd_calc_mean(kd490_indx);
        all_kd_e0.ed0_490(ista) = ed0_calc_mean(kd490_indx);
        all_kd_e0.Cluster(ista) = cluster(ista);
        all_kd_e0.kd0(ista) = {kd_calc_mean};
        all_kd_e0.kd0_sd(ista) = {kd_calc_sd};
        all_kd_e0.ed0(ista) = {ed0_calc_mean};
        all_kd_e0.ed0_sd(ista) = {ed0_calc_sd};
        all_kd_e0.es(ista) = {es_calc_mean};
        all_kd_e0.es_sd(ista) = {es_calc_sd};
        all_kd_e0.ed_lambda(ista) = {plt_lambda'};
        
       
    end
    
    if save_output
        save([inpath,'all_hyperpro_kd_ed0.mat'],'all_kd_e0');
    end
    
    %% Plot average kd, ed0 and es for all stations
    
    % Plot of mean kd
    % 
    figure(3);
    clf
    hold on
    
    % Plot average ed0 and es for all stations
    for ikd = 1:height(all_kd_e0)
        if isempty(all_kd_e0.kd0{ikd})
            continue
        else
            hkd=plot(all_kd_e0.ed_lambda{ikd},all_kd_e0.kd0{ikd},'-b','Linewidth',2);
            if ikd==1
                hold on
            end
        end
    end
    
    xlabel('Wavelength (nm)','FontWeight','bold');
    zlabel('{\itKd(\lambda)} (m^{-1})','FontWeight','bold');
    box on
    
    hleg = legend(hkd,{'Kd'});
    
    % Plot of mean ed0 and es
    
    figure(4);
    clf
    hold on
    
    % Plot average ed0 and es for all stations

    for ied = 1:height(all_kd_e0)
        if isempty(all_kd_e0.ed0{ied})
            continue
        else
            if max(all_kd_e0.ed0{ied})>200 || max(all_kd_e0.ed0{ied})<10 %Filter bad records
                continue
            end
    
            hed0=plot(all_kd_e0.ed_lambda{ied},all_kd_e0.ed0{ied},'-b','Linewidth',2);
            if ied==1
                hold on
            end
            hes=plot(all_kd_e0.ed_lambda{ied},all_kd_e0.es{ied},':m','Linewidth',2);
        end
    end
    
    xlabel('Wavelength (nm)','FontWeight','bold');
    zlabel('{\itEd(z,\lambda)} or {\itEs(z,\lambda)} (\muW cm^{-2} nm^{-1})','FontWeight','bold');
    box on
    
    hleg = legend([hed0(1),hes(1)],{'Ed0','Es'});
    
    % Write compiled data to spreadsheet
    if save_output
        writetable(all_kd_e0(:,1:4),[inpath,cruise_txt,'_all_kd_eo.xlsx'],'Sheet',1,'Range','A1')
    end
    
    %% Compute and plot cluster means
    
    est_indx = find(strcmp(cluster,'estuary')); 
    in_indx = find(strcmp(cluster,'inner')); 
    mid_indx = find(strcmp(cluster,'mid')); 
    out_indx = find(strcmp(cluster,'outer')); 
    
    meankdest = mean(all_kd_calc_mean(est_indx,:),1,'omitnan');
    meankdin = mean(all_kd_calc_mean(in_indx,:),1,'omitnan');
    meankdmid = mean(all_kd_calc_mean(mid_indx,:),1,'omitnan');
    meankdout = mean(all_kd_calc_mean(out_indx,:),1,'omitnan');
    
    % plot mean k for the different water mass clusters
    figure(5);
    clf
    
    hkdest = plot(plt_lambda,meankdest./meankdest(plt_lambda==492.16),'-','Linewidth',2);
    hold on
    hkdin = plot(plt_lambda,meankdin./meankdin(plt_lambda==492.16),'-','Linewidth',2);
    hkdmid = plot(plt_lambda,meankdmid./meankdmid(plt_lambda==492.16),'-','Linewidth',2);
    hkdout = plot(plt_lambda,meankdout./meankdout(plt_lambda==492.16),'-','Linewidth',2);
    
    set(gca,'YScale','log','Fontsize',14,'FontWeight','bold','YLim',[0.02,10],'YTick',[0.01,0.03,0.1,0.3,1,3,6],...
        'XTick',400:100:800,'XLim',[lambda_min,lambda_max]);
    % set(gca,'YTickLabels',{}); %,'YTickLabels',{}); 
    xlabel('Wavelength (nm)','FontWeight','bold');
    ylabel('{\itK_d(\lambda)} (m^{-1})','FontWeight','bold');
    
    % text(450,2,'Coastal','Fontsize',14);
    % text(450,0.4,'Mid-Shelf','Fontsize',14);
    % text(450,0.075,'Slope','Fontsize',14);
    box on
    
    % hleg = legend([hkdest,hkdin,hkdmid,hkdout],{'Estuary','Inner','Mid','Outer'},'Location','SouthEast');
    
    pbaspect([4 3 1]);
    
    if save_output
        set(gca,'YLim',[0.3,20]);
        print([inpath,'kd_cluster_averages_normalized_',crz,'.tif'],'-dtiff','-r600');
        %save([inpath,crz,'_kd_cluster_averages.mat'],'plt_lambda','meankdest','meankdin','meankdmid','meankdout');
    end

end

disp('Execution complete');

%------------- END OF CODE --------------
%Please send suggestions for improvement
% to Steven Lohrenz at this email address: slohrenz@umassd.edu        
  
% pp_daily_algorithm.m 
% 
% Absorption-based daily PP calculation - Program to estimate primary 
% production profile from from HyperPro-derived ed and kd values and 
% filterpad aph values and extrapolation of light level over photoperiod
%
% Syntax:  pp_daily_algorithm
%
% Inputs:
%   1) Folder location with photosynthetic parameter spreadsheet location
%       and PUR and aph data
%   2) Photosynthetic parameter data for a given station imported from
%       'PECurves_GC-phyto_GC#.xlsx' spreadsheet
%   3) 'Stn#_PUR_aph_dat.mat' files for a given station derived from 
%       generated using 'ed_calculation_daily_extrapolation.m'
%   4) 'Stn#_aph_dat*m.mat' files for filterpad aph slope generated using
%       'aph_slope_aph_mean_calculations.m'
%
% Outputs:
%    1) 'Stn#_PP_dat.mat' files with PAR and PUR profiles for primary production
%       calculations with units converted to mol Q or mol photons m-2 h-1
%    2) Plot of primary production profiles using station photosynthetic 
%       parameter data and model estimated parameters
%   
% Other m-files required:
% ed_calculation_daily_extrapolation.m
%
% MAT-files required: 
%    1) 'Stn#_PUR_aph_dat.mat' files for a given station derived from 
%       generated using 'ed_calculation_with_unit_conversion_allstation.m'
%    2) 'Stn#_aph_dat*m.mat' files for filterpad aph slope generated using
%       'aph_slope_aph_mean_calculations.m'
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 20 Oct 2024

%% ------------- BEGIN CODE --------------%

clf
clc
clearvars

folderpath='\';

% Create blank table for results
all_IPP_tab = table();

for icr = 2:5
    cruisename = ['GC',num2str(icr)];
    
    switch cruisename
        case 'GC2'
            aphdatpath='Apr2009\Filterpad\';
            % aphdatpath='Filterpad\GC2\';
            pe_dat=readtable([folderpath,'P-E\PECurves_GC-phyto_GC2.xlsx']);
        case 'GC3'
            aphdatpath='Jul2009\Filterpad\';
            % aphdatpath='Filterpad\GC3\';
            pe_dat=readtable([folderpath,'P-E\PECurves_GC-phyto_GC3.xlsx']);
        case 'GC4'
            aphdatpath='Nov2009\Filterpad\';
            % aphdatpath='Filterpad\GC4\';
            pe_dat=readtable([folderpath,'P-E\PECurves_GC-phyto_GC4.xlsx']);
        case 'GC5'
            aphdatpath='Mar2010\Filterpad\';
            % aphdatpath='Filterpad\GC5\';
            pe_dat=readtable([folderpath,'P-E\PECurves_GC-phyto_GC5.xlsx']);
    end
    
    sta_txt=pe_dat.Station;
    % sta_txt = {'b2'};
    
    % Flag for saving and printing output
    save_output = true;
    
    % Create a depth variable (NOT USED - DEPTH OF INTEGRATION BASED ON
    % EUPHOTIC DEPTH)
    % 
    % z = 0:1:100;
    % z = 0:2:38;
    
    P_int = zeros(size(sta_txt,1),1);
    P_est_int = zeros(size(sta_txt,1),1);
    P_std_int = zeros(size(sta_txt,1),1);
    P_constk_int = zeros(size(sta_txt,1),1);
    cluster = cell(size(sta_txt,1),1);

    for istn=1:length(sta_txt)
    
        sta_indx=find(strcmpi(pe_dat.Station,sta_txt(istn)));

        dt = pe_dat.Date(sta_indx);  % Date
        utc = pe_dat.UTC(sta_indx);  % UTC time
        lat = pe_dat.LatN(sta_indx);
        lon = pe_dat.LonW(sta_indx);

        chl=pe_dat.Tchla(sta_indx); % mg m-3
        pbmax=pe_dat.PBmax(sta_indx); % mgC mgChl-1 h-1
        phi=pe_dat.phi(sta_indx);  % mol C (mol quanta)-1
        alphab=pe_dat.alpha_mgC_mgChl_permolQ_perm2(sta_indx); %mgC mgChl-1 (mol quanta m-2)-1
        aph_440 = pe_dat.aph_440(sta_indx);
        % kd490 = pe_dat.Kd490(sta_indx);  % kd490 is retreived from
        %   HyperPro profile data
        cluster(istn)=pe_dat.Cluster(sta_indx);

        % Load light profile data for station
        % load([folderpath,aphdatpath,sta_txt{istn},'_PUR_aph_dat.mat']);
    
        % Load aph data for given station
        aph_file_list = dir([folderpath,aphdatpath,sta_txt{istn},'_aph_dat*m.mat']);
        if length(aph_file_list)>1
            for iaph=1:length(aph_file_list)
                strn1=strfind(aph_file_list(iaph).name,'dat_')+4;
                strn2=strfind(aph_file_list(iaph).name,'m.mat')-1;
                pad_depth=str2num(aph_file_list(iaph).name(strn1:strn2));
                if pad_depth<7
                    load([folderpath,aphdatpath,aph_file_list(1).name]);
                else
                    continue
                end
            end
        elseif isempty(aph_file_list)
            continue
        else        
            load([folderpath,aphdatpath,aph_file_list.name]);
        end

        % Full spectral daily PP

        % Get solar elevation and azimuth for photoperiod
        % mDatetm = dt + days((0:1440)'./1440);  % Increment days in minute intervals (60*24)
        mDatetm = dt + days((0:23)'./24);  % Increment days in hour intervals
        [sAz,sEl] = SolarAzEl_revised(mDatetm,zeros(size(mDatetm,1),1)+lat,zeros(size(mDatetm,1),1)-lon,zeros(size(mDatetm,1),1));

        % Find indices and times of sunrise and sunset
        sol_zen1 = 90 - sEl;
        cos_solzen_air = cosd(sol_zen1);
        utc_h_from_midnght = mod((0:23.98) - 5,24);
        cos_indx = find(cos_solzen_air>0); % Determine daytime indices

        % Apply Snell's law to get transmitted angle across air-water
        %  interface
        sol_zen2 = asind(sind(sol_zen1)./1.34);   % Divide sin of theta by water refractive index (1.34)
        cos_solzen_wtr = cosd(sol_zen2);

        % Optional to plot cos_solzen_wtr vs utc_h_from_midnght
        % clf
        % plot(utc_h_from_midnght(cos_solzen_wtr>0),cos_solzen_wtr(cos_solzen_wtr>0),'bo')

        % Adjust irradiance values by cosine of solar zenith angle (cos_solzen_wtr), referencing time surface
        %   irradiance measurement was acquired

        % Find time closest to irradiance profile
        [~,~,~,h,mi,s] = datevec(utc);
        pro_hr = h + mi/60 + s/3600;
        pro_tm_indx = floor(pro_hr);
        if isnan(pro_tm_indx)
            continue;
        end

        % Get cos_solzen for irradiance profile time
        cos_solzen_pro_tm = cos_solzen_wtr(pro_tm_indx);
        cos_solzen = cos_solzen_wtr./cos_solzen_pro_tm; % Normalize cosine 
        %  to cosine at time of profile and then extrapolate that over 
        %  photoperiod

        %% Calculate light profiles

        % Approximate value for Ed/Eo
        avg_cos = 0.8;
        
        % This calls a script to generate PUR_z_micro, PUR_z_pico, PAR_z
        ed_calculation_daily_extrapolation;

        % Skip station if no hyperpro data
        if isnan(kd490)
            continue
        end

        if log(abs(aph_slope_var))>-8.3
            PUR_z=PUR_z_micro; % mol Q or mol photons m-2 h-1
            aph440_mean = mean(aph440_micro);
        else 
            PUR_z=PUR_z_piconano; % mol Q or mol photons m-2 h-1
            aph440_mean = mean(aph440_piconano);
        end
        
        % Calculate pmax_aph440_spec
        pmax_aph440_spec = (pbmax.*chl./aph_440)./12000;   % mol C m-2 h-1

        % Estimated pmax_aph440_spec from pmax_box_plot.m
        pmax_aph440_spec_est = 0.0215*(1 - exp(-5.83*kd490));   % mol C m-2 h-1 
        
        % Estimated phi based on correlation with pmax_aph440_spec (see pmax_aph_spec_vs_phi.m)
        phi_est = -60.9.*pmax_aph440_spec_est.^2 + 5.19.*pmax_aph440_spec_est - 0.00369;
    
        % Use to determine +-50% changes from original WRM values

        % pmax_aph440_spec_fiftyplus=pmax_aph440_spec+(pmax_aph440_spec*.5);
        % pmax_aph440_spec_fiftyminus=pmax_aph440_spec-(pmax_aph440_spec*.5);
        % aph_440_fiftyplus=aph_440+(aph_440*.5);
        % aph_440_fiftyminus=aph_440-(aph_440*.5);
        % phi_fiftyplus=phi+(phi*.5); phi(phi>0.125) = 0.125;
        % phi_fiftyminus=phi-(phi*.5);
        % kd490_fiftyplus=kd490+(kd490*.5);
        % kd490_fiftyminus=kd490-(kd490*.5);
    
        % Run pp algorithm on profile for each hour, adjusting irradiance
        %   for cosine of solar zenith angle

        clear P_z_day P_z_std_day P_z_est_day P_z_PAR_day

        for icos = 1:length(cos_indx)
            % cos_zen = cos_solzen(cos_indx(icos))./cos_solzen_pro_tm;  % This adjusts the cosine to account for the time of day the profile was done

            P_z = aph_440.*pmax_aph440_spec.*(1 - exp(-phi.*(PUR_z(:,icos))./(pmax_aph440_spec))); % mol C m-3 h-1
        
            % for sensitivity analysis modify the required variables each time
            % P_z = aph_440_fiftyminus.*pmax_aph440_spec.*(1 - exp(-aph_440_fiftyminus.*phi.*(PUR_z./avg_cos)./(aph_440_fiftyminus.*pmax_aph440_spec))); % mol C m-3 h-1
            
            % PAR-based PP using conventional P-E parameters
            P_z_std = chl.*pbmax.*(1 - exp(-alphab.*(PAR_z(:,icos))./pbmax))./12000; % mol C m-3 h-1
            
            % PAR-based PP using conventional P-E parameters and constant KPAR
            P_z_PAR = chl.*pbmax.*(1 - exp(-alphab.*(PAR_z_constK(:,icos))./pbmax))./12000; % mol C m-3 h-1
            
            % Full spectral PP estimated using predicted pmax_aph440_spec
            P_z_est = aph_440.*pmax_aph440_spec_est.*(1 - exp(-phi_est.*(PUR_z(:,icos))./(pmax_aph440_spec_est))); % mol C m-3 h-1
        
            if exist('P_z_day','var')
                P_z_day = P_z_day + P_z;
                P_z_std_day = P_z_std_day + P_z_std;
                P_z_PAR_day = P_z_PAR_day + P_z_PAR;
                P_z_est_day = P_z_est_day + P_z_est;
            else
                P_z_day = P_z;
                P_z_std_day =  P_z_std;
                P_z_PAR_day = P_z_PAR;
                P_z_est_day = P_z_est;
            end
            
        end

        % Calculate water-column integrated primary production (mol C m-2 d-1)
        
        P_int(istn) = sum(P_z_day.*z_incr);
        P_est_int(istn) = sum(P_z_est_day.*z_incr);
        P_std_int(istn) = sum(P_z_std_day.*z_incr);
        P_constk_int(istn) = sum(P_z_PAR_day.*z_incr);  % Constant KPAR assumed
    
        %% Plot figure
        
        figure(1);
        clf

        hmodel = plot(P_z_day,z,'-b','Linewidth',5);%2.5);
        hold on

        set(gca,'YDir','reverse','XAxisLocation','top','XScale','log','FontSize',16);
        xlim auto

        % hest = plot(P_z_est,z,'-+c','Linewidth',1.5);
        hest = plot(P_z_est_day,z,'--','Color',[0.5 0.5 0.5],'Linewidth',4);

        hstd = plot(P_z_std_day,z,'-m','Linewidth',5);
        % hconstK = plot(P_z_PAR,z,'-*','Color','magenta','Linewidth',1.5);
        % hconstK = plot(P_z_PAR,z,'--','Color',[0.7 0.7 0.7],'Linewidth',1.5);

        ylabel('Depth (m)','FontWeight','bold');
        xlabel('PP(z) (mol C m^{-3} d^{-1})','FontWeight','bold');
        title(sta_txt{istn});
        title_text = sprintf('%s - %s\n%s', cruisename, sta_txt{istn}, char(cluster(istn)));
        t = title(title_text, 'FontSize', 16, 'FontWeight', 'normal');

        % Get the X and Y limits of the plot
        xLimits = xlim;
        yLimits = ylim;
        t.Position = [xLimits(1) * 3, yLimits(1) + (yLimits(2) - yLimits(1)) * 0.17, 0];

        grid on
        ax=gca;
        ax.XMinorTick="off";
        ax.YMinorTick="off";
        % ax.XMinorGrid="off";
        % ax.YMinorGrid="off";
        ax.GridColor='k';
        % ax.MinorGridColor='b';
        % ax.GridLineStyle="--";
        ax.MinorGridLineStyle="none";
        ax.LineWidth=1;

        ax.FontSmoothing="on";
        ax.XColor='k';
        ax.YColor='k';
        % ax.XTickLabel = {}; % Uncomment this if don't want X-axis tick_label
        % set(gcf,'Color','W','InvertHardcopy', 'off'); % Turn off image
        %   background color

        hlgnd = legend([hmodel,hest,hstd],{'Wavelength-Resolved Model','Wavelength-Resolved Model Estimate', ...
            'Wavelength-Integrated'},...
            'FontSize',14,'Location','Southeast');


%% Save output for individual stations
        PP_dat = struct('z',{z'},'PUR_z',{PUR_z},'P_z',{P_z},'P_z_std',{P_z_std},'P_z_est',{P_z_est},'P_z_PAR',{P_z_PAR}); %,'PAR_z',{PAR_z'},'PAR_z_constK',{PAR_z_constK'},'P_z_PAR',{P_z_PAR});

        if save_output
            save([folderpath,'P-E\',cruisename,'_',sta_txt{istn},'_PP_dat.mat'],'PP_dat');

            disp(['Printing ','pp_',sta_txt{istn},'.tif']);
            print([folderpath,'P-E\',cruisename,'_pp_',sta_txt{istn},'.tif'],'-dtiff','-r600');

        end

    end

    IPP_tab = table(repmat({cruisename},size(sta_txt(P_int~=0 & ~isnan(P_int)))),...
        sta_txt(P_int~=0 & ~isnan(P_int)),cluster(P_int~=0 & ~isnan(P_int)),P_int(P_int~=0 & ~isnan(P_int)),...
        P_est_int(P_int~=0 & ~isnan(P_int)),P_std_int(P_int~=0 & ~isnan(P_int)),P_constk_int(P_int~=0 & ~isnan(P_int)));
    IPP_tab.Properties.VariableNames = ["Cruise","Station","Cluster","WRM","WRME","WIM","WIMConstk"];
    IPP_tab = sortrows(IPP_tab,1);

    all_IPP_tab = [all_IPP_tab;IPP_tab];

end

all_IPP_tab = sortrows(all_IPP_tab,["Cruise","Station"]);

%% Save output for compiled IPP results
 
if save_output
    save([folderpath,'P-E\all_IPP_table.mat'],'all_IPP_tab');
    
    filename = 'all_IPP_tab.xlsx';
    writetable(all_IPP_tab,[folderpath,'P-E\',filename],'Sheet',1,'Range','A1','WriteMode','overwritesheet')
end

%% Plot WRM vs WIM results

figure(2);
clf

hpp = plot(all_IPP_tab.WRM,all_IPP_tab.WIM,'ob','Linewidth',1.5);
hold on
h1to1 = plot([0.03,2],[0.03,2],':k','LineWidth',1.5);

[R_IPP,P] = corrcoef(all_IPP_tab.WRM,all_IPP_tab.WIM);
rtmeansq = rmse(all_IPP_tab.WRM,all_IPP_tab.WIM);
rsquared = R_IPP(1,2).^2;

hpe = plot(all_IPP_tab.WRM,all_IPP_tab.WRME,'+m','Linewidth',1.5);

ylabel('WIM or WRME IPP (mol C m-2 h^{-1})','FontWeight','bold');
xlabel('WRM IPP (mol C m-2 h^{-1})','FontWeight','bold');

set(gca,'XScale','log','YScale','log','XLim',[0.03,2],'YLim',[0.03,2]);

legend([hpp,hpe,h1to1],{'WRM vs WIM','WRM vs WRME','1:1 Line'},'Location','northwest');

disp('Completed');



% pmax_bot_plot.m
% Syntax: pmax_box_plot
%
% This script retrieves the Pmax and Kd490 data from a spreadsheet and
% generates box plots for both and then computes a regression fit for the
% averages for the given clusters
%
% Inputs:
%    1) Folder location with 'Pmax_aph_spec_vs_phi.xlsx' spreadsheet
%    4) Spreadsheet with stations and depths for sample collection
%    
% Outputs:
%    output - Figures with plotted results
%   
% Other m-files required: None
% 
% MAT-files required: None
%
% Author: Steven E. Lohrenz
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 15 Apr 2024

%% ------------- BEGIN CODE --------------%% 

clc
clearvars
close all

outpath ='\P-E\';

% Read in data to table
for icr = 1:5
    crz = ['GC',num2str(icr)];
    
    switch crz
        case 'GC1'
            inpath = '\GC1_HyperPro_data_GC1\';
            cruise_txt = 'GC1';
            filename = [cruise_txt,'_all_kd_eo.xlsx'];
            kdata1 =  readtable([inpath,filename],'Sheet','Sheet1');
        case 'GC2'
            inpath = '\\HyperPro_data_GC2\';
            cruise_txt = 'GC2';
            filename = [cruise_txt,'_all_kd_eo.xlsx'];
            kdata2 =  readtable([inpath,filename],'Sheet','Sheet1');
        case 'GC3'
            inpath = '\HyperPro_data_GC3\';
            cruise_txt = crz;
            filename = [cruise_txt,'_all_kd_eo.xlsx'];    
            kdata3 =  readtable([inpath,filename],'Sheet','Sheet1');
        case 'GC4'
            inpath ='\HyperPro_data_GC4\';
            cruise_txt = 'GC4';  % Specify cruise for file loading and naming
            filename = [cruise_txt,'_all_kd_eo.xlsx'];    
            kdata4 =  readtable([inpath,filename],'Sheet','Sheet1');
        case 'GC5'
            inpath ='\HyperPro_data_GC5\';
            cruise_txt = 'GC5';
            filename = [cruise_txt,'_all_kd_eo.xlsx'];    
            kdata5 =  readtable([inpath,filename],'Sheet','Sheet1');
    end
end

% Prepare datasets
% Kd data
kdata = [kdata1;kdata2;kdata3;kdata4;kdata5];  % Combine Kd data for cruises
kdata = kdata(~(kdata.kd0_490==0 | strcmp(kdata.Cluster,'Unknown')),:);  % Eliminate missing entries
kgroupdata = categorical(kdata.Cluster); % Set water mass Cluster as categorical variable

% Photosynthesis-irradiance (P-E) data 
pefile = 'Pmax_aph_spec_vs_phi_table.xlsx';
pepath = '\P-E\';
pedata =  readtable([pepath,pefile],'Sheet','Data');
pegroupdata = categorical(pedata.Cluster); % Set water mass Cluster as categorical variable


%% Plot box chart of Pmax_aph440_spec for different clusters
figure(1);
clf
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size

boxchart(pegroupdata,pedata.Pmax_aph440_spec);
hold on
box on

xlabel('Cluster','FontSize',14,'FontWeight','bold');
ylabel('{\itP_{max}^{aph440}} (mol C m^{-2} h^{-1})','FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12);

print([outpath,'cluster_pmax_box_plot.tif'],'-dtiff','-r600');

% ANOVA for P-E data

% Anova comparison of group means
% [p_pmax,t_pmax,stats_pmax] = anova1(pedata.Pmax_aph440_spec,pedata.Cluster);
% [c_pmax,m_pmax,h_pmax,gnames_pmax] = multcompare(stats_pmax);


%% Plot box chart of Kd490 for different clusters
figure(2);
clf
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size

boxchart(kgroupdata,kdata.kd0_490);
hold on
box on
xlabel('Cluster','FontSize',14,'FontWeight','bold');
ylabel('\it{K_d(490)} (m^{-1})','FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12);

print([outpath,'cluster_kd490_box_plot.tif'],'-dtiff','-r600');

% ANOVA for Kd490

% Anova comparison of group means
% [p_kd,t_kd,stats_kd] = anova1(kdata.kd0_490,kdata.Cluster);
% [c,m,h,gnames] = multcompare(stats_kd);

%% Apply fit of cluster means of Pmax_aph440_spec vs Kd490

Pmax_group_mean = groupsummary(pedata(:,["Cluster","Pmax_aph440_spec","Kd490"]),"Cluster",{"mean","std"});
[r,p] = corrcoef(Pmax_group_mean.mean_Kd490,Pmax_group_mean.mean_Pmax_aph440_spec);

% Linear least squares fit to determine relationship

kdata_mean = Pmax_group_mean.mean_Kd490;
ydata_mean = Pmax_group_mean.mean_Pmax_aph440_spec;

yneg = Pmax_group_mean.std_Pmax_aph440_spec;
ypos = Pmax_group_mean.std_Pmax_aph440_spec;
xneg = Pmax_group_mean.std_Kd490;
xpos = Pmax_group_mean.std_Kd490;

% Plot relationship between cluster means of Pmax_aph440_spec vs Kd490
figure(3);
clf
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size

hpmax = errorbar(kdata_mean,ydata_mean,yneg,ypos,xneg,xpos,'ko');
hold on
box on
set(gca,'FontSize',12,'XLim',[-0.08,2]);

% Function for curve fit
func = @(x,kdata)x(1)*(1 - exp(-x(2)*kdata));

x0 = [0.25,0.2]; % Initial estimates
[x,resnorm,residual,exitflag,output] = lsqcurvefit(func,x0,kdata_mean,ydata_mean);

% Determine fit statistics
yest = x(1)*(1 - exp(-x(2).*kdata_mean));
[R,P] = corrcoef(ydata_mean, yest);
rtmeansq = rmse(ydata_mean,yest);

% Plot fitted curve
pltx = 0.005:0.01:1.75;
yest_plt = x(1)*(1 - exp(-x(2).*pltx));
hest = plot(pltx,yest_plt,':m','LineWidth',1.75);
xlabel('{\itK_d(490)} (m^{-1})','FontSize',14,'FontWeight','bold');
ylabel('{\itP_{max}^{aph440}} (mol C m^{-2} h^{-1})','FontSize',14,'FontWeight','bold');

% Annotate plot
fit_txt1 = ['{\itP_{max}^{aph440}} = ',num2str(x(1),'%5.4f'),'(1 - exp[-',num2str(x(2),'%4.2f'),'\it{K_d(490)}])'];
ht1 = text(0.7,0.005,fit_txt1,"FontSize",10);
ht2 = text(0.7,0.00275,['N = ',num2str(4),'; \it{r^2} = ',num2str(R(2,1).^2,'%5.3f'),...
    '; RMSE = ',num2str(rtmeansq,'%6.5f')],'FontSize',10);

print([outpath,'pmax_aph440_vs_kd490_curve_fit.tif'],'-dtiff','-r600');

disp('Completed');


%End of code


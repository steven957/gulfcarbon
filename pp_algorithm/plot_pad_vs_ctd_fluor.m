% plot_pad_vs_ctd_fluor.m
% Script plots relationship between aph and total chlorophyll a with ctd fluorescence
%
% Syntax: plot_pad_vs_ctd_fluor
%
% Inputs:
%    1) Folder location with compiled filterpad data in 'all_pad.mat' file
%      generated with 'pad_extract_seabass_v2.m'
%    2) Folder location with compiled pigment data in 'all_hplc_gc#.mat' generated using
%      'hplc_seabass_output_gulfcarbon.m'
%    3) Folder location with compiled SBE911 .cnv data files ('ctd_alldat.mat') produced 
%      using 'readctd_nonseabass.m'
%
% Outputs:
%    1) Plot of filterpad absorption vs ctd fluorescence
%    2) Plot of TChla vs ctd fluorescence
%   
%
% MAT-files required: 
%    1) 'all_pat.mat' (Compiled filterpad data from 'pad_extract_seabass_v2.m')
%    2) 'all_hplc_gc#.mat' generated using hplc_seabass_output_gulfcarbon.m'
%    3) 'ctd_bottle_all.mat' produced using 'readctd_nonseabass_rosette_files.m'
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 15 Oct 2025

%% ------------- BEGIN CODE --------------%% 

clearvars
clc


% Specify single station if desired; otherwise comment out to do all
% sta_txt='F5_0m';

% Set to true to run all stations
run_flag = 1;

all_aph_ctd_fluor = table; % Create empty table to populate with data
cruise_number_selection = 1:5;
outpath = '\';

for icr = cruise_number_selection %1:5
    cruise_name = ['gulfcarbon',num2str(icr)];
    
    % Folder path and file information for each cruise
    switch cruise_name
        case 'gulfcarbon1'
            padpath='\Jan2009\Filterpad\To_SeaBASS\';
            ctdpath='\Jan2009\CTD\';
            pigpath='\Pigments\';
        case 'gulfcarbon2'
            padpath='\Apr2009\Filterpad\To_SeaBASS\';
            ctdpath='\Apr2009\CTD\';
            pigpath='\Pigments\';
        case 'gulfcarbon3'
            padpath='\Jul2009\Filterpad\To_SeaBASS\';
            ctdpath = '\Jul2009\CTD\';
            pigpath='\Pigments\';
        case 'gulfcarbon4'
            padpath='\Nov2009\Filterpad\To_SeaBASS\';
            ctdpath = '\Nov2009\CTD\';
            pigpath='\Pigments\';
        case 'gulfcarbon5'
            padpath='\Mar2010\Filterpad\To_SeaBASS\';
            ctdpath = '\Mar2010\CTD\';
            pigpath='\Pigments\';
    end
    
    % Load pad data
    padfile=[padpath,'all_pad.mat'];
    load(padfile);
    
    % Load the CTD data
    ctdfile = [ctdpath, 'ctd_bottle_all.mat'];
    load(ctdfile);
    
    % Find column indices for targeted variables
    depth_indx = find(contains(ctd_bottle_all(1).names,'Depth'));
    fluor_indx = find(contains(ctd_bottle_all(1).names,'Fluorescence'));

    % Load pigment data
    pigfile = [pigpath,'all_hplc_gc',num2str(icr),'.mat'];
    hplc_dat = load(pigfile);
    
    % Create tables for each input dataset
    ctd_data = table('Size',[length(all_pad.name),4],'VariableTypes',{'string','string','double','double'},'VariableNames',["Cruise","Station","Depth","Fluor"]);
    pad_data = table('Size',[length(all_pad.name),4],'VariableTypes',{'string','string','double','double'},'VariableNames',["Cruise","Station","Depth","aph440"]);
    pig_data = table('Size',[length(all_pad.name),4],'VariableTypes',{'string','string','double','double'},'VariableNames',["Cruise","Station","Depth","TChla"]);
    
    %% Loop through pad samples
    for ipad = 1:length(all_pad.name)
        
        if run_flag
            indx1 = ipad;
        else
            indx1=find(datetime(all_pad.date)>datetime('12/31/08') & ...
            contains(all_pad.name,sta_txt)); 
        end
        
        pad_hdr=all_pad.fields{1};
        pad_units=all_pad.units{1};
        %pad_hdr{2}='a_tot'; %Leave as "a_tot"
        
        var={'aph'};  %Select either "ap", "ad", or "aph"
        var_str=var;
        var_char_num=length(var);
        
        [~,nvar]=size(var);
        var_indx=zeros(nvar,1);
        for iv=1:nvar
            var_indx(iv)=find(strcmp(pad_hdr,var(iv))==1);
            lambda_indx(iv)=find(strcmp(pad_hdr,'wavelength'));
        end
        %pad_lambda=str2num(char(pad_lambda));
        
        pad_dat=cell(nvar,1);
        pad_depth=cell(nvar,1);
        pad_station=cell(nvar,1);
        pad_lambda=cell(nvar,1);
        max_dat=0;
        min_dat=100;
        
        % Loop through pad variables (if necessary)
        for ivar=1:nvar
            % Loop through reps
            for irep=1:length(indx1)
                jrep=indx1(irep);
                
                pad_depth(ivar,1)=all_pad.depth(jrep);
                pad_dat(ivar,1)={all_pad.data{jrep}(:,var_indx(ivar))};
                pad_lambda(ivar,1)={all_pad.data{jrep}(:,lambda_indx(ivar))};
                lambda_440_indx = find(pad_lambda{ivar,1}>=435 & pad_lambda{ivar,1}<445);
                pad_station(ivar,1)=all_pad.station(jrep);
                pad_440 = pad_dat{ivar,1}(lambda_440_indx);
                pad_data(ipad,:) = {cruise_name,all_pad.station{jrep},str2num(all_pad.depth{jrep}),mean(pad_440)};
 
                % Find max for plotting
                new_max_dat=max(all_pad.data{jrep}(:,var_indx(ivar)),[],1);
                new_min_dat=min(all_pad.data{jrep}(:,var_indx(ivar)),[],1);
                if new_max_dat > max_dat
                    max_dat=new_max_dat;
                end
                if new_min_dat < min_dat
                    min_dat=new_min_dat;
                end
            end
        end
        
    pad_data.aph440(pad_data.aph440<=0)=NaN; 
    
    %% Find matching CTD data
    % CTD fields (for GC1-3,5 only)
    %{
        {'scan: Scan Count'                                         }
        {'nbf: Bottles Fired'                                       }
        {'prDM: Pressure, Digiquartz [db]'                          }
        {'depSM: Depth [salt water, m]'                             }
        {'t090C: Temperature [ITS-90, deg C]'                       }
        {'c0S/m: Conductivity [S/m]'                                }
        {'sigma-t00: Density [sigma-t, Kg/m^3 ]'                    }
        {'t190C: Temperature, 2 [ITS-90, deg C]'                    }
        {'c1S/m: Conductivity, 2 [S/m]'                             }
        {'sigma-t11: Density, 2 [sigma-t, Kg/m^3 ]'                 }
        {'sal00: Salinity [PSU]'                                    }
        {'sal11: Salinity, 2 [PSU]'                                 }
        {'sbeox0Mg/L: Oxygen, SBE 43 [mg/l]'                        }
        {'bat: Beam Attenuation, Chelsea/Seatech/Wetlab CStar [1/m]'}
        {'flC: Fluorescence, Chelsea Aqua 3 Chl Con [ug/l]'         }
        {'par: PAR/Irradiance, Biospherical/Licor'                  }
        {'latitude: Latitude [deg]'                                 }
        {'longitude: Longitude [deg]'                               }
        {'timeS: Time, Elapsed [seconds]'                           }
        {'flag:  0.000e+00'                                         }
    %}
        for ictd = 1:length(ctd_bottle_all)
            % Extract matching CTD data for the current pad station
            if strcmpi(ctd_bottle_all(ictd).sta, pad_station{ivar})
                ctd_indx = find(abs(ctd_bottle_all(ictd).ctddat(:,depth_indx) - str2num(char(pad_depth))) < 5); % Find index of matching station
                ctd_data(ipad,:) = {cruise_name,ctd_bottle_all(ictd).sta,mean(ctd_bottle_all(ictd).ctddat(ctd_indx,depth_indx)),mean(ctd_bottle_all(ictd).ctddat(ctd_indx,fluor_indx))}; % Store matching CTD data
                continue
            end
        end
        % If station not found, continue to next
        if ctd_data.Fluor(ipad)==0
            continue
        end

    %% Find matching pigment data

        for ipig = 1:height(hplc_dat.all_hplc)

            % Extract matching CTD data for the current pad station
            if strcmpi(hplc_dat.all_hplc.Station(ipig), pad_station{ivar})
                pig_indx = find(abs(hplc_dat.all_hplc.SamplingDepth_meters_(ipig) - str2num(char(pad_depth))) < 5); % Find index of matching station
                pig_data(ipad,:) = {cruise_name,hplc_dat.all_hplc.Station(ipig),hplc_dat.all_hplc.SamplingDepth_meters_(ipig),hplc_dat.all_hplc.x_TChlA_(ipig)}; % Compile matching pig data
                continue
            end
        end
        % If station not found, continue to next
        if pig_data.TChla(ipad)==0
            continue
        end
       
    end
    
    % Clean data
    ctd_data.Fluor(ctd_data.Fluor<=0.01 | ctd_data.Fluor>=100) = NaN;
    pig_data.TChla(pig_data.TChla==0) = NaN;
    all_aph_ctd_fluor = [all_aph_ctd_fluor;table(ctd_data.Cruise,ctd_data.Station,ctd_data.Depth,pig_data.Depth,pad_data.Depth,pad_data.aph440,ctd_data.Fluor,pig_data.TChla)];


end

all_aph_ctd_fluor.Properties.VariableNames = [{'Cruise'},{'Station'},{'CTD Depth'},{'Pigment Depth'},{'Pad Depth'},{'aph440'},{'Fluor'},{'TChla'}];
all_aph_ctd_fluor = all_aph_ctd_fluor(~isnan(all_aph_ctd_fluor.aph440) & ~isnan(all_aph_ctd_fluor.Fluor),:);

%% Plot aph vs Fluor results
xlab='Filterpad Absorption (m^{-1})';

% Override selection of min and max if desired
% min_dat=0;
% max_dat=1;

f1 = figure(1);
clf
scrsz = get(groot,'ScreenSize');
set(f1,'Position',[scrsz(4).*.4 scrsz(3).*.03,scrsz(3).*.4 scrsz(4).*.85])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

hplot1=[];
for icr = cruise_number_selection
    for k=1:nvar
        marker_indx = strcmp(all_aph_ctd_fluor.Cruise,['gulfcarbon',num2str(icr)]);
        hplot1(icr)=plot(all_aph_ctd_fluor.Fluor(marker_indx),all_aph_ctd_fluor.aph440(marker_indx),'+','Linewidth',1.2);
    
        set(gca,'Linewidth',1.3,'fontname','times','fontweight','bold','fontsize',18,...
                'Yscale','log','XScale','log','XLim',[0.02,50],'YLim',[0.001,5]);
        % ylabel('aph440 (m^{-1})','fontname','times','fontweight','bold','fontsize',16);
        % xlabel('CTD Chl Fluorescence','fontname','times','fontweight','bold','fontsize',16);

        set(hplot1(icr),'Linewidth',1.5);
        hold on
    end
end

% Get correlation statistics
[R_Fluor,P1] = corrcoef(all_aph_ctd_fluor.Fluor,all_aph_ctd_fluor.aph440);
rtmeansq1 = 10.^rmse(log(all_aph_ctd_fluor.Fluor),log(all_aph_ctd_fluor.aph440));
r_square_Fluor = R_Fluor(1,2).^2;

% Fit relationship to power function and plot line
curvefit1 = fit(all_aph_ctd_fluor.Fluor,all_aph_ctd_fluor.aph440,'power1');
xdt = 0.08:0.5:50;
hft = plot(curvefit1);
set(hft,'Color','k','LineStyle','--');

set(hft,'DisplayName',"a_{ph}(440) vs Fluor");

hft_ax = get(gca);

hft_ax.YLabel.String = "a_{ph}(440) (m^{-1})";
hft_ax.XLabel.String = 'CTD Chl Fluorescence';

leg_labels = {'GC1','GC2','GC3','GC4','GC5'};
hleg = legend(hplot1(cruise_number_selection),leg_labels(cruise_number_selection),'Location','Northwest');

text(.09,3,['a_{ph}(440) = ' char(num2str(curvefit1.a,'%.4f')) '*Chl\_Fluor^{' num2str(curvefit1.b,'%.3f'),'}'],'fontname','times','fontweight','bold','fontsize',14);
text(.09,1.7,['r^{2} = ' num2str(r_square_Fluor,'%.3f')],'fontname','times','fontweight','bold','fontsize',14);

pbaspect([3 2 1]);
fig_pos=get(gca,'position');
fig_pos(2)=fig_pos(2)-0.0;
set(gca,'Position',fig_pos);

print([outpath,'ctd_fluor_vs_path_aph_',cell2mat(leg_labels(cruise_number_selection))],'-dtiff','-r600');

%% Plot TChla vs Fluor results
xlab='Total Chl a (mg m^{-3})';

% Override selection of min and max if desired
% min_dat=0;
% max_dat=1;

f2 = figure(2);
clf
scrsz = get(groot,'ScreenSize');
set(f2,'Position',[scrsz(4).*.4 scrsz(3).*.03,scrsz(3).*.4 scrsz(4).*.85])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

hplot2=[];
for icr = cruise_number_selection
    for k=1:nvar
        marker_indx = strcmp(all_aph_ctd_fluor.Cruise,['gulfcarbon',num2str(icr)]);
        hplot2(icr)=plot(all_aph_ctd_fluor.Fluor(marker_indx),all_aph_ctd_fluor.TChla(marker_indx),'o','Linewidth',1.2);
    
        set(gca,'Linewidth',1.3,'fontname','times','fontweight','bold','fontsize',18,...
                'Yscale','log','XScale','log','XLim',[0.006,70],'YLim',[0.002,100]);
        % ylabel('aph440 (m^{-1})','fontname','times','fontweight','bold','fontsize',16);
        % xlabel('CTD Chl Fluorescence','fontname','times','fontweight','bold','fontsize',16);

        set(hplot2(icr),'Linewidth',1.5);
        hold on
    end
end

% Get correlation statistics
[R_Fluor2,P2] = corrcoef(all_aph_ctd_fluor.Fluor(~isnan(all_aph_ctd_fluor.TChla)),all_aph_ctd_fluor.TChla(~isnan(all_aph_ctd_fluor.TChla)));
rtmeansq2 = 10.^rmse(log(all_aph_ctd_fluor.Fluor(~isnan(all_aph_ctd_fluor.TChla))),log(all_aph_ctd_fluor.TChla(~isnan(all_aph_ctd_fluor.TChla))));
r_square_Fluor2 = R_Fluor2(1,2).^2;

% Fit relationship to power function and plot line
curvefit2 = fit(all_aph_ctd_fluor.Fluor(~isnan(all_aph_ctd_fluor.TChla)),all_aph_ctd_fluor.TChla(~isnan(all_aph_ctd_fluor.TChla)),'power1');
xdt = 0.02:0.5:70;
hft = plot(curvefit2);
set(hft,'Color','k','LineStyle','--');

set(hft,'DisplayName',"TChla vs Fluor");

hft_ax = get(gca);

hft_ax.YLabel.String = "Total Chl-a (mg m^{-3})";
hft_ax.XLabel.String = 'CTD Chl Fluorescence';

leg_labels = {'GC1','GC2','GC3','GC4','GC5'};
hleg = legend(hplot2(cruise_number_selection),leg_labels(cruise_number_selection),'Location','Northwest');

text(.035,55,['TChla = ' char(num2str(curvefit2.a,'%.4f')) '*Chl\_Fluor^{' num2str(curvefit2.b,'%.3f'),'}'],'fontname','times','fontweight','bold','fontsize',14);
text(.035,32,['r^{2} = ' num2str(r_square_Fluor2,'%.3f')],'fontname','times','fontweight','bold','fontsize',14);

pbaspect([3 2 1]);
fig_pos=get(gca,'position');
fig_pos(2)=fig_pos(2)-0.0;
set(gca,'Position',fig_pos);

print([outpath,'ctd_fluor_vs_chl_',cell2mat(leg_labels(cruise_number_selection))],'-dtiff','-r600');

disp('Completed');
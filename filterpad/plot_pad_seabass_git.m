% plot_pad_seabass_git.m - Plot selected filterpad data compiled from SeaBASS files 
% Syntax: plot_pad_seabass_git
%
% Inputs:
%    1) Folder location with compiled filterpad data in 'all_pad.mat' file
%    generated with 'pad_extract_seabass_git.m'
%
% Outputs:
%    output - Plots of fitlerpad absorption 
%   
%
% MAT-files required: 'all_pad.mat' (Compiled filterpad data from 'pad_extract_seabass_git.m')
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 21 November 2025

%% ------------- BEGIN CODE --------------%% 

clearvars
clc

% Folder path and file information for each cruise

cruise_name = 'gulfcarbon1';

switch cruise_name
    case 'gulfcarbon1'
        inpath='YOUR FOLDER WITH "all_pad.mat" FILE';
        file_pref='GulfCarbon1_FilterPad_';
         
    case 'gulfcarbon2'
        inpath='YOUR FOLDER WITH "all_pad.mat" FILE';
        file_pref='GulfCarbon2_FilterPad_';

    case 'gulfcarbon3'
        inpath='YOUR FOLDER WITH "all_pad.mat" FILE';
        file_pref='GulfCarbon3_FilterPad_';
        
    case 'gulfcarbon4'
        inpath='YOUR FOLDER WITH "all_pad.mat" FILE';
        file_pref='GulfCarbon4_FilterPad_';
        
    case 'gulfcarbon5'
        inpath='YOUR FOLDER WITH "all_pad.mat" FILE';
        file_pref='GulfCarbon5_FilterPad_';
end


% Specify single station as below if desired; otherwise set run_flag to 1 to plot all files
sta_txt='F5_0m';
run_flag = 0;  

infile=[inpath,'all_pad.mat'];

load(infile);

% Loop through stations
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
    
    [m,n]=size(var);
    var_indx=zeros(n,1);
    for ip=1:n
        var_indx(ip)=find(strcmp(pad_hdr,var(ip))==1);
        lambda_indx(ip)=find(strcmp(pad_hdr,'wavelength'));
    end
    
    pad_dat=cell(m,1);
    pad_depth=cell(m,1);
    pad_station=cell(m,1);
    pad_lambda=cell(m,1);
    max_dat=0;
    min_dat=100;
    
    % Loop through variables to plot multiple if desired
    for i=1:n
        for irep=1:length(indx1)
            j=indx1(irep);
            
            %Clean data
            %bad_indx=find(all_pad.data{j}(:,var_indx(1))<all_pad.data{j}(:,var_indx(6)));
            %all_pad.data{j}(bad_indx,var_indx)=NaN; 
            
            pad_depth(i,1)={all_pad.depth{j}(:,1)};
            pad_dat(i,1)={all_pad.data{j}(:,var_indx(i))};
            pad_lambda(i,1)={all_pad.data{j}(:,lambda_indx(i))};
            pad_station(i,1)=strrep(strrep(all_pad.name(j),file_pref,''),'_1m.dat','');
            
            %Find max for plotting
            new_max_dat=max(all_pad.data{j}(:,var_indx(i)),[],1);
            new_min_dat=min(all_pad.data{j}(:,var_indx(i)),[],1);
            if new_max_dat > max_dat
                max_dat=new_max_dat;
            end
            if new_min_dat < min_dat
                min_dat=new_min_dat;
            end
        end
    end
    
    ylab='Filterpad Absorption (m^{-1})';
    %Override selection of min and max
    %min_dat=0;
    %max_dat=1;
    %O
    hplot=[];
    figure(1);
    clf
    for k=1:n
        hplot(k)=plot(pad_lambda{k},pad_dat{k},'-','Linewidth',1.2);
    
        if k==1
            set(gca,'Linewidth',1.3,'fontname','times','fontweight','bold','fontsize',16,...
            'xlim',[300 800],'xtick',300:100:800,'ylim',[min_dat-0.01.*abs(min_dat) max_dat+0.1.*max_dat]);
            xlabel('Wavelength (nm)','fontname','times','fontweight','bold','fontsize',16);
            ylabel(ylab,'fontname','times','fontweight','bold','fontsize',16);
            
            pad_data=pad_lambda{k};
    
         end
        
        set(hplot(k),'Linewidth',1.5);
        
        hold on
        
        pad_data=[pad_data,pad_dat{k}];
    end
    
    
    hleg=legend(hplot,strrep(var_str,'_',''));
    set(hleg,'Fontsize',10);
    
    stn_indx = strfind(all_pad.name{indx1},'_St');
    if ~isempty(stn_indx)
        sta_ttl = all_pad.name{indx1}(stn_indx+3:end - 6);
        text(500,max_dat*.9,['Stn ',strrep(sta_ttl,'_',' ')],'fontname','times','fontweight','bold','fontsize',14);
    else
        uw_indx = strfind(all_pad.name{indx1},'_UWS');
        sta_ttl = all_pad.name{indx1}(uw_indx+4:end - 6);
        text(500,max_dat.*.9,['UWS ',strrep(sta_ttl,'_',' ')],'fontname','times','fontweight','bold','fontsize',14);
    end
    
    %text(375,.2,'Shelf','fontname','times','fontweight','bold','fontsize',[18]);
    hold on
    fig_pos=get(gca,'position');
    fig_pos(2)=fig_pos(2)-0.0;
    set(gca,'Position',fig_pos);

    if ~run_flag
        continue
    end

    % save([inpath,file_pref,sta_txt,'_',var_str,'.mat pad_data']);
end

disp('Completed');
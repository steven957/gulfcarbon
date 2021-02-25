%Plot hyper spectra

clear all
clc
clf all reset 

%List of station names to select
station_names={'A1','A1b','A2','A3','B4','C1','C5','E0','E1','E3','D5','F4','F5','MR2','MR3',...
'G2','G3','G4','G5','H1','H2','H3','UWS1','UWS2','Run52'};

sta_txt='Run52';
file_pref='GulfCarbon2_gulfcarbon2_hypersas_St';

num_char=length(file_pref)+length(sta_txt);

inpath='C:\Users\slohrenz\Documents\Steve\DATA\NSF\GulfCarbon\GC_hypersensor_data\GC2_Hyper_sensors\GC2_HyperSAS_data\Output\To_SeaBASS\';
infile=[inpath,'all_hyper_April2009.mat'];

load(infile);

indx1=find(datenum(all_hyper.start_date)>datenum('12/31/08') & ...
    contains(all_hyper.name,sta_txt)); 
[m,n]=size(indx1);

hyper_hdr=all_hyper.fields{1};
hyper_units=all_hyper.units{1};

%Variable options are 'Lt','Lsky','Lw','Es' and 'Rrs'
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
        
var_indx=find(strncmp(hyper_hdr,var,var_char_num)==1); 
hyper_lambda=str2double(strrep(hyper_hdr(var_indx),var,''));

hyper_dat=cell(m,1);
hyper_depth=cell(m,1);
hyper_station=cell(m,1);
max_dat=0;
min_dat=100;

for ih=1:m
    j=indx1(ih);
    
    hyper_dat(ih,1)={all_hyper.data{j}(:,var_indx-2)}; %Adjust var_indx to account
    % for data and time fields not included in data array
    hyper_station(ih,1)=all_hyper.station(j);
    
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

hfig1=figure(1);
%Set up axes
hrrs=axes('YLimMode','auto','XTick',350:100:800,'Xlim',[340 810],'Position',[0.25 0.2 0.6 0.7]);

xtick_lab=get(gca,'Xtick');
set(hrrs,'Visible','on','Color','none','YAxisLocation','left','FontSize',14,...
    'FontWeight','normal','FontName','Arial','Linewidth',1.5,'XtickLabel',xtick_lab); 
box on
hold on

%Axis labels
hylab=ylabel(ylabel_txt,'FontSize',14,'FontWeight','bold','FontName','Arial',...
    'Rotation',90);
%ylab_pos=get(hylab,'Position');
%ylab_pos=[ylab_pos(1)+0.1.*time_range ylab_pos(2) ylab_pos(3)];
%set(hylab,'Position',ylab_pos);
xlabel('Wavelength(nm)','FontSize',14,'FontWeight','bold','FontName','Arial');
hold on

%Override selection of min and max if desired
%min_dat=0;
%max_dat=1;

for k=1:m
    hplot=plot(hyper_lambda,hyper_dat{k},'-','Linewidth',1.2);
    set(hplot,'Linewidth',1.5);
    hold on
end

%hleg=legend(hplot,strrep(var_str,'_',''));
%set(hleg,'Fontsize',10);

ylimits=get(gca,'YLim');
text(500,ylimits(2)+0.05.*ylimits(2),['Station ',strrep(sta_txt,'_',' ')],'fontname','times','fontweight','bold','fontsize',[18]);
hold on
fig_pos=get(gca,'position');
fig_pos(2)=fig_pos(2)-0.0;
set(gca,'Position',fig_pos);

%eval(['save ',inpath,'GC2_hypersas_',sta_txt,'_',var_str,'.mat hyper_data']);

disp('Completed');
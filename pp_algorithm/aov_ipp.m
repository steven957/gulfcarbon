%  aov_ipp.m
%
%   This script has two parts: 
%    Part-1: performs 1-way ANOVA for three models (WRM, WIM, WRME)
%    Part-2: perform M-way ANOVA

%  Syntax: aov_ipp

%  Inputs:
%    1) Spreadsheet with different model IPP results (e.g.,'all_IPP_tab.xlsx') 
%        with cruise, station, cluster, and IPP output
%
%  Output:
%    1) Part-1: 1- way ANOVA tables for WRM, WIM & WRME
%    2) Part-2: M-way ANOVA table

%  Other m-files required: None 
%  MAT-files required: None

%  Author: Israt Jahan Mili
%  School for Marine Science and Technology, University of Massachusetts Dartmouth
%  email address: mmili@umassd.edu
%  Website: https://www.umassd.edu/smast/
%  Last revision: 9 Dec 2025

%% ------------- BEGIN CODE --------------%
clc
clearvars

inpath='\';
xdata=readtable([inpath,'all_IPP_tab.xlsx'],'Sheet','Sheet1');
new_data = [xdata(:,["Cruise","Station","Cluster","WRM"]),repmat({"WRM"},size(xdata(:,"WRM"),1),1)];
new_data2 = [xdata(:,["Cruise","Station","Cluster","WIM"]),repmat({"WIM"},size(xdata(:,"WIM"),1),1)];
new_data3 = [xdata(:,["Cruise","Station","Cluster","WRME"]),repmat({"WRME"},size(xdata(:,"WRME"),1),1)];

new_data = renamevars(new_data,"WRM","IPP");
new_data2 = renamevars(new_data2,"WIM","IPP");
new_data3 = renamevars(new_data3,"WRME","IPP");

ipp_data = [new_data;new_data2;new_data3];
ipp_data = renamevars(ipp_data,"Var5","Model");

% 1st part: 1-way ANOVA for WRM, WIM, & WRME

% % Assign variables name
% wrm=xdata.WRM;
% wrme=xdata.WRME;
% wim=xdata.WIM;
% cluster=categorical(xdata.cluster);
% 
% % 1-way ANOVA_WRM, WIM, & WRME
% ipp_tab=table(wrm,wim,wrme,cluster);
% aov_wrm = anova(ipp_tab(:,[1,4]),'wrm');
% aov_wim = anova(ipp_tab(:,[2,4]),'wim');
% aov_wrme=anova(ipp_tab(:,[3,4]),'wrme');

% -------------------------------------%
% 2nd part: M-way ANOVA
% data=readtable([inpath,'all_IPP_tab.xlsx'],'Sheet','Sheet1');

% Convert 'Model' & 'Cluster' to categorical column
ipp_data.Model=categorical(ipp_data.Model);
ipp_data.Cluster=categorical(ipp_data.Cluster);
ipp_data.Cruise=categorical(ipp_data.Cruise);

% Defining table and selected column (If don't want to add 'Cruise' as
% variable then select (:,2:5))
tbl = ipp_data(:,1:5);

% If want variability of different clusters and different models
% aov=anova(tbl,'IPP ~ Model + Cluster + Model:Cluster');

% If want variability of different clusters and different models and
% different cruises

%aov=anova(tbl(tbl.Model~="WRME" & tbl.Cluster~="estuary" & tbl.Cluster~="outer",:),'IPP ~ Model'); 

p = signtest(tbl.IPP(tbl.Model=="WRM"),tbl.IPP(tbl.Model=="WIM"));
aov=anova(tbl,'IPP ~ Model + Cluster + Cruise + Model:Cruise:Cluster'); 
MBE = mean(tbl.IPP(tbl.Model=="WIM")-tbl.IPP(tbl.Model=="WRM"));

% If want variability of different clusters and different models and
% different cruises
aov2=anova(tbl,'IPP ~ Cruise + Cluster + Cluster:Cruise'); 
% aov2=anova(xdata,'WRM ~ Cruise + Cluster + Cluster:Cruise'); 

hbx = boxplot(xdata.WRM,xdata.Cluster);

%  ===========
%  END OF CODE  %
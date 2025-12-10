% pmax_aph440_vs_phi_fit.m
%
% Script plots relationship between maximum photosynthetic rate normalized
%   to phytoplankton absorption at 440 nm (Pmax_aph440_spec) vs maximum photosynthetic quantum
%   yield for carbon fixation (phi)
%
% Inputs:
%    1) Folder location with spreadsheet with phytoplankton absorption at 
%      440 nm (Pmax_aph440_spec) vs maximum photosynthetic quantum yield
%      for carbon fixation (phi)
%    2) Spreadsheet with data ('Pmax_aph_spec_vs_phi_table.xlsx') 
%
% Outputs:
%    1) Plot of Pmax_aph440_spec vs phi with correlation statistics
%   
% Other m-files required: None 
%
% MAT-files required: None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 9 Dec 2025

%% ------------- BEGIN CODE --------------%


input_folder = '\P-E\';

sheet_file = 'Pmax_aph_spec_vs_phi_table.xlsx';

pe_tab = readtable([input_folder,sheet_file]);

% Omit NaNs
Pmax_aph440_spec_2 = pe_tab.Pmax_aph440_spec(~isnan(pe_tab.Pmax_aph440_spec) & ~isnan(pe_tab.phi));
phi_2 = pe_tab.phi(~isnan(pe_tab.Pmax_aph440_spec) & ~isnan(pe_tab.phi));

% ft = fittype('a*x^b'); % Power fit
ft = fittype('poly2');  % Second order polynomial fit

% Model fit 
[curve,gof,output] = fit(Pmax_aph440_spec_2,phi_2,ft);

%% Plot curve
f1=figure(1);
clf

hp2 = plot(curve,Pmax_aph440_spec_2,phi_2,'o');
hold on
hp2(2).LineWidth = 1.5;
hp2(1).LineWidth = 1;

set(gca,'LineWidth',1.5,'Fontsize',12);
xlabel('{\itP_{max}^{aph440}est}','Fontsize',13,'FontWeight','bold');
ylabel('\bf{\phi}^{C}_{max}est','Fontsize',13,'FontWeight','bold');

% Annotation
fit_txt1 = ['{\phi^{C}_{max}est} = ',num2str(curve.p1,'%4.1f'),...
    '{\itP_{max}^{aph440}est}','^2 + '];
fit_txt2 = [num2str(curve.p2,'%4.2f'),'{\itP_{max}^{aph440}est} - ',num2str(-curve.p3,'%7.5f')];

txt1 = text(0.018,0.03,fit_txt1,"FontSize",10);
txt2 = text(0.025,0.02,fit_txt2,"FontSize",10);

txt3 = text(0.018,0.01,['N = ',num2str(output.numobs),'; {\itr^2} = ',num2str(gof.adjrsquare,'%5.3f'),...
    '; RMSE = ',num2str(gof.rmse,'%5.3f')],'FontSize',10);

% Print plot
print([input_folder,'pmax_aph440_vs_phi_curve_fit.tif'],'-dtiff','-r300');

disp('Completed');

% From earlier script, modified by Darren Pilcher 03/17/15
% Uses .mat output from calc_michigan_gpp_avg.m 

close all
clear

% Load output from each simulation 

% ********** Past Nutrients ***************

load('.../past/phygrow1_depth_int_FA_0710.mat');
model_phygrow1.control = phygrow1; clear phygrow1
load('.../past/phygrow2_depth_int_FA_0710.mat');
model_phygrow2.control = phygrow2; clear phygrow2
load('.../past/phy1_depth_int_FA_0710.mat');
model_phy1.control = phy1; clear phy1
load('.../past/phy2_depth_int_FA_0710.mat');
model_phy2.control = phy2; clear phy2
load('.../past/zoo_depth_int_FA_0710.mat');
model_zoo.control = zoo; clear zoo
load('.../past/lakewide_surf_temp_0710.mat');
model_surf_temp.control = lakewide_surf_temp; clear lakewide_surf_temp
load('.../past/lakewide_bottom_temp_0710.mat');
model_bottom_temp.control = lakewide_bottom_temp; clear lakewide_bottom_temp


% ********** Present Nutrients **************

load('.../present/phygrow1_depth_int_FA_0710.mat');
model_phygrow1.pres = phygrow1; clear phygrow1
load('.../present/phygrow2_depth_int_FA_0710.mat');
model_phygrow2.pres = phygrow2; clear phygrow2
load('.../present/phy1_depth_int_FA_0710.mat');
model_phy1.pres = phy1; clear phy1
load('.../present/phy2_depth_int_FA_0710.mat');
model_phy2.pres = phy2; clear phy2
load('.../present/zoo_depth_int_FA_0710.mat');
model_zoo.pres = zoo; clear zoo



% *********** Quagga Present Nutrients ******************

load('.../quagga_present/phygrow1_depth_int_FA_0710.mat');
model_phygrow1.qpres = phygrow1; clear phygrow1
load('.../quagga_present/phygrow2_depth_int_FA_0710.mat');
model_phygrow2.qpres = phygrow2; clear phygrow2
load('.../quagga_present/phy1_depth_int_FA_0710.mat');
model_phy1.qpres = phy1; clear phy1
load('.../quagga_present/phy2_depth_int_FA_0710.mat');
model_phy2.qpres = phy2; clear phy2
load('.../quagga_present/zoo_depth_int_FA_0710.mat');
model_zoo.qpres = zoo; clear zoo

% ************ Quagga Past Nutrients ****************

load('.../quagga_past/phygrow1_depth_int_FA_0710.mat');
model_phygrow1.qpast = phygrow1; clear phygrow1
load('.../quagga_past/phygrow2_depth_int_FA_0710.mat');
model_phygrow2.qpast = phygrow2; clear phygrow2
load('.../quagga_past/phy1_depth_int_FA_0710.mat');
model_phy1.qpast = phy1; clear phy1
load('.../quagga_past/phy2_depth_int_FA_0710.mat');
model_phy2.qpast = phy2; clear phy2
load('.../quagga_past/zoo_depth_int_FA_0710.mat');
model_zoo.qpast = zoo; clear zoo



% For GPP, integrate with depth, and convert units to mgC/m^2/day and find mean across both stations 
model_npp.control = squeeze(mean(((squeeze(sum(model_phygrow1.control,3)) + squeeze(sum(model_phygrow2.control,3))) .* (24*60*60) .* (200) .* (1/1000) .* (12.01) .* 1000),2)); 
model_npp.pres = squeeze(mean(((squeeze(sum(model_phygrow1.pres,3)) + squeeze(sum(model_phygrow2.pres,3))) .* (24*60*60) .* (200) .* (1/1000) .* (12.01) .* 1000),2));  
model_npp.qpres = squeeze(mean(((squeeze(sum(model_phygrow1.qpres,3)) + squeeze(sum(model_phygrow2.qpres,3))) .* (24*60*60) .* (200) .* (1/1000) .* (12.01) .* 1000),2)); 
model_npp.qpast = squeeze(mean(((squeeze(sum(model_phygrow1.qpast,3)) + squeeze(sum(model_phygrow2.qpast,3))) .* (24*60*60) .* (200) .* (1/1000) .* (12.01) .* 1000),2)); 

% Take mean over desired years and multiply by model NPP to Fahnenstiel GPP conversion factor of 1.4 
model_gpp.control = squeeze(mean(model_npp.control(1:2,:),1)) * 1.4;
model_gpp.pres = squeeze(mean(model_npp.pres(1:2,:),1)) * 1.4;
model_gpp.qpres = squeeze(mean(model_npp.qpres(1:2,:),1)) * 1.4;
model_gpp.qpast = squeeze(mean(model_npp.qpast(1:2,:),1)) * 1.4;

% For plotting purpose, make a version of model_gpp.qpres that is based on the high and low GPP:NPP conversion ratios 
model_gpp.qpres_low = squeeze(mean(model_npp.qpres(1:2,:),1)) * 1.18;
model_gpp.qpres_high = squeeze(mean(model_npp.qpres(1:2,:),1)) * 1.86;

model_std.qpres = std((model_npp.qpres(1:2,:) * 1.4),0,1);

% Take model monthly average
mon = [1 32 60 91 121 152 182 213 244 274 305 336 366];
for m = 1:12 
	model_gpp_mon(1,m) = mean(model_gpp.control(mon(m):mon(m+1)-1)); 
	model_gpp_mon(2,m) = mean(model_gpp.qpres(mon(m):mon(m+1)-1)); 

	model_gpp_mon_alt(1,m) = mean(model_gpp.pres(mon(m):mon(m+1)-1));
        model_gpp_mon_alt(2,m) = mean(model_gpp.qpast(mon(m):mon(m+1)-1)); 

% This complicated line is finding the standard deviation as the interannual variability of the monthly mean model GPP 
	model_gpp_mon_std(1,m) = std(mean((model_npp.control(1:2,mon(m):mon(m+1)-1) * 1.4),2),0,1);
	model_gpp_mon_std(2,m) = std(mean((model_npp.qpres(1:2,mon(m):mon(m+1)-1) * 1.4),2),0,1);

	model_gpp_mon_std_alt(1,m) = std(mean((model_npp.pres(1:2,mon(m):mon(m+1)-1) * 1.4),2),0,1);
        model_gpp_mon_std_alt(2,m) = std(mean((model_npp.qpast(1:2,mon(m):mon(m+1)-1) * 1.4),2),0,1);
end  

% Now load annual depth integrated NPP in units of mmolP/m^2/sec and convert to units of mgC/m^2/day 
load('.../lakewide.depth.integrated.npp.past.mat');
lakewide_npp.control = npp * (24*60*60)*(200)*(12.01); clear npp 
load('.../lakewide.depth.integrated.npp.pres.mat');
lakewide_npp.pres = npp * (24*60*60)*(200)*(12.01); clear npp
load('.../lakewide.depth.integrated.npp.qpast.mat');
lakewide_npp.qpast = npp * (24*60*60)*(200)*(12.01); clear npp 
load('.../lakewide.depth.integrated.npp.qpres.mat');
lakewide_npp.qpres = npp * (24*60*60)*(200)*(12.01); clear npp 

for yr = 1:4 
	lakewide_avg_npp.control(yr,:) = lake_surface_avg_michigan(squeeze(lakewide_npp.control(yr,:,:,:)),2); 
	lakewide_avg_npp.pres(yr,:) = lake_surface_avg_michigan(squeeze(lakewide_npp.pres(yr,:,:,:)),2); 
	lakewide_avg_npp.qpast(yr,:) = lake_surface_avg_michigan(squeeze(lakewide_npp.qpast(yr,:,:,:)),2); 
	lakewide_avg_npp.qpres(yr,:) = lake_surface_avg_michigan(squeeze(lakewide_npp.qpres(yr,:,:,:)),2); 
end


% Estimates of PP from Fahenstiel 2010 Fig. 2  

fah_pp_past = [nan nan nan nan 870 710 950 600 nan nan nan nan; nan nan nan nan 1160 1020 1200 690 nan nan nan nan; nan nan 475 1110 1160 990 990 440 nan 405 nan nan; nan nan 490 925 1315 725 1200 735 nan 500 315 nan;...
    nan nan 375 825 1450 1550 1250 590 490 850 300 190; nan nan nan 715 1100 1190 1140 1190 300 225 205 235; nan nan 325 720 1035 1000 890 635 590 525 250 150; nan nan 365 525 965 1200 910 nan nan nan nan nan];
fah_pp_past_mean = nanmean(fah_pp_past,1);
fah_pp_past_std = nanstd(fah_pp_past,0,1); 

fah_pp = [nan nan 215 210 200 675 775 635 435 315 300 175; nan nan 215 250 235 440 925 645 470 405 325 nan]; 
fah_pp_mean = nanmean(fah_pp,1);
fah_pp_std = nanstd(fah_pp,0,1);
fah_pp_plot = nan(1,365);
fah_pp_plot_std = nan(1,365);

fah_pp_plot(75) = fah_pp_mean(3); fah_pp_plot(106) = fah_pp_mean(4); fah_pp_plot(136) = fah_pp_mean(5); fah_pp_plot(167) = fah_pp_mean(6); fah_pp_plot(197) = fah_pp_mean(7);
fah_pp_plot(228) = fah_pp_mean(8); fah_pp_plot(259) = fah_pp_mean(9); fah_pp_plot(289) = fah_pp_mean(10); fah_pp_plot(320) = fah_pp_mean(11); fah_pp_plot(351) = fah_pp_mean(12);

fah_pp_plot_std(75) = fah_pp_std(3); fah_pp_plot_std(106) = fah_pp_std(4); fah_pp_plot_std(136) = fah_pp_std(5); fah_pp_plot_std(167) = fah_pp_std(6); fah_pp_plot_std(197) = fah_pp_std(7);
fah_pp_plot_std(228) = fah_pp_std(8); fah_pp_plot_std(259) = fah_pp_std(9); fah_pp_plot_std(289) = fah_pp_std(10); fah_pp_plot_std(320) = fah_pp_std(11); fah_pp_plot_std(351) = fah_pp_std(12);

fah_to_interp = [fah_pp_mean(end) fah_pp_mean(3:end) fah_pp_mean(3)];
fah_pp_interp = interp1([1 89 120 150 181 211 242 273 303 334 365 454],fah_to_interp,[16:380]);

plot_gpp = [fah_pp_past_mean; fah_pp_mean; model_gpp_mon];
plot_gpp_err = [fah_pp_past_std; fah_pp_std; model_gpp_mon_std]; 

plot_gpp_alt = [fah_pp_mean;  model_gpp_mon_alt];
plot_gpp_err_alt = [fah_pp_std; model_gpp_mon_std_alt];

%%

% Now to plot up the results
x = [1:365];

% Compute RMSE values for pre-quagga and post-quagga Fahnenstiel comparisons 
% RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error
RMSE_gpp_past = sqrt(mean((plot_gpp(3,3:end) - plot_gpp(1,3:end)).^2)); 
RMSE_gpp = sqrt(mean((plot_gpp(4,3:end) - plot_gpp(2,3:end)).^2));

RMSE_gpp_pres = sqrt(mean((plot_gpp_alt(2,3:end) - plot_gpp_alt(1,3:end)).^2));
RMSE_gpp_qpast = sqrt(mean((plot_gpp_alt(3,3:end) - plot_gpp_alt(1,3:end)).^2));

% Plot figures
figure(1);
set(gcf, 'units','centimeters','position',[10 10 19 12]);
set(gcf,'Color',[1 1 1])
subplot(2,1,1)
[h,herr] = barwitherr(plot_gpp_err([1 3],:)',[1:12],plot_gpp([1 3],:)');
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
axis([ 0 13 0 1750])
ylabel('GPP (mgC/m^2/day)','Fontsize',10)
set(gca,'xticklabel',{,'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
title('"Past Nutrients" Productivity','Fontsize',12)
h = legend('Observed','Model');
set(h,'Fontsize',10)

subplot(2,1,2)
[h,herr] = barwitherr(plot_gpp_err([2 4],:)',[1:12],plot_gpp([2 4],:)');
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
axis([ 0 13 0 1750])
ylabel('GPP (mgC/m^2/day)','Fontsize',10)
set(gca,'xticklabel',{,'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
title('"Quagga Present Nutrients" Productivity','Fontsize',12)
h = legend('Observed','Model');
set(h,'Fontsize',10)

figure(2); 
set(gcf, 'units','centimeters','position',[10 10 19 14]);
set(gcf,'Color',[1 1 1])
subplot(2,1,1)
plot(nanmean(lakewide_avg_npp.control,1),'k','Linewidth',2)
hold on
plot(nanmean(lakewide_avg_npp.pres,1),'r','Linewidth',2)
plot(nanmean(lakewide_avg_npp.qpast,1),'c','Linewidth',2);
plot(nanmean(lakewide_avg_npp.qpres,1),'b','Linewidth',2)
title('Lakewide Averaged NPP 2007-2010','Fontsize',12)
ylabel('NPP (mgC/m^2/day)','Fontsize',10)
set(gca,'xtick',[1 32 60 91 121 152 182 213 244 274 305 336]);
set(gca,'xticklabel',{'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
axis([0 365 0 1500])
h = legend('Past Nutrients','Present Nutrients','Quagga Past Nutrients','Quagga Present Nutrients');
set(h,'FontSize',10);

subplot(2,1,2) 
plot(model_surf_temp.control,'k','Linewidth',2)
hold on
plot(model_bottom_temp.control,'k--','Linewidth',2)
title('Lakewide Averaged Temperature 2007-2010','Fontsize',12)
ylabel('Temperature ({\circ}C)','Fontsize',10)
set(gca,'xtick',[1 32 60 91 121 152 182 213 244 274 305 336]);
set(gca,'xticklabel',{'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
axis([0 365 0 25])
h = legend('Surface','Bottom');
set(h,'FontSize',10);






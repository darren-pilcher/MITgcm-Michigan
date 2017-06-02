% Script used to generate the 25m average variables of NPP, PHY, and ZOO to use for ...
% plotting in the plot_depth_int_npp_zoo.m

% * NOTE: First need to generate 25m average files using calc_lakewide_var.m  
% Created by Darren Pilcher 05/2017

ncload('.../LMgrid1.kronos.nc','HFacC','X','Y','Z','Zl')

% Set the C:P ratio used for converting NPP,PHY, and ZOO into units of carbon
CP1 = 200;
CP2 = 200;

% Convert phy productivity from mmolP/m^2/s to mgC/m^2/day 
load('.../lakewide.25m.avg.npp1.past.mat');
lakewide_npp1.control = npp1 * (24*60*60)*(CP1)*(12.01); clear npp1
load('.../lakewide.25m.avg.npp1.pres.mat');
lakewide_npp1.pres = npp1 * (24*60*60)*(CP1)*(12.01); clear npp1
load('.../lakewide.25m.avg.npp1.qpast.mat');
lakewide_npp1.qpast = npp1 * (24*60*60)*(CP1)*(12.01); clear npp1
load('.../lakewide.25m.avg.npp1.qpres.mat');
lakewide_npp1.qpres = npp1 * (24*60*60)*(CP1)*(12.01); clear npp1

load('.../lakewide.25m.avg.npp2.past.mat');
lakewide_npp2.control = npp2 * (24*60*60)*(CP2)*(12.01); clear npp2
load('.../lakewide.25m.avg.npp2.pres.mat');
lakewide_npp2.pres = npp2 * (24*60*60)*(CP2)*(12.01); clear npp2
load('.../lakewide.25m.avg.npp2.qpast.mat');
lakewide_npp2.qpast = npp2 * (24*60*60)*(CP2)*(12.01); clear npp2
load('.../lakewide.25m.avg.npp2.qpres.mat');
lakewide_npp2.qpres = npp2 * (24*60*60)*(CP2)*(12.01); clear npp2

% Now also load annual 25m average Phy1 in units of mmolP/m^2 and convert to units of mgC/m^2
load('.../lakewide.25m.avg.phy1.past.mat');
lakewide_phy1.control = phy1 * CP1 * 12.01; clear phy1
load('.../lakewide.25m.avg.phy1.pres.mat');
lakewide_phy1.pres = phy1 * CP1 * 12.01; clear phy1
load('.../lakewide.25m.avg.phy1.qpast.mat');
lakewide_phy1.qpast = phy1 * CP1 * 12.01; clear phy1
load('.../lakewide.25m.avg.phy1.qpres.mat');
lakewide_phy1.qpres = phy1 * CP1 * 12.01; clear phy1

% Now also load annual 25m average Phy2 in units of mmolP/m^2 and convert to units of mgC/m^2
load('.../lakewide.25m.avg.phy2.past.mat');
lakewide_phy2.control = phy2 * CP2 * 12.01; clear phy2
load('.../lakewide.25m.avg.phy2.pres.mat');
lakewide_phy2.pres = phy2 * CP2 * 12.01; clear phy2
load('.../lakewide.25m.avg.phy2.qpast.mat');
lakewide_phy2.qpast = phy2 * CP2 * 12.01; clear phy2
load('.../lakewide.25m.avg.phy2.qpres.mat');
lakewide_phy2.qpres = phy2 * CP2 * 12.01; clear phy2

% Now also load annual 25m average Zoo in units of mmolP/m^2 and convert to units of mgC/m^2
load('.../lakewide.25m.avg.zoo.past.mat');
lakewide_zoo.control = zoo * 200 * 12.01; clear zoo
load('.../lakewide.25m.avg.zoo.pres.mat');
lakewide_zoo.pres = zoo * 200 * 12.01; clear zoo
load('.../lakewide.25m.avg.zoo.qpast.mat');
lakewide_zoo.qpast = zoo * 200 * 12.01; clear zoo
load('.../lakewide.25m.avg.zoo.qpres.mat');
lakewide_zoo.qpres = zoo * 200 * 12.01; clear zoo

% Take 4 year mean for each
lakewide_mean_npp1.control = squeeze(nanmean(lakewide_npp1.control,1));
lakewide_mean_npp1.pres = squeeze(nanmean(lakewide_npp1.pres,1));
lakewide_mean_npp1.qpast = squeeze(nanmean(lakewide_npp1.qpast,1));
lakewide_mean_npp1.qpres = squeeze(nanmean(lakewide_npp1.qpres,1));

lakewide_mean_npp2.control = squeeze(nanmean(lakewide_npp2.control,1));
lakewide_mean_npp2.pres = squeeze(nanmean(lakewide_npp2.pres,1));
lakewide_mean_npp2.qpast = squeeze(nanmean(lakewide_npp2.qpast,1));
lakewide_mean_npp2.qpres = squeeze(nanmean(lakewide_npp2.qpres,1));

lakewide_mean_phy1.control = squeeze(nanmean(lakewide_phy1.control,1));
lakewide_mean_phy1.pres = squeeze(nanmean(lakewide_phy1.pres,1));
lakewide_mean_phy1.qpast = squeeze(nanmean(lakewide_phy1.qpast,1));
lakewide_mean_phy1.qpres = squeeze(nanmean(lakewide_phy1.qpres,1));

lakewide_mean_phy2.control = squeeze(nanmean(lakewide_phy2.control,1));
lakewide_mean_phy2.pres = squeeze(nanmean(lakewide_phy2.pres,1));
lakewide_mean_phy2.qpast = squeeze(nanmean(lakewide_phy2.qpast,1));
lakewide_mean_phy2.qpres = squeeze(nanmean(lakewide_phy2.qpres,1));

lakewide_mean_zoo.control = squeeze(nanmean(lakewide_zoo.control,1));
lakewide_mean_zoo.pres = squeeze(nanmean(lakewide_zoo.pres,1));
lakewide_mean_zoo.qpast = squeeze(nanmean(lakewide_zoo.qpast,1));
lakewide_mean_zoo.qpres = squeeze(nanmean(lakewide_zoo.qpres,1));

% Save variables
save(['.../lakewide.25m.mean.npp1.CP',num2str(CP1),'.mat'],'lakewide_mean_npp1')
save(['.../lakewide.25m.mean.npp2.CP',num2str(CP2),'.mat'],'lakewide_mean_npp2')
save(['.../lakewide.25m.mean.phy1.CP',num2str(CP1),'.mat'],'lakewide_mean_phy1')
save(['.../lakewide.25m.mean.phy2.CP',num2str(CP2),'.mat'],'lakewide_mean_phy2')
save('.../lakewide.25m.mean.zoo.mat','lakewide_mean_zoo')











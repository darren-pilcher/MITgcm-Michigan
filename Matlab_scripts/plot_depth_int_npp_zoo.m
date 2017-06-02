% Script used to make spatial plots of depth averaged NPP, PHY, and ZOO populations 
% Uses files that are generated in the calc_lakewide_var.m script

% Created by Darren Pilcher 04/24/15

close all
clear

ncload('.../LMgrid1.kronos.nc','HFacC','X','Y','Z','Zl','Depth')

CP1 = 200; 
CP2 = 200; 

load(['.../lakewide.25m.mean.npp1.CP',num2str(CP1),'.mat',]'lakewide_mean_npp1')
load(['.../lakewide.25m.mean.npp2.CP',num2str(CP2),'.mat',]'lakewide_mean_npp2')
load(['.../lakewide.25m.mean.phy1.CP',num2str(CP1),'.mat',]'lakewide_mean_phy1')
load(['.../lakewide.25m.mean.phy2.CP',num2str(CP2),'.mat',]'lakewide_mean_phy2')
load('.../lakewide.25m.mean.zoo.mat','lakewide_mean_zoo')

lakewide_mean_npp.control = lakewide_mean_npp1.control + lakewide_mean_npp2.control;
lakewide_mean_npp.pres = lakewide_mean_npp1.pres + lakewide_mean_npp2.pres;
lakewide_mean_npp.qpast = lakewide_mean_npp1.qpast + lakewide_mean_npp2.qpast;
lakewide_mean_npp.qpres = lakewide_mean_npp1.qpres + lakewide_mean_npp2.qpres;

lakewide_mean_total_phy.control = lakewide_mean_phy1.control + lakewide_mean_phy2.control; 
lakewide_mean_total_phy.pres = lakewide_mean_phy1.pres + lakewide_mean_phy2.pres; 
lakewide_mean_total_phy.qpast = lakewide_mean_phy1.qpast + lakewide_mean_phy2.qpast; 
lakewide_mean_total_phy.qpres = lakewide_mean_phy1.qpres + lakewide_mean_phy2.qpres; 

% Calculate spatial correlations with depth 
H = squeeze(HFacC(1,:,:));
fih = find(H==0);
Depth(fih) = nan; 

%%

% Make a plot in Color of Control NPP, PHY, and ZOO 
figure(1); 
set(gcf, 'units','centimeters','position',[10 10 19 11.5]);
set(gcf,'Color',[1 1 1])
subplot('position',[0.1 0.1 0.25 0.8])
pcolor(X,Y,squeeze(nanmean(lakewide_mean_npp.control,3)))
shading flat 
colormap(jet)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
ylabel('Past Nutrients','Fontsize',14)
h = colorbar('Location','SouthOutside');
caxis([0 150]) 
xlabel(h,'NPP (mgC/m^2/d)','Fontsize',10)
set(h,'xtick',[0:50:150])
title('NPP ','Fontsize',12)

subplot('position',[0.4 0.1 0.25 0.8])
pcolor(X,Y,squeeze(nanmean(lakewide_mean_total_phy.control,3)))
shading flat 
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
h = colorbar('Location','SouthOutside');
caxis([0 650]) 
xlabel(h,'Total PHY (mgC/m^2)','Fontsize',10)
set(h,'xtick',[0:200:600])
title('PHY','Fontsize',12)

subplot('position',[0.7 0.1 0.25 0.8])
pcolor(X,Y,squeeze(nanmean(lakewide_mean_zoo.control,3)))
shading flat 
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
h = colorbar('Location','SouthOutside');
caxis([0 300]) 
xlabel(h,'ZOO (mgC/m^2)','Fontsize',10)
set(h,'xtick',[0:100:300])
title('ZOO','Fontsize',12)


figure(2);  
set(gcf, 'units','centimeters','position',[10 10 19 11.5]);
set(gcf,'Color',[1 1 1])
subplot('position',[0.1 0.1 0.25 0.8])
pcolor(X,Y,squeeze(nanmean(lakewide_mean_npp.qpres,3)))
shading flat 
colormap(jet)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
ylabel('Quagga Present Nutrients','Fontsize',14)
h = colorbar('Location','SouthOutside');
caxis([0 150]) 
xlabel(h,'NPP (mgC/m^2/d)','Fontsize',10)
set(h,'xtick',[0:50:150])
title('NPP','Fontsize',12)

subplot('position',[0.4 0.1 0.25 0.8])
pcolor(X,Y,squeeze(nanmean(lakewide_mean_total_phy.qpres,3)))
shading flat 
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
h = colorbar('Location','SouthOutside');
caxis([0 650]) 
xlabel(h,'Total PHY (mgC/m^2)','Fontsize',10)
set(h,'xtick',[0:200:600])
title('PHY','Fontsize',12)

subplot('position',[0.7 0.1 0.25 0.8])
pcolor(X,Y,squeeze(nanmean(lakewide_mean_zoo.qpres,3)))
shading flat 
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
h = colorbar('Location','SouthOutside');
caxis([0 300]) 
xlabel(h,'ZOO (mgC/m^2)','Fontsize',10)
set(h,'xtick',[0:100:300])
title('ZOO','Fontsize',12)

% COMBINE ALL DIFFERENCE PLOTS INTO A SINGLE FIGURE 

V = [-1:.1:.2];
cm = [127/255 0 1; 76/255 0 153/255; 0 0 153/255; 0 0 204/255; 0 0 1; 0 102/255 204/255; 0 128/255 1; 51/255 153/255 1; 102/255 178/255 1; 204/255 229/255 1]; 
cm = [cm; 1 153/255 153/255; 1 51/255 51/255];
figure(6);
set(gcf, 'units','centimeters','position',[10 10 19 23]);
set(gcf,'Color',[1 1 1])

subplot('position',[0.1 0.675 0.2 0.275]) % Panel 1
[cs,h] = contourf(X,Y,((squeeze(nanmean(lakewide_mean_npp.pres,3)) - squeeze(nanmean(lakewide_mean_npp.control,3))) ./ squeeze(nanmean(lakewide_mean_npp.control,3))),V);
colormap(cm);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'','','',''},'Fontsize',8)
ylabel('Present Nutrients','Fontsize',12,'fontweight','bold')
caxis([-1 .2])
title('NPP','Fontsize',16)

subplot('position',[0.375 0.675 0.2 0.275]) % Panel 2
[cs,h] = contourf(X,Y,((squeeze(nanmean(lakewide_mean_total_phy.pres,3)) - squeeze(nanmean(lakewide_mean_total_phy.control,3))) ./ squeeze(nanmean(lakewide_mean_total_phy.control,3))),V);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'','','','',''},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'','','',''},'Fontsize',8)
caxis([-1 .2])
title('PHY','Fontsize',16)

subplot('position',[0.65 0.675 0.2 0.275]) % Panel 3
[cs,h] = contourf(X,Y,((squeeze(nanmean(lakewide_mean_zoo.pres,3)) - squeeze(nanmean(lakewide_mean_zoo.control,3))) ./ squeeze(nanmean(lakewide_mean_zoo.control,3))),V);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'','','','',''},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'','','',''},'Fontsize',8)
caxis([-1 .2])
title('ZOO','Fontsize',16)

axes('Position', [0.65 0.675 0.3 0.275], 'Visible', 'off');
h = colorbar('EastOutside','Fontsize',8);
caxis([-1 .2])
ylabel(h,'% Change','Fontsize',10)
set(h,'ytick',[-1:.1:.2])
set(h,'yticklabel',{'-100','','-80','','-60','','-40','','-20','','0','','20'})

subplot('position',[0.1 0.375 0.2 0.275]) % Panel 4
contourf(X,Y,((squeeze(nanmean(lakewide_mean_npp.qpast,3)) - squeeze(nanmean(lakewide_mean_npp.control,3))) ./ squeeze(nanmean(lakewide_mean_npp.control,3))),V)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'','','',''},'Fontsize',8)
ylabel('Quagga Past Nutrients','Fontsize',12,'fontweight','bold')
caxis([-1 .2])

subplot('position',[0.375 0.375 0.2 0.275]) % Panel 5
contourf(X,Y,((squeeze(nanmean(lakewide_mean_total_phy.qpast,3)) - squeeze(nanmean(lakewide_mean_total_phy.control,3))) ./ squeeze(nanmean(lakewide_mean_total_phy.control,3))),V)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'','','','',''},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'','','',''},'Fontsize',8)
caxis([-1 .2])

subplot('position',[0.65 0.375 0.2 0.275]) % Panel 6
contourf(X,Y,((squeeze(nanmean(lakewide_mean_zoo.qpast,3)) - squeeze(nanmean(lakewide_mean_zoo.control,3))) ./ squeeze(nanmean(lakewide_mean_zoo.control,3))),V)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'','','','',''},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'','','',''},'Fontsize',8)
caxis([-1 .2])

axes('Position', [0.65 0.375 0.3 0.275], 'Visible', 'off');
h = colorbar('EastOutside','Fontsize',8);
caxis([-1 .2])
ylabel(h,'% Change','Fontsize',10)
set(h,'ytick',[-1:.1:.2])
set(h,'yticklabel',{'-100','','-80','','-60','','-40','','-20','','0','','20'})

subplot('position',[0.1 0.075 0.2 0.275]) % Panel 7
contourf(X,Y,((squeeze(nanmean(lakewide_mean_npp.qpres,3)) - squeeze(nanmean(lakewide_mean_npp.control,3))) ./ squeeze(nanmean(lakewide_mean_npp.control,3))),V)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
ylabel('Quagga Present Nutrients','Fontsize',12,'fontweight','bold')
caxis([-1 .2])

subplot('position',[0.375 0.075 0.2 0.275]) % Panel 8
contourf(X,Y,((squeeze(nanmean(lakewide_mean_total_phy.qpres,3)) - squeeze(nanmean(lakewide_mean_total_phy.control,3))) ./ squeeze(nanmean(lakewide_mean_total_phy.control,3))),V)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'','','','',''},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
caxis([-1 .2])

subplot('position',[0.65 0.075 0.2 0.275]) % Panel 9
contourf(X,Y,((squeeze(nanmean(lakewide_mean_zoo.qpres,3)) - squeeze(nanmean(lakewide_mean_zoo.control,3))) ./ squeeze(nanmean(lakewide_mean_zoo.control,3))),V)
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'','','','',''},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
caxis([-1 .2])

axes('Position', [0.65 0.075 0.3 0.275], 'Visible', 'off');
h = colorbar('EastOutside','Fontsize',8);
caxis([-1 .2])
ylabel(h,'% Change','Fontsize',10)
set(h,'ytick',[-1:.1:.2])
set(h,'yticklabel',{'-100','','-80','','-60','','-40','','-20','','0','','20'})









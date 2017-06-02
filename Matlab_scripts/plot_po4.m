% Short script to make spatial plots of PO4 in the different model simulations, to illustrate that with quaggas, 
% PO4 is still in excess in nearshore during winter, suggesting that river tributaries won't matter then. 

% Created by Darren Pilcher 01/2017

close all 
clear 

DIR1 = '.../present/';
DIR2 = '.../quagga_present/'; 

ncload('.../LMgrid1.kronos.nc','HFacC','X','Y');
H = squeeze(HFacC(1,:,:)); clear HFacC
fi = find(H==0);

for yr = 2007:2010
        ncload([DIR1,num2str(yr),'/po4.',num2str(yr),'.nc']);
	phos.pres(yr-2006,:,:,:) = squeeze(po4(:,1,:,:)); clear po4 

	ncload([DIR2,num2str(yr),'/po4.',num2str(yr),'.nc']);
        phos.qpres(yr-2006,:,:,:) = squeeze(po4(:,1,:,:)); clear po4
end

% Take average over all 4 years

mean_phos.pres = squeeze(nanmean(phos.pres,1)); 
mean_phos.qpres = squeeze(nanmean(phos.qpres,1)); 

% Mask out land 
for t = 1:size(mean_phos.pres,1)
	tmp = squeeze(mean_phos.pres(t,:,:)); 
	tmp(fi) = nan; 
	mean_phos.pres(t,:,:) = tmp; clear tmp 

        tmp = squeeze(mean_phos.qpres(t,:,:)); 
        tmp(fi) = nan; 
        mean_phos.qpres(t,:,:) = tmp; clear tmp
end


% Generate Figure
V = [0:0.02:.2];
figure(1) 
set(gcf, 'units','centimeters','position',[10 10 19 11.5]);
set(gcf,'Color',[1 1 1])

subplot('position',[0.1 0.55 0.15 0.375]) % Panel 1
%[cs, h] = contourf(X,Y,squeeze(nanmean(mean_phos.pres([1:59 336:end],:,:),1)),V);
pcolor(X,Y,squeeze(nanmean(mean_phos.pres([1:59 336:end],:,:),1))); 
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('DJF','Fontsize',12)
caxis([0 .2])

subplot('position',[0.3 0.55 0.15 0.375]) % Panel 2
pcolor(X,Y,squeeze(nanmean(mean_phos.pres(60:151,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('MAM','Fontsize',12)
caxis([0 .2])

subplot('position',[0.5 0.55 0.15 0.375]) % Panel 3
pcolor(X,Y,squeeze(nanmean(mean_phos.pres(152:243,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('JJA','Fontsize',12)
caxis([0 .2])

subplot('position',[0.7 0.55 0.15 0.375]) % Panel 4
pcolor(X,Y,squeeze(nanmean(mean_phos.pres(244:335,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('SON','Fontsize',12)
caxis([0 .2])

axes('Position', [0.675 0.55 0.25 0.375], 'Visible', 'off');
h = colorbar('EastOutside','Fontsize',8);
caxis([0 .2])
set(h,'ytick',[0:.05:.2])
%set(h,'yticklabel',{'0','0.50','1.00','1.50','2.00'},'Fontsize',8)
ylabel(h, 'PO_4 (mmol/m^3)','Fontsize',8);


subplot('position',[0.1 0.075 0.15 0.375]) % Panel 5 
pcolor(X,Y,squeeze(nanmean(mean_phos.qpres([1:59 336:end],:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('DJF','Fontsize',12)
caxis([0 .2])

subplot('position',[0.3 0.075 0.15 0.375]) % Panel 6 
pcolor(X,Y,squeeze(nanmean(mean_phos.qpres(60:151,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('MAM','Fontsize',12)
caxis([0 .2])

subplot('position',[0.5 0.075 0.15 0.375]) % Panel 7
pcolor(X,Y,squeeze(nanmean(mean_phos.qpres(152:243,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('JJA','Fontsize',12)
caxis([0 .2])

subplot('position',[0.7 0.075 0.15 0.375]) % Panel 8
pcolor(X,Y,squeeze(nanmean(mean_phos.qpres(244:335,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42{\circ} N','','44{\circ} N','','46{\circ} N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88{\circ} W','','86{\circ} W',''},'Fontsize',8)
title('SON','Fontsize',12)
caxis([0 .2])

axes('Position', [0.675 0.075 0.25 0.375], 'Visible', 'off');
h = colorbar('EastOutside','Fontsize',8);
caxis([0 .2])
ylabel(h, 'PO_4 (mmol/m^3)','Fontsize',8);








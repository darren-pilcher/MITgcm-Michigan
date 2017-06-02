% Script used to plot difference plots of pCO2 between Past Nutrients and Quagga Present Nutrients simulations
% Created by Darren Pilcher 03/19/15 

close all 
clear 

% Past Nutrients simulation  
dir1 = '.../past/';

% Quagga Present Nutrients simulation
dir2 = '.../quagga_present/';

% Load Past Nutrients simulation pCO2 and CO2flux 
for yr = 2007:2010 
	ncload([dir1,num2str(yr),'/pco2.',num2str(yr),'.nc']);
	pco2_spatial.bar(yr-2006,:,:,:) = npzd_pCO2_ave;
	% Lakewide average
	lakewide_pco2.bar(yr-2006,:) = lake_surface_avg_michigan(npzd_pCO2_ave,1) * 10^6; clear npzd_pCO2_ave
end


% Load Quagga Present Nutrients simulation pCO2 and CO2flux
for yr = 2007:2010 
    	ncload([dir2,num2str(yr),'/pco2.',num2str(yr),'.nc']);
	pco2_spatial.new(yr-2006,:,:,:) = npzd_pCO2_ave;
	% Lakewide average
	lakewide_pco2.new(yr-2006,:) = lake_surface_avg_michigan(npzd_pCO2_ave,1) * 10^6; clear npzd_pCO2_ave
end

%%

% Take average over all 4 years and convert units to uatm 
pco2_spatial_avg.bar = squeeze(nanmean(pco2_spatial.bar,1)) * 10^6;
pco2_spatial_avg.new = squeeze(nanmean(pco2_spatial.new,1)) * 10^6; 

% Mask Land 
ncload('.../LMgrid1.kronos.nc','HFacC');
H = squeeze(HFacC(1,:,:)); clear HFacC
fi = find(H==0);

for t = 1:size(pco2_spatial_avg.bar,1)
	tmp = squeeze(pco2_spatial_avg.bar(t,:,:)); 
	tmp(fi) = nan; 
	pco2_spatial_avg.bar(t,:,:) = tmp; clear tmp

        tmp = squeeze(pco2_spatial_avg.new(t,:,:)); 
        tmp(fi) = nan; 
        pco2_spatial_avg.new(t,:,:) = tmp; clear tmp
end

% Plot Results 

V = [-200:50:200];
figure(3)
set(gcf, 'units','centimeters','position',[10 10 19 11.5]);
set(gcf,'Color',[1 1 1])

subplot('position',[0.1 0.55 0.15 0.375]) % [left bottom width height]
pcolor(X,Y,squeeze(nanmean(pco2_spatial_avg.new([1:59 336:end],:,:),1)) - squeeze(nanmean(pco2_spatial_avg.bar([1:59 336:end],:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
caxis([-200 200])
title('DJF','Fontsize',12)

subplot('position',[0.3 0.55 0.15 0.375])
pcolor(X,Y,squeeze(nanmean(pco2_spatial_avg.new(60:151,:,:),1)) - squeeze(nanmean(pco2_spatial_avg.bar(60:151,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
caxis([-200 200])
title('MAM','Fontsize',12)

subplot('position',[0.5 0.55 0.15 0.375])
pcolor(X,Y,squeeze(nanmean(pco2_spatial_avg.new(152:243,:,:),1)) - squeeze(nanmean(pco2_spatial_avg.bar(152:243,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
caxis([-200 200])
title('JJA','Fontsize',12)

subplot('position',[0.7 0.55 0.15 0.375])
pcolor(X,Y,squeeze(nanmean(pco2_spatial_avg.new(244:335,:,:),1)) - squeeze(nanmean(pco2_spatial_avg.bar(244:335,:,:),1)));
shading flat
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42^{\circ}N','','44^{\circ}N','','46^{\circ}N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88^{\circ}W','','86^{\circ}W',''},'Fontsize',8)
caxis([-200 200])
title('SON','Fontsize',12)

axes('Position', [0.675 0.55 0.25 0.375], 'Visible', 'off');
h = colorbar('EastOutside','Fontsize',8);
caxis([-200 200])
xlabel(h, 'pCO_2 ({\mu}atm)','Fontsize',8);
set(h,'Xtick',[-200:50:200])
hold off

subplot('position',[0.1 0.075 0.8 0.325]) % Panel 5
plot(nanmean(lakewide_pco2.bar,1),'k','Linewidth',3)
hold on
plot(nanmean(lakewide_pco2.new,1),'b','Linewidth',3)
plot(385*ones(1,365),'k--');
set(gca,'xtick',[1 32 60 91 121 152 182 213 244 274 305 336]);
set(gca,'xticklabel',{'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'},'Fontsize',8)
axis([0 365 200 600])
ylabel('pCO_2 ({\mu}atm)','Fontsize',8)
hleg = legend('Past Nutrients','Quagga Present Nutrients','Orientation','horizontal');
set(hleg,'Position', [0.1675 0.325 0.7 0.05])
set(hleg,'Fontsize',9);
title('Lakewide Averaged pCO_2 2007-2010','Fontsize',12)




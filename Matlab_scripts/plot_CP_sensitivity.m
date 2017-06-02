% Script used to plot NPP using different C:P ratios 
% Uses files created in the calc_lakewide_npp.m script

% Created by Darren Pilcher 05/2017

close all 
clear 

ncload('.../LMgrid1.kronos.nc','HFacC','X','Y','Z','Zl','Depth')

% Load depth integrated NPP from Past Nutrients simulation
load('.../lakewide.depth.integrated.npp1.past.mat','npp1')
load('.../lakewide.depth.integrated.npp2.past.mat','npp2')

npp1_CP200 = npp1 * (24*60*60)*(200)*(12.01); 
npp1_CP300 = npp1 * (24*60*60)*(300)*(12.01);
npp2_CP200 = npp2 * (24*60*60)*(200)*(12.01); 
npp2_CP150 = npp2 * (24*60*60)*(150)*(12.01);

npp_CP1.control = npp1_CP200 + npp2_CP200;
npp_CP2.control = npp1_CP200 + npp2_CP150;
npp_CP3.control = npp1_CP300 + npp2_CP150;

% Load depth integrated NPP from Quagga Present Nutrients simulation
load('.../lakewide.depth.integrated.npp1.qpres.mat','npp1')
load('.../lakewide.depth.integrated.npp2.qpres.mat','npp2')
npp1_CP200 = npp1 * (24*60*60)*(200)*(12.01);
npp1_CP300 = npp1 * (24*60*60)*(300)*(12.01);
npp2_CP200 = npp2 * (24*60*60)*(200)*(12.01);
npp2_CP150 = npp2 * (24*60*60)*(150)*(12.01);

npp_CP1.qpres = npp1_CP200 + npp2_CP200;
npp_CP2.qpres = npp1_CP200 + npp2_CP150;
npp_CP3.qpres = npp1_CP300 + npp2_CP150;

% Re-order dimensions for surface averaging function 
npp_CP1.control = permute(npp_CP1.control,[1 3 2 4]);
npp_CP2.control = permute(npp_CP2.control,[1 3 2 4]);
npp_CP3.control = permute(npp_CP3.control,[1 3 2 4]);

npp_CP1.qpres = permute(npp_CP1.qpres,[1 3 2 4]);
npp_CP2.qpres = permute(npp_CP2.qpres,[1 3 2 4]);
npp_CP3.qpres = permute(npp_CP3.qpres,[1 3 2 4]);

for yr = 1:4
	lakewide_npp_CP1.control(yr,:) = lake_surface_avg_michigan(squeeze(npp_CP1.control(yr,:,:,:)),2); 
	lakewide_npp_CP2.control(yr,:) = lake_surface_avg_michigan(squeeze(npp_CP2.control(yr,:,:,:)),2);
	lakewide_npp_CP3.control(yr,:) = lake_surface_avg_michigan(squeeze(npp_CP3.control(yr,:,:,:)),2);

        lakewide_npp_CP1.qpres(yr,:) = lake_surface_avg_michigan(squeeze(npp_CP1.qpres(yr,:,:,:)),2);
        lakewide_npp_CP2.qpres(yr,:) = lake_surface_avg_michigan(squeeze(npp_CP2.qpres(yr,:,:,:)),2);
        lakewide_npp_CP3.qpres(yr,:) = lake_surface_avg_michigan(squeeze(npp_CP3.qpres(yr,:,:,:)),2);
end


% Generate figure 
figure(1); 
set(gcf, 'units','centimeters','position',[10 10 19 14]);
set(gcf,'Color',[1 1 1])
subplot(2,1,1)
plot(mean(lakewide_npp_CP1.control,1),'k','Linewidth',2)
hold on 
plot(mean(lakewide_npp_CP2.control,1),'b','Linewidth',2)
plot(mean(lakewide_npp_CP3.control,1),'r','Linewidth',2)
title('Past Nutrients NPP','Fontsize',12)
ylabel('NPP (mgC/m^2/day)','Fontsize',10)
set(gca,'xtick',[1 32 60 91 121 152 182 213 244 274 305 336]);
set(gca,'xticklabel',{'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
axis([0 365 0 1500])
h = legend('C:P 200','C:P 200,150','C:P 300,150');
set(h,'FontSize',10);

subplot(2,1,2)
plot(mean(lakewide_npp_CP1.qpres,1),'k','Linewidth',2)
hold on 
plot(mean(lakewide_npp_CP2.qpres,1),'b','Linewidth',2)
plot(mean(lakewide_npp_CP3.qpres,1),'r','Linewidth',2)
title('Quagga Present Nutrients NPP','Fontsize',12)
ylabel('NPP (mgC/m^2/day)','Fontsize',10)
set(gca,'xtick',[1 32 60 91 121 152 182 213 244 274 305 336]);
set(gca,'xticklabel',{'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'})
axis([0 365 0 1500])
h = legend('C:P 200','C:P 200,150','C:P 300,150');
set(h,'FontSize',10);










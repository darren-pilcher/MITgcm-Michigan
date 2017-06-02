% Script used to compare model output phytoplankton biomass to EPA collected biovolume phyto data 

% Created by Darren Pilcher 09/2016

close all 
clear 

% Model Run Directory 
RUN_DIR = '.../quagga_present/';

% Load EPA data from Reavie et al. [2014] J. Great Lakes Res. 
[year month taxa station depth biovol] = textread('/home/disk/clark/dpilcher/michigan/obs/epa_phycount_michigan_0710.txt','%*n%n%*s%n%n%n%n%*n%*n%n%*n','headerlines',3,'delimiter','\t','emptyvalue',NaN);

% Sort EPA data into years, stations, and taxa. 
% Taxa are defined below as: 1 - cenric diatoms, 2 - pennate diatoms, 3 - chlorophytes (green algae), 4 - chrysophytes, 5 - cryptophytes, 6 - cyaanobacteria (blue-green algae)...
% 7 - pyrrophytes (dinoflagellates), 8 - unidentified entries
 
for yr = 2007:2010
	for s = 1:11
		for m = 1:2 
			fi_dia1 = find(year==yr&station==s&depth==1&month==m&taxa==1);
			fi_dia2 = find(year==yr&station==s&depth==1&month==m&taxa==2);
			fi_dino = find(year==yr&station==s&depth==1&month==m&taxa==7);
			fi_chl = find(year==yr&station==s&depth==1&month==m&taxa==3); 
			fi_chr = find(year==yr&station==s&depth==1&month==m&taxa==4);
			fi_cry = find(year==yr&station==s&depth==1&month==m&taxa==5);
			fi_cya = find(year==yr&station==s&depth==1&month==m&taxa==6);
			fi_total = find(year==yr&station==s&depth==1&month==m);
 
			epa_dia(yr-2006,s,m) = nansum(biovol(fi_dia1)) + nansum(biovol(fi_dia2)); 
			epa_dino(yr-2006,s,m) = nansum(biovol(fi_dino)) + nansum(biovol(fi_chl)) + nansum(biovol(fi_chr))+ nansum(biovol(fi_cry))+ nansum(biovol(fi_cya));
			epa_total(yr-2006,s,m) = nansum(biovol(fi_total)); 
		end
	end
end

% Replace any 0 values with nan due to no data 

for m = 1:2 
	tmp_dia = squeeze(epa_dia(:,:,m)); fi = find(tmp_dia==0); tmp_dia(fi)=nan; epa_dia(:,:,m) = tmp_dia; clear tmp_dia
	tmp_dino = squeeze(epa_dino(:,:,m)); fi = find(tmp_dino==0); tmp_dino(fi)=nan; epa_dino(:,:,m) = tmp_dino; clear tmp_dino
end

clear fi*

% Convert EPA observed values of biovolume (um^-3) to biomass (mgC/m^3) using a 0.20 conversion factor 
% from Rocha and Duncan 1985

epa_dia = epa_dia .* .20 ./ 1000; 
epa_dino = epa_dino .* .20 ./ 1000; 
epa_total = epa_total .* .20 ./ 1000; 

epa_dates(1,:,:) = [90 90 91 90 91 91 91 91 92 92 92; 214 213 213 214 214 214 215 215 215 215 216]';
epa_dates(2,:,:) = [106 105 106 106 106 108 108 108 109 109 109; 214 215 215 215 215 215 216 216 216 216 217]';
epa_dates(3,:,:) = [98 98 98 98 99 99 99 99 99 99 100; 217 216 216 217 217 217 218 218 218 218 219]';
epa_dates(4,:,:) = [92 91 92 92 92 92 93 92 93 93 93; 214 213 213 214 215 215 215 216 216 216 217]'; 


% Load model output.  Each file is a seperate EPA sampling location 

for yr = 2007:2010
	for s = 1:11
		ncload([RUN_DIR,num2str(yr),'/phy1.',num2str(yr),'.station',num2str(s),'.nc'])
		ncload([RUN_DIR,num2str(yr),'/phy2.',num2str(yr),'.station',num2str(s),'.nc'])
		
		model_phy1(yr-2006,s,1) = nanmean(phy1(epa_dates(yr-2006,s,1),1:4),2);
                model_phy1(yr-2006,s,2) = nanmean(phy1(epa_dates(yr-2006,s,2),1:4),2);
                model_phy2(yr-2006,s,1) = nanmean(phy2(epa_dates(yr-2006,s,1),1:4),2);
                model_phy2(yr-2006,s,2) = nanmean(phy2(epa_dates(yr-2006,s,2),1:4),2);

		clear phy1 phy2
	end
end

% Convert Model output from units of mmolP/m^3 to mgC/m^3 

model_phy1 = model_phy1 .* 200 .* 12; 
model_phy2 = model_phy2 .* 200 .* 12; 

model_total = model_phy1 + model_phy2; 

plot_dino1 = [squeeze(nanmean(epa_dino(:,:,1),2)) squeeze(nanmean(model_phy1(:,:,1),2))];
plot_dino2 = [squeeze(nanmean(epa_dino(:,:,2),2)) squeeze(nanmean(model_phy1(:,:,2),2))]; 

plot_dino1_err = [squeeze(nanstd(epa_dino(:,:,1),0,2)) squeeze(nanstd(model_phy1(:,:,1),0,2))]; 
plot_dino2_err = [squeeze(nanstd(epa_dino(:,:,2),0,2)) squeeze(nanstd(model_phy1(:,:,2),0,2))];

plot_dia1 = [squeeze(nanmean(epa_dia(:,:,1),2)) squeeze(nanmean(model_phy2(:,:,1),2))];
plot_dia2 = [squeeze(nanmean(epa_dia(:,:,2),2)) squeeze(nanmean(model_phy2(:,:,2),2))];

plot_dia1_err = [squeeze(nanstd(epa_dia(:,:,1),0,2)) squeeze(nanstd(model_phy2(:,:,1),0,2))]; 
plot_dia2_err = [squeeze(nanstd(epa_dia(:,:,2),0,2)) squeeze(nanstd(model_phy2(:,:,2),0,2))]; 

plot_total1 = [squeeze(nanmean(epa_total(:,:,1),2)) squeeze(nanmean(model_total(:,:,1),2))]; 
plot_total2 = [squeeze(nanmean(epa_total(:,:,2),2)) squeeze(nanmean(model_total(:,:,2),2))];

plot_total1_err = [squeeze(nanstd(epa_total(:,:,1),0,2)) squeeze(nanstd(model_total(:,:,1),0,2))];
plot_total2_err = [squeeze(nanstd(epa_total(:,:,2),0,2)) squeeze(nanstd(model_total(:,:,2),0,2))];

plot_time_avg_total1 = [squeeze(nanmean(epa_total(:,:,1),1))' squeeze(nanmean(model_total(:,:,1),1))']; 
plot_time_avg_total2 = [squeeze(nanmean(epa_total(:,:,2),1))' squeeze(nanmean(model_total(:,:,2),1))'];

plot_time_avg_total1_err = [squeeze(nanstd(epa_total(:,:,1),0,1))' squeeze(nanstd(model_total(:,:,1),0,1))'];
plot_time_avg_total2_err = [squeeze(nanstd(epa_total(:,:,2),0,1))' squeeze(nanstd(model_total(:,:,2),0,1))'];

% Compute RMSE values for spring and summer total phytoplankton 
% RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

RMSE_total1 = sqrt(mean((plot_total1(:,2) - plot_total1(:,1)).^2)); 
RMSE_total1 = sqrt(mean((plot_total2(:,2) - plot_total2(:,1)).^2));  
 
% Plot results 
figure(1)
subplot(2,1,1) 
[h, herr] = barwitherr(plot_dino1_err,[2007:1:2010],plot_dino1);
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
ylim([0 90])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
title('Other PHY Spring','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

subplot(2,1,2)
[h, herr] = barwitherr(plot_dino2_err,[2007:1:2010],plot_dino2);
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
ylim([0 90])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
title('Other PHY Summer','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

figure(2)
subplot(2,1,1)
[h, herr] = barwitherr(plot_dia1_err,[2007:1:2010],plot_dia1);
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
ylim([0 90])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
title('Diatom Spring','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

subplot(2,1,2)
[h, herr] = barwitherr(plot_dia2_err,[2007:1:2010],plot_dia2);
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
ylim([0 90])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
title('Diatom Summer','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

% For total PHY
figure(3) 
subplot(2,1,1)
[h,herr] = barwitherr(plot_total1_err,[2007:1:2010],plot_total1);
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
ylim([0 110])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
title('Spring PHY','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

subplot(2,1,2)
[h2,herr2] = barwitherr(plot_total2_err,[2007:1:2010],plot_total2);
set(h2(1),'FaceColor',[0.1 0.1 0.1]);
set(h2(2),'FaceColor',[0.8 0.8 0.8]);
set(herr2(:), 'LineWidth', 3,'Color',[0.4 0.4 0.4])
ylim([0 110])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
title('Summer PHY','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

figure(4)
set(gcf, 'units','centimeters','position',[10 10 19 13]); %17
subplot(2,1,1)
[h,herr] = barwitherr(plot_time_avg_total1_err,[1:11],plot_time_avg_total1);
set(h(1),'FaceColor',[0.1 0.1 0.1]);
set(h(2),'FaceColor',[0.8 0.8 0.8]);
set(herr(:), 'LineWidth', 1,'Color',[0.4 0.4 0.4])
ylim([0 110])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
xlabel('Station Number','Fontsize',10)
set(gca,'xticklabel',{'MI 11' 'MI 17' 'MI 18M' 'MI 19' 'MI 23' 'MI 27' 'MI 32' 'MI 34' 'MI 40' 'MI 41' 'MI 47'},'Fontsize',8)
title('Spring PHY','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')

subplot(2,1,2)
[h2,herr2] = barwitherr(plot_time_avg_total2_err,[1:11],plot_time_avg_total2);
set(h2(1),'FaceColor',[0.1 0.1 0.1]);
set(h2(2),'FaceColor',[0.8 0.8 0.8]);
set(herr2(:), 'LineWidth', 1,'Color',[0.4 0.4 0.4])
ylim([0 110])
ylabel('Biomass (mgC/m^3)','Fontsize',10)
xlabel('Station Number','Fontsize',10)
set(gca,'xticklabel',{'MI 11' 'MI 17' 'MI 18M' 'MI 19' 'MI 23' 'MI 27' 'MI 32' 'MI 34' 'MI 40' 'MI 41' 'MI 47'},'Fontsize',8)
title('Summer PHY','Fontsize',12)
legend('OBS','MODEL','Location','Northwest')



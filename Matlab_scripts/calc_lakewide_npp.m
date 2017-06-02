% Script used to calculate lakewide depth integrated NPP files that are used in the ...
% determine_michigan_gpp_multiple_years.m and plot_CP_sensitivity scripts  

% Created by Darren Pilcher 03/18/15

clear 
% Change directory for each simulation 
dir = '.../past/';

for yr = 2007:2010
    ncload([dir,num2str(yr),'/phygrow1.',num2str(yr),'.nc']);
    model_phygrow1(yr-2006,:,:,:,:) = squeeze(npzd_phygrow1_ave); clear npzd_phygrow1_ave
    ncload([dir,num2str(yr),'/phygrow2.',num2str(yr),'.nc']);
    model_phygrow2(yr-2006,:,:,:,:) = squeeze(npzd_phygrow2_ave); clear npzd_phygrow2_ave
end

load('.../volume_michigan.mat');
ncload('.../LMgrid1.kronos.nc','rA','HFacC')

NPP1 = model_phygrow1;
NPP2 = model_phygrow2;
clear model_phygrow*

for d = 1:29
	thick(:,:,d) = (squeeze(volume(:,:,d)) ./ rA') *-1;
end

thick = permute(thick,[3 2 1]);

for yr = 1:4
        for t = 1:365
                tmp1(yr,:,:,:,t) = squeeze(NPP1(yr,t,:,:,:)) .* thick;
                tmp2(yr,:,:,:,t) = squeeze(NPP2(yr,t,:,:,:)) .* thick;
        end
end

clear NPP1
clear NPP2

% Sum over all depths 

npp1 = squeeze(nansum(tmp1,2));
npp2 = squeeze(nansum(tmp2,2));

clear tmp1 tmp2

save('.../lakewide.depth.integrated.npp1.past.mat','npp1')
save('.../lakewide.depth.integrated.npp2.past.mat','npp2')






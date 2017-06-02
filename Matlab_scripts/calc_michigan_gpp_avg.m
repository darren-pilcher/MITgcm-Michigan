% Script used to calculate a 2007-2010 average GPP based on 2 Fahnenstiel
% station locations 
% Also calculates 4 year, depth integrated average for Phy1, Phy2, and Zoo
% Saves output as a .mat files that are used in determine_michigan_multiple_years_quagga.m  
% Created by Darren Pilcher 03/17/15

clear
% Change directory for each simulation
dir = ('.../past/');

for yr = 2007:2010 
    ncload([dir,num2str(yr),'/phygrow1.92x.99y.',num2str(yr),'.nc']);
    model_phygrow1(yr-2006,1,:,:) = squeeze(npzd_phygrow1_ave); clear npzd_phygrow1_ave
    ncload([dir,num2str(yr),'/phygrow2.92x.99y.',num2str(yr),'.nc']);
    model_phygrow2(yr-2006,1,:,:) = squeeze(npzd_phygrow2_ave); clear npzd_phygrow2_ave
    ncload([dir,num2str(yr),'/phy1.92x.99y.',num2str(yr),'.nc']);
    model_phy1(yr-2006,1,:,:) = squeeze(phy1); clear phy1
    ncload([dir,num2str(yr),'/phy2.92x.99y.',num2str(yr),'.nc']);
    model_phy2(yr-2006,1,:,:) = squeeze(phy2); clear phy2
    ncload([dir,num2str(yr),'/zoo.92x.99y.',num2str(yr),'.nc']);
    model_zoo(yr-2006,1,:,:) = squeeze(zoo); clear zoo
end

for yr = 2007:2010 
    ncload([dir,num2str(yr),'/phygrow1.88x.88y.',num2str(yr),'.nc']);
    model_phygrow1(yr-2006,2,:,:) = squeeze(npzd_phygrow1_ave); clear npzd_phygrow1_ave
    ncload([dir,num2str(yr),'/phygrow2.88x.88y.',num2str(yr),'.nc']);
    model_phygrow2(yr-2006,2,:,:) = squeeze(npzd_phygrow2_ave); clear npzd_phygrow2_ave
    ncload([dir,num2str(yr),'/phy1.88x.88y.',num2str(yr),'.nc']);
    model_phy1(yr-2006,2,:,:) = squeeze(phy1); clear phy1
    ncload([dir,num2str(yr),'/phy2.88x.88y.',num2str(yr),'.nc']);
    model_phy2(yr-2006,2,:,:) = squeeze(phy2); clear phy2
    ncload([dir,num2str(yr),'/zoo.88x.88y.',num2str(yr),'.nc']);
    model_zoo(yr-2006,2,:,:) = squeeze(zoo); clear zoo
end        

ncload('.../LMgrid1.kronos.nc','rA');
load('.../volume_michigan.mat');


thick1 = squeeze(volume(92,99,:)/rA(92,99))*-1; % -1 for sign change
thick2 = squeeze(volume(88,88,:)/rA(88,88))*-1; % -1 for sign change

nt = size(model_phygrow1,4);
nyr = size(model_phygrow1,1); 

for t = 1:nt 
    for yr = 1:nyr 
        phygrow1(yr,1,:,t) = squeeze(model_phygrow1(yr,1,:,t)) .*thick1;
        phygrow1(yr,2,:,t) = squeeze(model_phygrow1(yr,2,:,t)) .*thick2;

        phygrow2(yr,1,:,t) = squeeze(model_phygrow2(yr,1,:,t)) .*thick1;
        phygrow2(yr,2,:,t) = squeeze(model_phygrow2(yr,2,:,t)) .*thick2;
        
        phy1(yr,1,:,t) = squeeze(model_phy1(yr,1,:,t)) .*thick1;
        phy1(yr,2,:,t) = squeeze(model_phy1(yr,2,:,t)) .*thick2;
        
        phy2(yr,1,:,t) = squeeze(model_phy2(yr,1,:,t)) .*thick1;
        phy2(yr,2,:,t) = squeeze(model_phy2(yr,2,:,t)) .*thick2;
        
        zoo(yr,1,:,t) = squeeze(model_zoo(yr,1,:,t)) .*thick1;
        zoo(yr,2,:,t) = squeeze(model_zoo(yr,2,:,t)) .*thick2;

    end
end

clear model* 

save([dir,'phygrow1_depth_int_FA_0710.mat'],'phygrow1')
save([dir,'phygrow2_depth_int_FA_0710.mat'],'phygrow2')
save([dir,'phy1_depth_int_FA_0710.mat'],'phy1')
save([dir,'phy2_depth_int_FA_0710.mat'],'phy2')
save([dir,'zoo_depth_int_FA_0710.mat'],'zoo')








% Script used to generate a .mat file that contains all 4 years of a single variable, in space, but in units of /m^2
% Files are used in plot_depth_int_npp_zoo.m script 

% Created by Darren Pilcher 03/21/16

clear 

% Select which variable to create file for
var = 'phy';
%var = 'zoo';
%var = 'temp';

% Change through simulations
dir = '.../quagga_jan28/';

load('.../volume_michigan.mat');
ncload('.../LMgrid1.kronos.nc','rA','HFacC','Depth')


switch lower(var)
	case 'phy'

for yr = 2007:2010
    ncload([dir,num2str(yr),'/phy1.',num2str(yr),'.nc']);
    model_phy1(yr-2006,:,:,:,:) = squeeze(phy1); clear phy1
    ncload([dir,num2str(yr),'/phy2.',num2str(yr),'.nc']);
    model_phy2(yr-2006,:,:,:,:) = squeeze(phy2); clear phy2
end 

% Since we're just looking at the top 25m, squeeze all variables to just this depth to save computational time 
model_phy1 = squeeze(model_phy1(:,:,1:5,:,:)); 
model_phy2 = squeeze(model_phy2(:,:,1:5,:,:)); 

volume = squeeze(volume(:,:,1:5)); 

for d = 1:5
thick(:,:,d) = (squeeze(volume(:,:,d)) ./ rA') *-1;
end

thick = permute(thick,[3 2 1]);

for yr = 1:4
        for t = 1:365
                tmp1(yr,:,:,:,t) = squeeze(model_phy1(yr,t,:,:,:)) .* thick;
		tmp2(yr,:,:,:,t) = squeeze(model_phy2(yr,t,:,:,:)) .* thick;
        end
end

clear model_phy*

% Need to mask out land so that averaging over 25m does not include below lake bottom points that are outputed with a value of 0 
for d = 1:5
	H = squeeze(HFacC(d,:,:));
	fi = find(H==0); 
	for t = 1:365
		for yr = 1:4 
			temp = squeeze(tmp1(yr,d,:,:,t)); 
			temp(fi) = nan; 
			tmp1(yr,d,:,:,t) = temp; clear temp

			temp = squeeze(tmp2(yr,d,:,:,t)); 
                        temp(fi) = nan; 
                        tmp2(yr,d,:,:,t) = temp; clear temp	
		end
	end
end
			

% Find 25m average value
phy1 = squeeze(nanmean(tmp1,2));
phy2 = squeeze(nanmean(tmp2,2));

% Change name based on simulation
save('.../lakewide.25m.avg.phy1.qpres.mat','phy1')
save('.../lakewide.25m.avg.phy2.qpres.mat','phy2')

	case 'zoo'

for yr = 2007:2010
    ncload([dir,num2str(yr),'/zoo.',num2str(yr),'.nc']);
    model_zoo(yr-2006,:,:,:,:) = squeeze(zoo); clear zoo
end

% Since we're just looking at the top 25m, squeeze all variables to just this depth to save computational time
model_zoo = squeeze(model_zoo(:,:,1:5,:,:));
volume = squeeze(volume(:,:,1:5));

for d = 1:5
thick(:,:,d) = (squeeze(volume(:,:,d)) ./ rA') *-1;
end

thick = permute(thick,[3 2 1]);

for yr = 1:4
        for t = 1:365
                tmp(yr,:,:,:,t) = squeeze(model_zoo(yr,t,:,:,:)) .* thick;
        end
end

clear model_zoo

% Need to mask out land so that averaging over 25m does not include below lake bottom points that are outputed with a value of 0
for d = 1:5
        H = squeeze(HFacC(d,:,:));
        fi = find(H==0);
        for t = 1:365
                for yr = 1:4
                        temp = squeeze(tmp(yr,d,:,:,t));
                        temp(fi) = nan;
                        tmp(yr,d,:,:,t) = temp; clear temp

                end
        end
end


% Find 25m average value
zoo = squeeze(nanmean(tmp,2));

% Change name based on simulation
save('.../lakewide.25m.avg.zoo.qpres.mat','zoo')


	case 'temp'


for yr = 2007:2010
    ncload([dir,num2str(yr),'/Ttave.',num2str(yr),'.nc']);
    model_temp(yr-2006,:,:,:,:) = Ttave; clear Ttave
end

% Calculate 4-year average values 
model_temp = squeeze(nanmean(model_temp,1)); 

% Calculate lakewide averaged surface temperature 
[lakewide_surf_temp] = lake_surface_avg_michigan(squeeze(model_temp(:,1,:,:)),1); 

nx = size(model_temp,3); 
ny = size(model_temp,4); 

model_bottom_temp = nan(365,nx,ny); 

% Find the lake bottom point and place into a new matrix 

for y = 1:ny
	for x = 1:nx 
		fi = find(squeeze(HFacC(:,x,y))>0);
		if isempty(fi) == 0  
                	model_bottom_temp(:,x,y) = model_temp(:,fi(end),x,y); 
		else
		end
	end 
end

[lakewide_bottom_temp] = lake_surface_avg_michigan(model_bottom_temp,1); 

save(.../lakewide_surf_temp_0710.mat','lakewide_surf_temp');
save(.../lakewide_bottom_temp_0710.mat'],'lakewide_bottom_temp');


end







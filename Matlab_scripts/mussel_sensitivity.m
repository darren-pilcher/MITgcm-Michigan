% Script used to teast out the model sensitivity to the mussel population numbers, 
% as compared to the variability in the Past Nutrients simulation 

% Created by Darren Pilcher 01/2017

close all 
clear 

load('.../lakewide.depth.integrated.npp.past.mat')
NPP.past = npp; clear npp
load('.../lakewide.depth.integrated.npp.qpast.mat')
NPP.qpast = npp; clear npp

for t = 1:4                                                              
	lakewide_avg_NPP.past(t,:) = lake_surface_avg_michigan(squeeze(NPP.past(t,:,:,:)),2);
	lakewide_avg_NPP.qpast(t,:) = lake_surface_avg_michigan(squeeze(NPP.qpast(t,:,:,:)),2);
end

ncload('.../LMgrid1.kronos.nc','HFacC','X','Y','Z','Zl','Depth','drF')
H = squeeze(HFacC(1,:,:));

mask1 = find(Depth' > 30); 
mask2 = find(Depth' < 30 | Depth' > 50);
mask3 = find(Depth' < 50 | Depth' > 90);
mask4 = find(Depth' < 90);

% Separate lake NPP into the different mussel depth bins
for t = 1:4 
	for d = 1:365
		tmp = squeeze(NPP.past(t,:,:,d)); 
		tmp(mask1) = nan; 
		NPP_depth1.past(t,:,:,d) = tmp; 

                tmp = squeeze(NPP.qpast(t,:,:,d));
                tmp(mask1) = nan;
                NPP_depth1.qpast(t,:,:,d) = tmp;  		

	        tmp = squeeze(NPP.past(t,:,:,d));
                tmp(mask2) = nan;
                NPP_depth2.past(t,:,:,d) = tmp;

                tmp = squeeze(NPP.qpast(t,:,:,d));
                tmp(mask2) = nan;
                NPP_depth2.qpast(t,:,:,d) = tmp;
	
                tmp = squeeze(NPP.past(t,:,:,d));
                tmp(mask3) = nan;
                NPP_depth3.past(t,:,:,d) = tmp;

                tmp = squeeze(NPP.qpast(t,:,:,d));
                tmp(mask3) = nan;
                NPP_depth3.qpast(t,:,:,d) = tmp;

                tmp = squeeze(NPP.past(t,:,:,d));
                tmp(mask4) = nan;
                NPP_depth4.past(t,:,:,d) = tmp;

                tmp = squeeze(NPP.qpast(t,:,:,d));
                tmp(mask4) = nan;
                NPP_depth4.qpast(t,:,:,d) = tmp;

	end
end	

% Take lakewide average for each depth bin 
for t = 1:4
	lakewide_NPP_depth1.past(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth1.past(t,:,:,:)),2); 
	lakewide_NPP_depth1.qpast(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth1.qpast(t,:,:,:)),2);	

	lakewide_NPP_depth2.past(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth2.past(t,:,:,:)),2);
        lakewide_NPP_depth2.qpast(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth2.qpast(t,:,:,:)),2);

        lakewide_NPP_depth3.past(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth3.past(t,:,:,:)),2);
        lakewide_NPP_depth3.qpast(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth3.qpast(t,:,:,:)),2);

        lakewide_NPP_depth4.past(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth4.past(t,:,:,:)),2);
        lakewide_NPP_depth4.qpast(t,:) = lake_surface_avg_michigan(squeeze(NPP_depth4.qpast(t,:,:,:)),2);
end	

% Take difference
lakewide_NPP_diff_depth1 = lakewide_NPP_depth1.qpast - lakewide_NPP_depth1.past; 
lakewide_NPP_diff_depth2 = lakewide_NPP_depth2.qpast - lakewide_NPP_depth2.past; 
lakewide_NPP_diff_depth3 = lakewide_NPP_depth3.qpast - lakewide_NPP_depth3.past; 
lakewide_NPP_diff_depth4 = lakewide_NPP_depth4.qpast - lakewide_NPP_depth4.past; 


% Combine into bar plot variable
bar_diff_depth4_mix = [(mean(lakewide_NPP_diff_depth1(1,:)) / mean(lakewide_NPP_depth1.past(1,:))) ...
			(mean(lakewide_NPP_diff_depth1(2,:)) / mean(lakewide_NPP_depth1.past(2,:))) ...
			(mean(lakewide_NPP_diff_depth1(3,:)) / mean(lakewide_NPP_depth1.past(3,:))) ...
			(mean(lakewide_NPP_diff_depth1(4,:)) / mean(lakewide_NPP_depth1.past(4,:))); ...
			(mean(lakewide_NPP_diff_depth2(1,:)) / mean(lakewide_NPP_depth2.past(1,:))) ...
                        (mean(lakewide_NPP_diff_depth2(2,:)) / mean(lakewide_NPP_depth2.past(2,:))) ...
                        (mean(lakewide_NPP_diff_depth2(3,:)) / mean(lakewide_NPP_depth2.past(3,:))) ...
                        (mean(lakewide_NPP_diff_depth2(4,:)) / mean(lakewide_NPP_depth2.past(4,:))); ...
			(mean(lakewide_NPP_diff_depth3(1,:)) / mean(lakewide_NPP_depth3.past(1,:))) ...
                        (mean(lakewide_NPP_diff_depth3(2,:)) / mean(lakewide_NPP_depth3.past(2,:))) ...
                        (mean(lakewide_NPP_diff_depth3(3,:)) / mean(lakewide_NPP_depth3.past(3,:))) ...
                        (mean(lakewide_NPP_diff_depth3(4,:)) / mean(lakewide_NPP_depth3.past(4,:))); ...
			(mean(lakewide_NPP_diff_depth4(1,:)) / mean(lakewide_NPP_depth4.past(1,:))) ...
                        (mean(lakewide_NPP_diff_depth4(2,:)) / mean(lakewide_NPP_depth4.past(2,:))) ...
                        (mean(lakewide_NPP_diff_depth4(3,:)) / mean(lakewide_NPP_depth4.past(3,:))) ...
                        (mean(lakewide_NPP_diff_depth4(4,:)) / mean(lakewide_NPP_depth4.past(4,:)))];




% From mussel depth placement code 

% For timing of sampling, Nalepa (2010) says that the sample dates for
% 2006-2008 are late summer/fall, so use September 1st as the date

% This vector runs from September 1st, 2006 - September 1st, 2008
% Note that leap day is accounted for here
dates = [1 366 732 1462];
new_dates = [1:1462];

% This is data estimated from Nalepa et al. (2010).  2010 Data is from the NOAA tecnical memorandum, i.e. Nalepa et al. [2014]
% The first row corresponds to 2006,2007,2008, and 2010 mussel population density at <30 m
% Second row is 31-50m, third row is 51-90m and last row is >90m

data = [12500 5000 20000 8679; 13000 13000 12000 8867; 4000 11000 15000 15122; 400 800 2000 4223]; % This is including ONLY SOUTHERN BASIN station locations in the Nalepa et al. [2014] technical memorandum.  These are the same stations used in Nalepa et al. [2010]

for d = 1:4
    interp_data(d,:) = interp1(dates,data(d,:),new_dates);
end

mussel = zeros(1462,276,200);

for x = 1:200
    for y = 1:276
        if Depth(y,x) >0 && Depth(y,x) < 30
            mussel(:,y,x) = interp_data(1,:);
        elseif Depth(y,x) >= 30 && Depth(y,x) < 50
            mussel(:,y,x) = interp_data(2,:);
        elseif Depth(y,x) >=50 && Depth(y,x) < 90
            mussel(:,y,x) = interp_data(3,:);
        elseif Depth(y,x) >=90
            mussel(:,y,x) = interp_data(4,:);
        end
    end
end

mussel_thick = zeros(1462,276,200);
for x = 1:200
    for y = 1:276
        if H(y,x)>0
            fi = find(HFacC(:,y,x)==0);
            mussel_thick(:,y,x) = mussel(:,y,x) ./ drF(fi(1)-1);
        else
        end
    end
end


% This filtration rate is equal to 5.4 L/day
% This is in units of /second
mussel_units = mussel_thick .* 225 ./ 3600 ./ 1000 ./ 1000;  % 225 ml/hr/mussel / 3600 s / 1000 ml  / 1000 L

% So, now we have a file that is September 1st, 2006 - September 1st, 2008
% But, we want to make a forcing file that's January 1st 2007 - January 1st
% 2011 (need extra date for MITgcm to finish 2010)

% Thus, cut out the 2006 part the we don't need and repeat September 1st
% 2008 for rest of timeframe since we don't have any mussel data for after
% 2008

% September 1st is Julian day 244 in non leap year, so we want to start
% with point 123 as this is January 1st 2007

mussel_tmp = mussel_units(123:end,:,:);

% Now add last day
end_point = mussel_units(1462,:,:);

last_part = repmat(end_point,[122 1 1]);  % Should be 1462 days since this accounts for the 1 extra day (Jan 1, 2011) that is needded for the model to complete Dec 31 2010

mussel_new = cat(1,mussel_tmp,last_part);

% Make a variable to plot for thesis that shows the average mussel population for each model year

plot_mussel = cat(1,mussel(123:end,:,:),repmat(mussel(end,:,:),[122 1 1]));

% Mask land as NaN for plotting
H = squeeze(HFacC(1,:,:));
fih = find(H==0);

for t = 1:size(plot_mussel,1)
    tmp = squeeze(plot_mussel(t,:,:));
    tmp(fih) = nan;
    plot_mussel(t,:,:) = tmp;
end

% Plot figure 
figure(1)
bar(bar_diff_depth4_mix)
ylabel('% Change','Fontsize',10)
axis([0.5 4.5 -.7 0])
set(gca,'yticklabel',{'-70' '-60' '-50' '-40' '-30' '-20' '-10' '0'})
set(gca,'xticklabel',{'{\leq} 30m' '31-50m' '51-90m' '> 90m'})
title('NPP Decrease from Mussels','Fontsize',12)
h = legend('2007','2008','2009','2010','Location','SouthEast');
set(h,'Fontsize',10)




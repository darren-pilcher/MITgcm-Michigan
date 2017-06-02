%Created by James Kralj April 2013
% Updated by Darren Pilcher Feb 2015 

% This forcing file will be made in terms of the number of times that the
% grid cell is completely filtered through by mussels per second.  The per
% second unit is neccessary so that it can be directly inputed into the
% model.

close all
clear

ncload('.../LMgrid1.kronos.nc','X','Y','Z','HFacC','Depth','drF'); 
H = squeeze(HFacC(1,:,:));

% For timing of sampling, Nalepa (2010) says that the sample dates for
% 2006-2008 are late summer/fall, so use September 1st as
% the date 

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

% Save File for using in MITgcm.michigan
write_binary('/Volumes/Darren2/Michigan/forcing/','quagga.5.4L.south.stations.2007-2010.bin',mussel_new)

% Make a figure showing yearly mussel population density
plot_mussel = cat(1,mussel(123:end,:,:),repmat(mussel(end,:,:),[122 1 1])); 

% Mask land as NaN for plotting 
H = squeeze(HFacC(1,:,:));
fih = find(H==0);

for t = 1:size(plot_mussel,1)
    tmp = squeeze(plot_mussel(t,:,:));
    tmp(fih) = nan;
    plot_mussel(t,:,:) = tmp;
end 

%%


V = [0:2000:18000]; 
% Generate Figure 
figure(1);  
set(gcf, 'units','centimeters','position',[10 10 19 11.5]);
set(gcf,'Color',[1 1 1])
subplot('position',[0.1 0.55 0.15 0.375]) % [left bottom width height]
[cs,h] = contourf(X,Y,squeeze(nanmean(plot_mussel(1:365,:,:),1)),V);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42° N','','44° N','','46° N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88° W','','86° W',''},'Fontsize',8)
title('2007','Fontsize',12)
caxis([0 18000])
 

subplot('position',[0.3 0.55 0.15 0.375])
[cs,h] = contourf(X,Y,squeeze(nanmean(plot_mussel(366:731,:,:),1)),V);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42° N','','44° N','','46° N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88° W','','86° W',''},'Fontsize',8)
title('2008','Fontsize',12)
caxis([0 18000])

subplot('position',[0.5 0.55 0.15 0.375])
[cs,h] = contourf(X,Y,squeeze(nanmean(plot_mussel(732:1096,:,:),1)),V);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42° N','','44° N','','46° N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88° W','','86° W',''},'Fontsize',8)
title('2009','Fontsize',12)
caxis([0 18000])

subplot('position',[0.7 0.55 0.15 0.375])
[cs,h] = contourf(X,Y,squeeze(nanmean(plot_mussel(1097:1461,:,:),1)),V);
set(gca,'ytick',[42:1:46]);
set(gca,'yticklabel',{'42° N','','44° N','','46° N'},'Fontsize',8)
set(gca,'xtick',[-88:-85])
set(gca,'xticklabel',{'88° W','','86° W',''},'Fontsize',8)
title('2010','Fontsize',12)
caxis([0 18000])

axes('Position', [0.1 0.425 0.75 0.1], 'Visible', 'off');
h = colorbar('SouthOutside','Fontsize',12);
caxis([0 18000])
set(h,'Xtick',[0:2000:18000])
xlabel(h, 'Mussel Density (# mussels/m^2)','Fontsize',12);







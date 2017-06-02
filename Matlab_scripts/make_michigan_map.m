% script used to create a map of Lake Michigan that includes contoured
% bathymetry, EPA station locations, and NDBC buoys. 

% created by Darren Pilcher 08/22/2103

close all
clear

% Load bathymetry data 
ncload('.../LMgrid1.kronos.nc','HFacC','Depth','X','Y')

% Mask land as NaN

H_land = squeeze(HFacC(1,:,:));
fi_land = find(H_land==0);
Depth(fi_land)=nan;

% Second mask to mask water as NaN
H_water = squeeze(HFacC(1,:,:));
fi_water = find(H_water==1);
H_water(fi_water)=nan;

% EPA station locations and model equivalent locations: 

% Station 1  (MI 11)   42.38 N  -87.00 E  (67,50)
% Station 2  (MI 17)   42.73 N  -87.42 E  (41,71) 
% Station 3  (MI 18M)  42.73 N  -87.00 E  (67,71) 
% Station 4  (MI 19)   42.73 N  -86.58 E  (92,71)
% Station 5  (MI 23)   43.13 N  -87.00 E  (67,95)
% Station 6  (MI 27M)  43.60 N  -86.92 E  (71,123)
% Station 7  (MI 32)   44.14 N  -87.23 E  (53,155)
% Station 8  (MI 34)   44.09 N  -86.77 E  (80,152)
% Station 9  (MI 40)   44.76 N  -86.97 E  (68,192)
% Station 10 (MI 41M)  44.74 N  -86.72 E  (83,191)
% Station 11 (MI 47)   45.18 N  -86.37 E  (104,218)

% NDBC Buoy Locations and model equivalent locations: 

% South basin buoy    42.674N 87.026W   (65,67)
% North basin buoy    45.344N 86.411W   (102,228)

% Create matrix of NaNs to then fill in with numbers
stations = nan(200,276);

% Mark grid cell points with locations.  These are just placeholders for using a later graphical program to make larger points (e.g. Powerpoint or Illustrator).
% EPA stations 
stations(67,50) = 50; stations(66,50) = 1; stations(68,50) = 1; stations(67,49) = 1; stations(67,51) = 1; to make a cross    
% stations(41,71) = 50;
% stations(67,71) = 50;
% stations(92,71) = 50;
% stations(67,95) = 50;
% stations(71,123) = 50;
% stations(53,155) = 50;
% stations(80,152) = 50;
% stations(68,192) = 50;
% stations(83,191) = 50;
% stations(104,218) = 50;
% NDBC stations
stations(65,67) = 100;
stations(102,228) = 100;
% Fahnenstiel Locations 
stations(88,88) = 200;
stations(92,99) = 200;
% Lake Express Ferry Transect:
stations(15,87) = 300;
stations(105,100) = 300;

% Contour interval for bathymetry 
V = [0:50:300];

% Put everything together into 1 Figure 
figure(1); cla;
set(gcf, 'units','centimeters','position',[10 10 9.5 11.5]);
set(gcf,'Color',[1 1 1])
h = pcolor(X,Y,H_water);
shading flat
hold on 
contour(X,Y,Depth,V,'LineColor','k')
pcolor(X,Y,stations')
shading flat
title('Lake Michigan Observation Locations','Fontsize',12)
set(gca,'xtick',[-88:1:-85]);
set(gca,'xticklabel',{'88° W','87° W','86° W','85° W'},'Fontsize',8)
set(gca,'ytick',[42:1:46])
set(gca,'yticklabel',{'42° N','43° N','44° N','45° N','46° N'},'Fontsize',8)
colormap([0.6 0.6 0.6])






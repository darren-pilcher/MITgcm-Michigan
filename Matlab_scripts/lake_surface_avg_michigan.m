% find lake wide average of surface variable
% written by val bennington 
% does NOT work when you input a 2D or 4D variable
% updated by Darren Pilcher 8/4/11
% changed so that sst(x,y,t)
% Modified for a netcdf and .mat file version by DJP 01/2017

% type 1 = netcdf, type 2 = .mat

function [lakewide] = lake_surface_avg_michigan(sst,type)


if type == 1
	ncload('.../LMgrid1.kronos.nc','HFacC');
	H = squeeze(HFacC(1,:,:)); clear HFacC
	fi = find(H==0);

	nt = size(sst,1);
	ny = size(sst,2);
	nx = size(sst,3);

	for t = 1:nt
		day = squeeze(sst(t,:,:));
		day(fi)=nan;
		tmp(t,:) = reshape(day,1,nx*ny);
	end

elseif type == 2
	ncload('.../LMgrid1.kronos.nc','HFacC');
	H = squeeze(HFacC(1,:,:))'; clear HFacC
	fi = find(H==0);

        nt = size(sst,3);
        ny = size(sst,2);
        nx = size(sst,1);

        for t = 1:nt
                day = squeeze(sst(:,:,t));
                day(fi)=nan;
                tmp(t,:) = reshape(day,1,nx*ny);
        end

end 

lakewide = nanmean(tmp,2);

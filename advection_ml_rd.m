% Compute the magnitude of advection given a velocity and tracer field,
% using a Reynolds decomposition
% November 3, 2013
% Sam Stevenson
% (modified from advection_ml.m)
% FULLY DEBUGGED AND VALIDATED JANUARY 2021 - see
% ~/RAPIDpaper/isoratio_budget_test.pdf for details

% Function assumes that the input u, v, and S (tracer) arrays have dimensions (time
% x lat x lon) and that velocities are given in units of m/s; time
% is assumed to have units of seconds
% lat, lon are assumed to follow ROMS grid conventions: 2D matrices, units
% of degrees
% and the overall grid structure is assumed to be an Arakawa C-grid

function [umdSmdx,updSmdx,umdSpdx,updSpdx,vmdSmdy,vpdSmdy,vmdSpdy,vpdSpdy,dSpdt,mnupdSpdx,mnvpdSpdy]=advection_ml_rd(salt,uvel,vvel,lat,lon,time,yr,mon,yrclim)
    RE=6378.e3;

    % Calculate horizontal tracer gradient 
    % dS/dx
    latarr=lat*pi/180.;
    lonarr=lon.*pi/180.;

    dSdx=zeros(size(salt));
    dSdy=zeros(size(salt));
    dSpdt=zeros(size(salt));

    for tt=1:size(salt,1)
        dS=squeeze(salt(tt,:,[2 2:end end])-salt(tt,:,[1 1:(end-1) (end-1)]));
        dStmp=dS.*cos(latarr(:,[1 1:(end-1) (end-1)]))./(lonarr(:,[2 2:end end]) - lonarr(:,[1 1:(end-1) (end-1)]));
        dSdx(tt,:,:)=0.5*(dStmp(:,1:(end-1))+dStmp(:,2:end))./RE;

        dStmp=squeeze(salt(tt,[2 2:end end],:)-salt(tt,[1 1:(end-1) (end-1)],:))./(latarr([2 2:end end],:) - latarr([1 1:(end-1) (end-1)],:));
        dSdy(tt,:,:)=0.5*(dStmp(1:(end-1),:)+dStmp(2:end,:))./RE;
    end
    
    Sm=zeros(12,size(dSdx,2),size(dSdx,3));
    dSmdx=zeros(12,size(dSdx,2),size(dSdx,3));
    dSmdy=zeros(12,size(dSdy,2),size(dSdy,3));
    um=zeros(12,size(uvel,2),size(uvel,3));
    vm=zeros(12,size(vvel,2),size(vvel,3));
    umdSmdx=zeros(12,size(uvel,2),size(uvel,3));
    vmdSmdy=zeros(12,size(vvel,2),size(vvel,3));
    
    
    for mm=1:12
        thism=find(mon == mm & yr >= yrclim(1) & yr <= yrclim(2));
        
        % Climatological mean salinity 
        Sm(mm,:,:)=squeeze(nanmean(salt(thism,:,:),1));

        % Climatological mean currents
        um(mm,:,:)=squeeze(nanmean(uvel(thism,:,:),1));
        vm(mm,:,:)=squeeze(nanmean(vvel(thism,:,:),1));
    end
    
    
    
    for mm=1:12        
        thism=find(mon == mm);

        % Gradient of climatological-mean salinity
        dS=squeeze(Sm(mm,:,[2 2:end end])-Sm(mm,:,[1 1:(end-1) (end-1)]));
        dStmp=dS.*cos(latarr(:,[1 1:(end-1) (end-1)]))./(lonarr(:,[2 2:end end]) - lonarr(:,[1 1:(end-1) (end-1)]));
        dSmdx(mm,:,:)=0.5*(dStmp(:,1:(end-1))+dStmp(:,2:end))./RE;

        dStmp=squeeze(Sm(mm,[2 2:end end],:)-Sm(mm,[1 1:(end-1) (end-1)],:))./(latarr([2 2:end end],:) - latarr([1 1:(end-1) (end-1)],:));
        dSmdy(mm,:,:)=0.5*(dStmp(1:(end-1),:)+dStmp(2:end,:))./RE;
        
        % Climatological mean advection
        umdSmdx(thism,:,:)=repmat(um(mm,:,:).*dSmdx(mm,:,:),[length(thism) 1 1]);
        vmdSmdy(thism,:,:)=repmat(vm(mm,:,:).*dSmdy(mm,:,:),[length(thism) 1 1]);
    end
    
    % Anomalous advection of climatological mean gradient
    updSmdx=zeros(size(salt));
    vpdSmdy=zeros(size(salt));
    
    for tt=1:length(time)
        mym=mon(tt);
        updSmdx(tt,:,:)=(uvel(tt,:,:)-um(mym,:,:)).*dSmdx(mym,:,:);
        vpdSmdy(tt,:,:)=(vvel(tt,:,:)-vm(mym,:,:)).*dSmdy(mym,:,:);
    end
    
    % Climatological mean advection of anomalous gradient
    umdSpdx=zeros(size(salt));
    vmdSpdy=zeros(size(salt));
    
    for tt=1:length(time)
       mym=mon(tt);
       umdSpdx(tt,:,:)=um(mym,:,:).*(dSdx(tt,:,:)-dSmdx(mym,:,:)); 
       vmdSpdy(tt,:,:)=vm(mym,:,:).*(dSdy(tt,:,:)-dSmdy(mym,:,:)); 
    end
    
    % Anomalous advection of anomalous gradient
    updSpdx=zeros(size(salt));
    vpdSpdy=zeros(size(salt));
    
    for tt=1:length(time)
        mym=mon(tt);
        updSpdx(tt,:,:)=(uvel(tt,:,:)-um(mym,:,:)).*(dSdx(tt,:,:)-dSmdx(mym,:,:));
        vpdSpdy(tt,:,:)=(vvel(tt,:,:)-vm(mym,:,:)).*(dSdy(tt,:,:)-dSmdy(mym,:,:));
    end

    % Climatological mean of eddy terms
    mnupdSpdx=zeros(size(salt));
    mnvpdSpdy=zeros(size(salt));
    for mm=1:12
        mym=find(mon == mm);
        mnupdSpdx(mym,:,:)=repmat(nanmean(updSpdx(mym,:,:),1),[length(mym) 1 1]);
        mnvpdSpdy(mym,:,:)=repmat(nanmean(vpdSpdy(mym,:,:),1),[length(mym) 1 1]);
    end    
    
    % Time derivative: use anomaly
    for tt=1:size(salt,1)
       mym=mon(tt);
       salt(tt,:,:)=salt(tt,:,:)-Sm(mym,:,:);
    end
    
    if size(salt,1) > 1
        for la=1:size(salt,2)
            for lo=1:size(salt,3)
                dStmp=squeeze(salt([2 2:end end],la,lo)-salt([1 1:(end-1) (end-1)],la,lo))./(time([2 2:end end])-time([1 1:(end-1) (end-1)]));
                dSpdt(:,la,lo)=0.5*(dStmp(1:(end-1))+dStmp(2:end));
            end
        end
    end
    
end
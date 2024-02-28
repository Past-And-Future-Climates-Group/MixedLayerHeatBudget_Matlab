% Compute terms in heat budget related to vertical advection, assuming a
% spatially and temporally varying mixed-layer depth

function [w_entr,wmdTmdz,wpdTmdz,wmdTpdz,wpdTpdz,mnwpdTpdz]=vertadv_ml_rd(time,lat,lon,mld,umld,vmld,Tmld,Tsub,wsub,mon,yr,yrclim)
    RE=6378.e3;

    % Compute entrainment velocity
    dHdt=zeros(size(mld));
    dHdx=zeros(size(mld));
    dHdy=zeros(size(mld));
    latarr=lat*pi/180.;
    lonarr=lon.*pi/180.;

    % Time, horizontal derivatives of MLD
    for la=1:size(mld,2)
        for lo=1:size(mld,3)
            dHtmp=squeeze(mld([2 2:end end],la,lo)-mld([1 1:(end-1) (end-1)],la,lo))./(time([2 2:end end])-time([1 1:(end-1) (end-1)]));
            dHdt(:,la,lo)=0.5*(dHtmp(1:(end-1))+dHtmp(2:end));  % m/day
        end
    end

    for tt=1:size(mld,1)
        dH=squeeze(mld(tt,:,[2 2:end end])-mld(tt,:,[1 1:(end-1) (end-1)]));
        dHtmp=dH.*cos(latarr(:,[1 1:(end-1) (end-1)]))./(lonarr(:,[2 2:end end]) - lonarr(:,[1 1:(end-1) (end-1)]));
        dHdx(tt,:,:)=0.5*(dHtmp(:,1:(end-1))+dHtmp(:,2:end))./RE;   % m/m

        dHtmp=squeeze(mld(tt,[2 2:end end],:)-mld(tt,[1 1:(end-1) (end-1)],:))./(latarr([2 2:end end],:) - latarr([1 1:(end-1) (end-1)],:));
        dHdy(tt,:,:)=0.5*(dHtmp(1:(end-1),:)+dHtmp(2:end,:))./RE;
    end        

    % entrainment velocity
    w_entr=dHdt/86400.+umld.*dHdx+vmld.*dHdy+wsub;


    %%%
    % Apply Reynolds decomposition to entrainment velocity,
    % cross-MLD temperature difference
    dTdz=(Tmld-Tsub)./mld;
    dTm=zeros(12,size(dTdz,2),size(dTdz,3));
    Tm_mld=zeros(12,size(dTdz,2),size(dTdz,3));
    Tm_sub=zeros(12,size(dTdz,2),size(dTdz,3));
    wm=zeros(12,size(dTdz,2),size(dTdz,3));
    mnwpdTpdz=zeros(12,size(dTdz,2),size(dTdz,3));

    % Compute climatological means
    for mm=1:12
        thism=find(mon == mm & yr >= yrclim(1) & yr <= yrclim(2));

        % Climatological mean vertical temperature gradient
        Tm_mld(mm,:,:)=squeeze(nanmean(Tmld(thism,:,:),1));
        Tm_sub(mm,:,:)=squeeze(nanmean(Tsub(thism,:,:),1));

        % Climatological mean entrainment velocity
        wm(mm,:,:)=squeeze(nanmean(w_entr(thism,:,:),1));

        % Climatological mean temperature difference
        dTm(mm,:,:)=(Tm_mld(mm,:,:)-Tm_sub(mm,:,:));
    end
    
    % Heaviside step function on climatological mean entrainment
    % velocity
    wsgn=wm;
    wsgn(wsgn <= 0)=0;   
    wsgn(wsgn ~= 0)=1;

    wmdTmdz=zeros(size(Tmld));
    wpdTmdz=zeros(size(Tmld));
    wmdTpdz=zeros(size(Tmld));
    wpdTpdz=zeros(size(Tmld));           


    for tt=1:length(time)
        mym=mon(tt);
        
        % Use climatological-mean temperature difference, along with temporally
        % varying MLD, to calculate climatological-mean vertical advection
        wmdTmdz(tt,:,:)=wsgn(mym,:,:).*wm(mym,:,:).*dTm(mym,:,:)./mld(tt,:,:);

        % Anomalous advection of climatological mean gradient
        wtmp=(w_entr(tt,:,:)-wm(mym,:,:));
        wpdTmdz(tt,:,:)=wsgn(mym,:,:).*wtmp.*(Tm_mld(mym,:,:)-Tm_sub(mym,:,:))./mld(tt,:,:);

        % Climatological mean advection of anomalous gradient
        tmp1=Tmld(tt,:,:)-Tm_mld(mym,:,:);
        tmp2=Tsub(tt,:,:)-Tm_sub(mym,:,:);
        wmdTpdz(tt,:,:)=wsgn(mym,:,:).*wm(mym,:,:,:).*(tmp1-tmp2)./mld(tt,:,:);

        % Anomalous advection of anomalous gradient
        wpdTpdz(tt,:,:)=wsgn(mym,:,:).*wtmp.*(tmp1-tmp2)./mld(tt,:,:);
    end

    for mm=1:12
       thism=find(mon == mm & yr >= yrclim(1) & yr <= yrclim(2)); 

        % Average of anomalous advection of anomalous gradient
        mnwpdTpdz(mm,:,:)=nanmean(wpdTpdz(thism,:,:),1);
    end

end
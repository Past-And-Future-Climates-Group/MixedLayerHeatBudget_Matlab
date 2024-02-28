% Compute mixed-layer heat budget for LME simulations, figure out why there
% is a change in the type of El Nino event
% March 2015
% Sam Stevenson

cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];  % colorbar
regbox=[-10,10,90,300];     % Region of interest
RE=6378.e3;
rho=1025;       % Mean density of seawater at 20C, 35 psu: kg/m^2
cp=3993;        % Specific heat of seawater at 20C, 35 psu, J/kg/K

% % CCSM4 control
% ensnames={'CCSM4 PI'};
% runnames={'b40.1850.track1.1deg.006'};
% dates={'080001-089912'};

%%%%%%% Last Millennium case %%%%%%
% Names of LM runs
ensnames={'Full','GHG','LULC','Orbital','Solar','Volcanic','OzoneAer','Control'};
runnames={'b.e11.BLMTRC5CN.f19_g16.001','b.e11.BLMTRC5CN.f19_g16.002','b.e11.BLMTRC5CN.f19_g16.003','b.e11.BLMTRC5CN.f19_g16.004', ...
   'b.e11.BLMTRC5CN.f19_g16.005','b.e11.BLMTRC5CN.f19_g16.006','b.e11.BLMTRC5CN.f19_g16.007', ...
   'b.e11.BLMTRC5CN.f19_g16.008','b.e11.BLMTRC5CN.f19_g16.009','b.e11.BLMTRC5CN.f19_g16.010', ...
   'b.e11.BLMTRC5CN.f19_g16.011','b.e11.BLMTRC5CN.f19_g16.012','b.e11.BLMTRC5CN.f19_g16.013';
   'b.e11.BLMTRC5CN.f19_g16.GHG.001','b.e11.BLMTRC5CN.f19_g16.GHG.002','b.e11.BLMTRC5CN.f19_g16.GHG.003', ...
   '','','','','','','','','','';
   'b.e11.BLMTRC5CN.f19_g16.LULC_HurttPongratz.001','b.e11.BLMTRC5CN.f19_g16.LULC_HurttPongratz.002','b.e11.BLMTRC5CN.f19_g16.LULC_HurttPongratz.003', ...
   '','','','','','','','','','';
   'b.e11.BLMTRC5CN.f19_g16.ORBITAL.001','b.e11.BLMTRC5CN.f19_g16.ORBITAL.002','b.e11.BLMTRC5CN.f19_g16.ORBITAL.003', ...
   '','','','','','','','','','';
   'b.e11.BLMTRC5CN.f19_g16.SSI_VSK_L.001','b.e11.BLMTRC5CN.f19_g16.SSI_VSK_L.003','b.e11.BLMTRC5CN.f19_g16.SSI_VSK_L.004','b.e11.BLMTRC5CN.f19_g16.SSI_VSK_L.005', ...
   '','','','','','','','','';
   'b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.001','b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.002','b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.003', ...
   'b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.004','b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.005','','','','','','','','';
   'b.e11.BLMTRC5CN.f19_g16.OZONE_AER.001','b.e11.BLMTRC5CN.f19_g16.OZONE_AER.002','b.e11.BLMTRC5CN.f19_g16.OZONE_AER.003', ...
   'b.e11.BLMTRC5CN.f19_g16.OZONE_AER.004','b.e11.BLMTRC5CN.f19_g16.OZONE_AER.005','','','','','','','','';
   'b.e11.B1850C5CN.f19_g16.0850cntl.001','','','','','','','','','','','','';
   };

% ensnames={'Full'};
% runnames={'b.e11.BLMTRC5CN.f19_g16.001','b.e11.BLMTRC5CN.f19_g16.011','b.e11.BLMTRC5CN.f19_g16.012','b.e11.BLMTRC5CN.f19_g16.013'};
dates={'085001-089912','090001-099912','100001-109912','110001-119912','120001-129912','130001-139912',...
   '140001-149912','150001-159912','160001-169912','170001-179912','180001-184912','185001-200512'};


%sstnc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/SST/b.e11.BLMTRC5CN.f19_g16.002.pop.h.SST.180001-184912.nc'));
sstnc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/b40.1850.track1.1deg.006.pop.h.HMXL.080001-089912.nc'));
tlat=sstnc{'TLAT'}(:,:);
tlon=sstnc{'TLONG'}(:,:);
mylat=find(tlat(:,150) >= regbox(1) & tlat(:,150) <= regbox(2));
mylon=find(tlon(150,:) >= regbox(3) & tlon(150,:) <= regbox(4));
tlat=tlat(mylat,mylon);
tlon=tlon(mylat,mylon);

%%%% Loop through each simulation
%for ee=1:length(ensnames)
for ee=6:6
    %for rr=1:size(runnames,2)
    for rr=5:5
        runname=runnames{ee,rr}
        if ~isempty(runname)
        
        for dd=1:length(dates)
            dates{dd}
            % Read in temperature, velocities in upper ocean
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/TEMP/',runname,'.pop.h.TEMP.',dates{dd},'.nc'));  % C
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.TEMP.',dates{dd},'.nc'));  % C
            temp=nc{'TEMP'}(:,1:20,mylat,mylon);
            time=nc{'time'}(:);
            [yr,mon,~]=datenumnoleap(time-29,[0 1 1]);
            
            z=nc{'z_t'}(1:20)/100.;
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/UVEL/',runname,'.pop.h.UVEL.',dates{dd},'.nc'));  % cm/s
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.UVEL.',dates{dd},'.nc'));  % cm/s
            uvel=nc{'UVEL'}(:,1:20,mylat,mylon)/100.;
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/VVEL/',runname,'.pop.h.VVEL.',dates{dd},'.nc'));
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.VVEL.',dates{dd},'.nc'));
            vvel=nc{'VVEL'}(:,1:20,mylat,mylon)/100.;            
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/WVEL/',runname,'.pop.h.WVEL.',dates{dd},'.nc'));
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.WVEL.',dates{dd},'.nc'));
            wvel=nc{'WVEL'}(:,1:20,mylat,mylon)/100.;    
            
            % Read in heat fluxes
            % SHF = "Total Surface Heat Flux, Including SW"
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/SHF/',runname,'.pop.h.SHF.',dates{dd},'.nc'));  % W/m^2
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.SHF.',dates{dd},'.nc'));  % W/m^2
            qnet=nc{'SHF'}(:,mylat,mylon);
            % SHF_QSW = "Solar Short-Wave Heat Flux"
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/SHF_QSW/',runname,'.pop.h.SHF_QSW.',dates{dd},'.nc'));  % W/m^2
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.SHF_QSW.',dates{dd},'.nc'));  % W/m^2
            qsw=nc{'SHF_QSW'}(:,mylat,mylon);
            
            % Use POP-computed mixed layer depth
            nc=netcdf(strcat('/glade/p/cesm0005/CESM-CAM5-LME/ocn/proc/tseries/monthly/HMXL/',runname,'.pop.h.HMXL.',dates{dd},'.nc'));  % cm
            %nc=netcdf(strcat('/glade/scratch/samantha/b40.1850.track1.1deg.006/',runname,'.pop.h.HMXL.',dates{dd},'.nc'));  % cm
            mld=nc{'HMXL'}(:,mylat,mylon)/100.;
            %mld=65+zeros(size(qsw));    % overwrite with fixed MLD
            
            % Heat flux term
            sfcflx=zeros(size(mld));
            qpen=zeros(size(mld));

            for tt=1:length(time)
                qpen(tt,:,:)=qsw(tt,:,:).*(0.58*exp(-1/0.35*mld(tt,:,:)) + 0.42*exp(-1/23.*mld(tt,:,:)));
                sfcflx(tt,:,:)=(1./(rho*cp*mld(tt,:,:))).*(qnet(tt,:,:)-qpen(tt,:,:));
            end
                                   
            % Integrate temperature, velocities above mixed layer
            pacz=repmat(z,[1 size(temp,3) size(temp,4)]);
            Tmld=mldavg_varytime(mld,temp,time,pacz);  
            Tmld(abs(Tmld) > 1e10)=0/0;

            umld=mldavg_varytime(mld,uvel,time,pacz);
            vmld=mldavg_varytime(mld,vvel,time,pacz);
            
            % Get temperature, velocities just below mixed layer
            Tsub=submld_varytime(mld,temp,time,pacz,'first');
            usub=submld_varytime(mld,uvel,time,pacz,'first');
            vsub=submld_varytime(mld,vvel,time,pacz,'first');
            wsub=submld_varytime(mld,wvel,time,pacz,'first');
            
            
            % Use Reynolds decomposition version of the horizontal advection formulation
            % NB: "mean" is a repeating 12-month climatology for all budget terms
            [umdTmdx,updTmdx,umdTpdx,updTpdx,vmdTmdy,vpdTmdy,vmdTpdy,vpdTpdy,dTdt,mnupdTpdx,mnvpdTpdy]=advection_ml_rd(Tmld,umld,vmld,tlat,tlon,time,yr,mon,[min(yr),max(yr)]);

            %%% Vertical advective terms
            [w_entr,wmdTmdz,wpdTmdz,wmdTpdz,wpdTpdz,mnwpdTpdz]=vertadv_ml_rd(time,tlat,tlon,mld,usub,vsub,Tmld,Tsub,wsub,mon,yr,[min(yr),max(yr)]);         
            
            clear temp
            clear uvel
            clear vvel
            clear wvel
            %clear mld
            clear qnet
            clear qsw
            clear vmld
            
            % Save to a .mat file for use later
            %save(strcat('/glade/scratch/samantha/LastMillennium/',runname,'_heatbudget_ml_rd_',dates{dd},'_varyMLD_test.mat'))
            save(strcat('/glade/p/cesm/palwg_dev/LME/proc/samantha/LME_heatbudget/',runname,'_heatbudget_ml_rd_',dates{dd},'_varyMLD.mat'))
            
        end
        end
    end
end

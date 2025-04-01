function [pre_S5_PSR,Tar_Vs_TopUC,Tar_Depth_TopUC]=sub_WRR_derived_PSRmodel(Tar_Lon,Tar_Lat,WRR,NZbor,Flag_Vs_TUC,Vs_TUC,Asp,rake);
%%%%%%
%     Width-Rupture-ratio (WRR) derived Probabilistic model of surface rupture (PSR)
%     *****Inputs*****
%     Tar_Lon:      Epicenter longitude; Tar_Lat: Epicenter latitute. (degree)
%     WRR:          The Width-Rupture-ratio is defined as the rupture width (W_{Rup}) 
%                   divided by the fault width (W_{Flt}). It can utilize the W_{Rup} 
%                   from an empirical scaling model that incorporates the width-limit 
%                   effect (for example, Huang et al., SRL, 95, 2352-2367, 2024).
%     NZbor:        The normalized depth of the bottom of the fault rupture within
%                   the seismogenic zone.
%     Flag_Vs_TUC:  A flag indicating the source of the source zone stiffness index Vs, 
%                   where 1 represents using the Vs at the Top Upper Crust from Simmons 
%                   et al. (2015) (default), and 0 indicates customizing the given Vs at 
%                   the Top Upper Crust, for example, manual_Vs_TUC=3.5 km/s. The [Vs_TUC] 
%                   (km/s) values are defaulted as grabbed from 
%                   "Global_Vs_Structure_grab_from_LLNL_G3D_JPS_by_Sea15.mat," which is 
%                   reformatted from LLNL_G3D_JPS.Interpolated files from 
%                   https://gs.llnl.gov/nuclear-threat-reduction/nuclear-explosion-monitoring/global-3d-seismic-tomography. 
%                   More information can be found in Simmons et al. (GRL, 42(21), 
%                   9270-9278, 2015).
%     Vs_TUC:       Manual given Vs_{TUC} opporated only if [Flag_Vs_TUC]=0. (km/s)
%     Asp:          The seismic fault aspect ratio, defined as the rupture length divided 
%                   by the rupture width can be obtained using an empirical source scaling 
%                   model that considers the width-limit effect (e.g., Huang et al., SRL, 95, 
%                   2352-2367, 2024).
%     rake:         Rupture rake angle.
%     *****Outputs*****
%     pre_S5_PSR:      Probability of surface rupture derived on the stage 5 model 
%                      provided in this study (Huang and Abrahamson, submitted)
%     Tar_Vs_TopUC:    The epicenter location's [Vs_TUC] (km/s) value grabbed in 
%                      Simmons et al. (GRL, 42(21), 9270-9278, 2015)
%     Tar_Depth_TopUC: The depth of the location of [Vs_TUC] is grabbed in Simmons et 
%                      al. (GRL, 42(21), 9270-9278, 2015). 
%                      These two values will be obtained only if [Flag_Vs_TUC]=1. 
%
%                                           by Bob J.Y. Huang 2025.03.31
%%%%%%
load Global_Vs_Structure_grab_from_LLNL_G3D_JPS_by_Sea15.mat; 
if(Flag_Vs_TUC==1)
  [Tar_Vs_TopUC,Tar_Depth_TopUC]=sub_getVs_TUC(Tar_Lon,Tar_Lat,Geodetic_longitude,Geodetic_latitute,Dep_L5,Vs_L5,Dep_L6,Vs_L6,Dep_L7,Vs_L7,Dep_L8,Vs_L8,Dep_L9,Vs_L9,Dep_L10,Vs_L10,Dep_L11,Vs_L11,Dep_L12,Vs_L12,Dep_L13,Vs_L13,Dep_L14,Vs_L14,Dep_L15,Vs_L15,Dep_L16,Vs_L16);
  disp(['Vs Top Upper Crust: ',num2str(Tar_Vs_TopUC),' km/s at depth of ',num2str(Tar_Depth_TopUC),' km collected by Simmons et al. (2015)']);
elseif(Flag_Vs_TUC==0)
  if(Vs_TUC<=0)
    pre_S5_PSR=-999;Tar_Vs_TopUC=-999;Tar_Depth_TopUC=-999;
    disp(['Error in Vs Top Upper Crust, which is customized given as: ',num2str(Vs_TUC),' km/s']);
    return
  elseif(Vs_TUC>2.8&Vs_TUC<3.8)
    Tar_Vs_TopUC=Vs_TUC;Tar_Depth_TopUC=-999;
    disp(['Vs Top Upper Crust is customized given as: ',num2str(Vs_TUC),' km/s when Flag_Vs_TUC=',num2str(Flag_Vs_TUC)]);
  else
    disp(['!!Noted!! Vs Top Upper Crust is suggested in the data range of 2.87 to 3.723 km/s, which is now customized given as: ',num2str(Vs_TUC),' km/s, the VsTUC scaling calculation is extrapolated.']);
    Tar_Vs_TopUC=Vs_TUC;Tar_Depth_TopUC=-999;
  end
else
  pre_S5_PSR=-999;Tar_Vs_TopUC=-999;Tar_Depth_TopUC=-999;
  disp(['!!Warning!! Flag_Vs_TUC is given in a wrong value. 1 = collected Vs at Top Upper Crust by Simmons et al. (2015); 0 = use customized value for Vs_TUC']);
  return
end   
% Coefficients
c0=-11.2346;
c1=11.9612;
c2=-8.5082;
c3=2.8917;
% Functional forms of logistic regression
S3_y=c0+c1*WRR+c2*NZbor+c3*Tar_Vs_TopUC;
pre_S3_PSR=(exp(S3_y))/(1+(exp(S3_y)));
% Functional forms for the Width limit term conducted from the fault aspect ratio
c4=-0.1297; 
c5=0.0701; 
c6=0.2172; 
Asp_start=0;Asp_end=5;
if(Asp>Asp_start&Asp<Asp_end)
  f1_Asp=c4+c5*(Asp);
elseif(Asp>=Asp_end)
  f1_Asp=c6;
end
% Style of Faulting term for normal (NML) fault only
%% classify style of fault by rake angles follows Ancheta et al. 2013 (see table 4.4 in PEER report 2020/02
%  "Mechanism"      Rake angles                  abbreviations
%  Strike-Slip      -180~-150; -30~30; 150~180   SS
%  Normal           -150~-30                     N
%  Reverse           30~150                      R
c7=0.1713;
if(rake>=-150&rake<-30)
  pre_S5_PSR=pre_S3_PSR+f1_Asp+c7;
else
  pre_S5_PSR=pre_S3_PSR+f1_Asp;
end
if(pre_S5_PSR>1)
  pre_S5_PSR=1;
end
if(pre_S5_PSR<0)
  pre_S5_PSR=0;
end


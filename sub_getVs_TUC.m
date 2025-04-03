function [Vs_TopUC,Depth_TopUC]=sub_getVs_TUC(Tar_Lon,Tar_Lat,Geodetic_longitude,Geodetic_latitute,Dep_L5,Vs_L5,Dep_L6,Vs_L6,Dep_L7,Vs_L7,Dep_L8,Vs_L8,Dep_L9,Vs_L9,Dep_L10,Vs_L10,Dep_L11,Vs_L11,Dep_L12,Vs_L12,Dep_L13,Vs_L13,Dep_L14,Vs_L14,Dep_L15,Vs_L15,Dep_L16,Vs_L16)
%----
%    Global Vs structure are from LLNL_G3D_JPS.Interpolated from https://gs.llnl.gov/nuclear-threat-reduction/nuclear-explosion-monitoring/global-3d-seismic-tomography
%    Detailed could be found in Simmons et al. (2015, denoted as Sea15)
%    Layer 5 : Top_of_upper_sediment
%          6 : Bottom_of_upper_sediment
%          7 : Top_of_middle_sediment
%          8 : Bottom_of_middle_sediment
%          9 : Top_of_lower_sediment
%          10: Bottom_of_lower_sediment
%          11: Top_of_upper_crust
%          12: Bottom_of_upper_crust
%          13: Top_of_middle_crust
%          14: Bottom_of_middle_crust
%          15: Top_of_lower_crust
%          16: Bottom_of_lower_crust
%    A Global_Vs_Structure_grab_from_LLNL_G3D_JPS_by_Sea15.csv and mat file are 
%      gathered and rearranged from the Sea15 for easier usage. Inputs parameters 
%      excepted target coordinates (Tar_Lon, Tar_Lat) are all imported from the mat 
%      file.
%    The reason to keep it complex for whole Vs structure is in case there are needs 
%      for results of other layers.
%                                           by Bob J.Y. Huang 2025.03.31
%----
diff_x=abs(Geodetic_longitude-Tar_Lon);
diff_y=abs(Geodetic_latitute-Tar_Lat);
index_TarRegion=find(diff_x==min(diff_x)&diff_y==min(diff_y));
icount=0;
for i=5:16
  icount=icount+1;
  eval(['div_depth(icount)=Dep_L',num2str(i),'(index_TarRegion(1));']);
  eval(['div_Vs(icount)=Vs_L',num2str(i),'(index_TarRegion(1));']);
end
Depth_TopUC=div_depth(7); % Layer of Top Upper Crust
Vs_TopUC=div_Vs(7); % Layer of Top Upper Crust
Vs_TopUC=round(Vs_TopUC*1000)/1000;

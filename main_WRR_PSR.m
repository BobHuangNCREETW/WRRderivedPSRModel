%%%%%
%    This is the example MATLAB program to run the WRR-derived PSR model developed in
%    Huang and Abrahamson (submitted) applying subroutine [sub_WRR_derived_PSRmodel.m]
%    Detaled explaination of parameter used is given in [sub_WRR_derived_PSRmodel.m].
%                                           by Bob J.Y. Huang 2025.03.31
%%%%%
% Parameter setting for case 1
example_Lon_CA=-122.235;example_Lat_CA=37.858; % A ramdom place in California
examepl_Wrup=15; % Assumed rupture width
example_dip_mean_CA=65.7; % Assumed fault dip angle
example_ZBSZ_CA=19.7; % Assumed bottom depth of the seismogenic zone
example_Wflt=example_ZBSZ_CA/sind(example_dip_mean_CA); % Assumed fault width considering fault dip and ZBSZ
example_WRR=examepl_Wrup/example_Wflt; % Assumed Width-Rupture ratio
example_NZbor=example_WRR+(1-example_WRR)/2; % Assumed normalized bottom depth of the fault rupture using Eq. (8) in Huang and Abrahamson 
Flag_Vs_TUC=1;manual_Vs_TUC=-999; % Grabbed the Vs_TUC from default Sea15
example_Asp=5; % Assumed aspect ratio
example_rake=0; % Assumed Strike-Slip fault
[pre_S5_PSR_case1,VsTUC_case1,DatTUC_case1]=sub_WRR_derived_PSRmodel(example_Lon_CA,example_Lat_CA,example_WRR,example_NZbor,Flag_Vs_TUC,manual_Vs_TUC,example_Asp,example_rake)
% Parameter setting for case 2
example_Lon_SEU=13.5;example_Lat_SEU=42; % A ramdom place in South Europe
examepl_Wrup=15;
example_dip_mean_SEU=50;
example_ZBSZ_SEU=16;
example_Wflt=example_ZBSZ_SEU/sind(example_dip_mean_SEU);
example_WRR=examepl_Wrup/example_Wflt;
example_NZbor=example_WRR+(1-example_WRR)/2;
Flag_Vs_TUC=1;manual_Vs_TUC=-999;
example_Asp=3;
example_rake=-90; % Assumed Normal fault
[pre_S5_PSR_case2,VsTUC_case2,DatTUC_case2]=sub_WRR_derived_PSRmodel(example_Lon_SEU,example_Lat_SEU,example_WRR,example_NZbor,Flag_Vs_TUC,manual_Vs_TUC,example_Asp,example_rake)
% Parameter settubg for case 3
example_Lon_TW=121.4;example_Lat_TW=24.9; % A ramdom place in Taiwan
examepl_Wrup=25;
example_dip_mean_TW=48;
example_ZBSZ_TW=26.5;
example_Wflt=example_ZBSZ_TW/sind(example_dip_mean_TW);
example_WRR=examepl_Wrup/example_Wflt;
example_NZbor=example_WRR+(1-example_WRR)/2;
Flag_Vs_TUC=0;manual_Vs_TUC=3.2; % Manually given the VsTUC parameter (km/s)
example_Asp=1.5;
example_rake=90; % Assumed Reverse fault
[pre_S5_PSR_case3,VsTUC_case3,DatTUC_case3]=sub_WRR_derived_PSRmodel(example_Lon_TW,example_Lat_TW,example_WRR,example_NZbor,Flag_Vs_TUC,manual_Vs_TUC,example_Asp,example_rake)





clc,clear
%%%计算收敛阶
%%%Example1_m1
del_Ex1m1=[1.5960,1.2201,0.7958];
err_Ex1m1_adp = [0.0649, 0.1463, 0.1404;
     0.0356, 0.0810, 0.0786;
     0.0038, 0.0171, 0.0169];
err_Ex1m1_opt= [0.1661, 0.4940, 0.4920;
     0.1384, 0.4158, 0.4151;
     0.0975, 0.2874, 0.2867];
r_Ex1m1_adp=showrate(del_Ex1m1,err_Ex1m1_adp(:,1),err_Ex1m1_adp(:,2),err_Ex1m1_adp(:,3))
r_Ex1m1_opt=showrate(del_Ex1m1,err_Ex1m1_opt(:,1),err_Ex1m1_opt(:,2),err_Ex1m1_opt(:,3))
%%%Example2_m1
del_Ex2m1=[0.8062,0.6001,0.3985];
err_Ex2m1_adp= [0.0277, 0.0597, 0.0587;
     0.0118, 0.0265, 0.0256;
     0.0021, 0.0093, 0.0083];
err_Ex2m1_opt= [0.0940, 0.2543, 0.2541;
     0.0785, 0.2047, 0.2047;
     0.0557, 0.1413, 0.1412];
r_Ex2m1_adp=showrate(del_Ex2m1,err_Ex2m1_adp(:,1),err_Ex2m1_adp(:,2),err_Ex2m1_adp(:,3))
r_Ex2m1_opt=showrate(del_Ex2m1,err_Ex2m1_opt(:,1),err_Ex2m1_opt(:,2),err_Ex2m1_opt(:,3))
%%%Example3_m1
del_Ex3m1=[0.3138,0.2081,0.1508];
err_Ex3m1_adp= [0.0118, 0.0241, 0.0240;
     0.0084, 0.0159, 0.0159;
     0.0011, 0.0022, 0.0022];
err_Ex3m1_opt= [0.0233, 0.0672, 0.0672;
     0.0168, 0.0487, 0.0487;
     0.0094, 0.0251, 0.0251];
r_Ex3m1_adp=showrate(del_Ex3m1,err_Ex3m1_adp(:,1),err_Ex3m1_adp(:,2),err_Ex3m1_adp(:,3))
r_Ex3m1_opt=showrate(del_Ex3m1,err_Ex3m1_opt(:,1),err_Ex3m1_opt(:,2),err_Ex3m1_opt(:,3))
%%%Example1_m2
del_Ex1m2=[1.6253,1.2002,0.7940];
err_Ex1m2_adp= [0.0640, 0.2216, 0.1479;
     0.0351, 0.1153, 0.0788;
     0.0060, 0.0316, 0.0249];
err_Ex1m2_opt= [0.1671, 0.7697, 0.5131;
     0.1382, 0.6623, 0.4196;
     0.0979, 0.4468, 0.2849];
r_Ex1m2_adp=showrate(del_Ex1m2,err_Ex1m2_adp(:,1),err_Ex1m2_adp(:,2),err_Ex1m2_adp(:,3))
r_Ex1m2_opt=showrate(del_Ex1m2,err_Ex1m2_opt(:,1),err_Ex1m2_opt(:,2),err_Ex1m2_opt(:,3))
%%%Example2_m2
del_Ex2m2=[0.8136,0.6162,0.3998];
err_Ex2m2_adp= [0.0276, 0.0935, 0.0602;
     0.0119, 0.0439, 0.0280;
     0.0038, 0.0160, 0.0118];
err_Ex2m2_opt= [0.0945, 0.3910, 0.2571;
     0.0786, 0.3225, 0.2103;
     0.0560, 0.2181, 0.1469];
r_Ex2m2_adp=showrate(del_Ex2m2,err_Ex2m2_adp(:,1),err_Ex2m2_adp(:,2),err_Ex2m2_adp(:,3))
r_Ex2m2_opt=showrate(del_Ex2m2,err_Ex2m2_opt(:,1),err_Ex2m2_opt(:,2),err_Ex2m2_opt(:,3))
%%%%Example3_m2
del_Ex3m2=[0.3187,0.2091,0.1507];
err_Ex3m2_adp = [
    0.0118, 0.0394, 0.0244;
    0.0083, 0.0258, 0.0166;
    0.0011, 0.0037, 0.0023
];
err_Ex3m2_opt= [
    0.0231, 0.1024, 0.0692;
    0.0169, 0.0707, 0.0477;
    0.0094, 0.0414, 0.0254
];
r_Ex3m2_adp=showrate(del_Ex3m2,err_Ex3m2_adp(:,1),err_Ex3m2_adp(:,2),err_Ex3m2_adp(:,3))
r_Ex3m2_opt=showrate(del_Ex3m2,err_Ex3m2_opt(:,1),err_Ex3m2_opt(:,2),err_Ex3m2_opt(:,3))
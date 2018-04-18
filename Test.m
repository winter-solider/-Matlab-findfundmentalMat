%% 实现了opencv中的八对坐标计算基本矩阵的方法；opencv 中的随机选择是伪随机方法，每次
%% 从均匀分布中选择八对，导致每次开始选择的都一样，

clear
clc
dbstop if error


kp1 = load('kp1.txt');
kp2 = load('kp2.txt');

[F, maxGoodCount] = ransac_fundamental_matrix(kp1, kp2,1,2000,0.9);
F
maxGoodCount






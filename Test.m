%% ʵ����opencv�еİ˶���������������ķ�����opencv �е����ѡ����α���������ÿ��
%% �Ӿ��ȷֲ���ѡ��˶ԣ�����ÿ�ο�ʼѡ��Ķ�һ����

clear
clc
dbstop if error


kp1 = load('kp1.txt');
kp2 = load('kp2.txt');

[F, maxGoodCount] = ransac_fundamental_matrix(kp1, kp2,1,2000,0.9);
F
maxGoodCount






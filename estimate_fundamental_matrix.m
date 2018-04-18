function [ F_matrix ] = estimate_fundamental_matrix(m1,m2)

%% 八组坐标计算基本矩阵
%% Mean Coordinates for Points a and b
C1 = mean(m1);  %C1 点的平均坐标，每一列的均值
Cu1 = C1(1,1);  % x坐标的均值
Cv1 = C1(1,2);  % y坐标的均值
    
C2 = mean(m2); %C2 点的平均坐标，每一列的均值
Cu2 = C2(1,1); % x坐标的均值
Cv2 = C2(1,2); % y坐标的均值
    
%Finding total Euclidean Distance of each point from Cu,Cv for each image
dist1 = 0; 
dist2 = 0;
for i = 1:size(m1,1)  %每个点离平均点的欧氏距离总和
    dist1 = dist1 + sqrt((m1(i,1)-Cu1)^2 + (m1(i,2)-Cv1)^2);
    dist2 = dist2 + sqrt((m2(i,1)-Cu2)^2 + (m2(i,2)-Cv2)^2);
end

dist1 = dist1/size(m1,1);  %%计算平均值
dist2 = dist2/size(m1,1);
 
%Finding scale transformation for point set a and b
s1 = sqrt(2)/(dist1);
s2 = sqrt(2)/(dist2);
    
%Transformation Matrix for Set of Points A and B
Ta = [s1, 0, -s1*Cu1; 
      0, s1, -s1*Cv1; 
      0, 0,   1];
Tb = [s2,  0,  -s2*Cu2; 
      0,  s2,  -s2*Cv2; 
      0,   0,   1]; 


%% Computation of Fundamental Matrix
A = ones(size(m1,1),9);

for i=1:size(m1,1)
    x1 = (m1(i, 1) - Cu1)*s1;
    y1 = (m1(i, 2)-Cv1)*s1;
    x2 = (m2(i, 1) - Cu2)*s2;
    y2 = (m2(i, 2) - Cv2)*s2;
    A(i, :) = [x2 * x1,  x2 * y1,  x2,  y2 * x1, y2*y1, y2, x1, y1,1]; 
end

[~, ~, V] = svd(A);
F_matrix = V(:,end);
F_matrix = reshape(F_matrix, [3 3])';
[U, S, V] = svd(F_matrix);
% set the smallest eigen value to zero to get rank 2 matrix

S(3,3) = 0;
S2 = S;

% re construct matrix with reduced rank
F_matrix = U * S2 * V';
F_matrix = Tb'*F_matrix*Ta;
F_matrix = F_matrix/F_matrix(3,3);  %可省略

end





















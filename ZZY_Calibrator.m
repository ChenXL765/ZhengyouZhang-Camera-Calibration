clear;clc;
N_image = 6;
size_square = 25;
N_point_x = 8;
N_point_y = 6;
N_point = N_point_x*N_point_y;

%% read image
%image coordinates
m = zeros(3,N_point,N_image);
for i = 1:N_image
    image = imread(['photo\',int2str(i),'.jpg']);
    m(1:2,:,i) = detectCheckerboardPoints(image)';
    m(3,:,i) = ones(1,N_point);
end
%world coordinates
M = zeros(3,N_point);
for i = 0:(N_point_x-1)
    for j = 0:(N_point_y-1)
        M(:,i*N_point_y+j+1) = [size_square*i;size_square*j;1];
    end
end

%% solve Homography
H = zeros(3,3,N_image);
for i = 1:N_image
    H(:,:,i) = Solve_Homography(m(:,:,i),M);
end

%% solve intrinsic matrix
V = zeros(2*N_image,6);
for i=1:N_image
    V(2*i-1,:)=[H(1,1,i)*H(1,2,i) H(1,1,i)*H(2,2,i)+H(2,1,i)*H(1,2,i) H(2,1,i)*H(2,2,i) H(3,1,i)*H(1,2,i)+H(1,1,i)*H(3,2,i) H(3,1,i)*H(2,2,i)+H(2,1,i)*H(3,2,i) H(3,1,i)*H(3,2,i)];
    v1=[H(1,1,i)^2 2*H(2,1,i)*H(1,1,i) H(2,1,i)^2 2*H(3,1,i)*H(1,1,i) 2*H(3,1,i)*H(2,1,i) H(3,1,i)^2];
    v2=[H(1,2,i)^2 2*H(2,2,i)*H(1,2,i) H(2,2,i)^2 2*H(3,2,i)*H(1,2,i) 2*H(3,2,i)*H(2,2,i) H(3,2,i)^2];
    V(2*i,:)=v1-v2;
end;
[u,s,v]=svd(V);
b = v(:,end);

v0    = (b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2);
lamda = b(6)-(b(4)^2+(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2)*(b(2)*b(4)-b(1)*b(5)))/b(1);
alpha = sqrt(lamda/b(1));
beta  = sqrt((lamda*b(1))/(b(1)*b(3)-b(2)^2));
c     = -(b(2)*alpha^2*beta)/lamda;
u0    = (c*(b(2)*v0))/alpha-(b(4)*alpha^2)/lamda;
A1 = [alpha c u0;0 beta v0;0 0 1];

param2 = [alpha c u0 beta v0];
options = optimset('Algorithm','levenberg-marquardt');
[x,res2] = lsqnonlin(@fun_2, param2 ,[],[], options, m, M, H);
res2 = res2/(N_image*N_point);
alpha = x(1);c = x(2);u0 = x(3);beta = x(4);v0 = x(5);
A2 = [alpha c u0;0 beta v0;0 0 1];
A = A2;

%% Radial Distortion
D = [];
d = [];
for i = 1:N_image
    RT = Solve_Extrinsic(A,H(:,:,i));
    XY = RT*M;
    UV = A*XY;
    XY=[XY(1,:)./XY(3,:); XY(2,:)./XY(3,:); XY(3,:)./XY(3,:)];
    UV=[UV(1,:)./UV(3,:); UV(2,:)./UV(3,:); UV(3,:)./UV(3,:)];
    for j = 1:N_point
        D = [D; [(UV(1,j)-u0)*( (XY(1,j))^2 + (XY(2,j))^2 ), (UV(1,j)-u0)*( (XY(1,j))^2 + (XY(2,j))^2 )^2] ;
                [(UV(2,j)-v0)*( (XY(1,j))^2 + (XY(2,j))^2 ), (UV(2,j)-v0)*( (XY(1,j))^2 + (XY(2,j))^2 )^2]];
        
        d = [d; (m(1,j,i)-UV(1,j)) ; (m(2,j,i)-UV(2,j))];
    end
end
k0 = inv(D'*D)*D'*d;
param3 = [alpha c u0 beta v0 k0'];
options = optimset('Algorithm','levenberg-marquardt');
[x,res3] = lsqnonlin(@fun_3, param3 ,[],[], options, m, M, H);
res3 = res3/(N_image*N_point);
alpha = x(1);c = x(2);u0 = x(3);beta = x(4);v0 = x(5);
A3 = [alpha c u0;0 beta v0;0 0 1];
A = A3;
k = [x(6); x(7)];
disp("相机内参为:");
disp(A);
disp("相机径向畸变系数为:");
disp(k);
disp("平均误差为：");
disp(res3);

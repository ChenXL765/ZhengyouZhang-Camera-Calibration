function f = fun_3(param, m, M, H)
A = [param(1) param(2) param(3); 0 param(4) param(5); 0 0 1];
k1 = param(6);
k2 = param(7);
u0 = param(3);
v0 = param(5);
N_point = size(m,2);
N_image = size(m,3);
f = [];
for i = 1:N_image
    lamda = (norm(inv(A)*H(:,1,i))+norm(inv(A)*H(:,2,i)))/2;
    r1 = 1/lamda*inv(A)*H(:,1,i);
    r2 = 1/lamda*inv(A)*H(:,2,i);
    r3 = cross(r1, r2);
    R = [r1 r2 r3];
    [u,s,v] = svd(R);
    R = u*v';
    t = 1/lamda*inv(A)*H(:,3,i);
    RT = [R(:,1) R(:,2) t];
    XY = RT*M;
    UV = A*XY;
    XY=[XY(1,:)./XY(3,:); XY(2,:)./XY(3,:); XY(3,:)./XY(3,:)];
    UV=[UV(1,:)./UV(3,:); UV(2,:)./UV(3,:); UV(3,:)./UV(3,:)];
    for j = 1:N_point
        UV(1,j) = UV(1,j) + (UV(1,j) - u0)*(k1*((XY(1,j))^2+(XY(2,j))^2) + k2*((XY(1,j))^2+(XY(2,j))^2)^2);
        UV(2,j) = UV(2,j) + (UV(2,j) - v0)*(k1*((XY(1,j))^2+(XY(2,j))^2) + k2*((XY(1,j))^2+(XY(2,j))^2)^2);
    end
    m_ = UV;
    m_ = [m_(1,:)./m_(3,:) ; m_(2,:)./m_(3,:); m_(3,:)./m_(3,:)];
    dm = m(:,:,i) - m_;
    f = [f dm(1,:) dm(2,:)];
end
end

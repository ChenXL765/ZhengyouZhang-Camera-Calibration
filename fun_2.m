function f = fun_2(param, m, M, H)
A = [param(1) param(2) param(3); 0 param(4) param(5); 0 0 1];
N_image = size(m,3);
f = [];
for i = 1:N_image
    lamda = 1/((norm(inv(A)*H(:,1,i))+norm(inv(A)*H(:,2,i)))/2);
    r1 = lamda*inv(A)*H(:,1,i);
    r2 = lamda*inv(A)*H(:,2,i);
    r3 = cross(r1, r2);
    R = [r1 r2 r3];
    [u,s,v] = svd(R);
    R = u*v';
    t = lamda*inv(A)*H(:,3,i);
    RT = [R(:,1) R(:,2) t];
    m_ = A*RT*M;
    m_ = [m_(1,:)./m_(3,:) ; m_(2,:)./m_(3,:); m_(3,:)./m_(3,:)];
    dm = m(:,:,i) - m_;
    f = [f dm(1,:) dm(2,:)];
end


function f = fun_1(H, m, M)
H0 = reshape(H,3,3);
H0 = H0';
m_ = H0*M;
m_ = [m_(1,:)./m_(3,:);m_(2,:)./m_(3,:);m_(3,:)./m_(3,:)];
d_m = m - m_;
f = [d_m(1,:) d_m(2,:)];
end
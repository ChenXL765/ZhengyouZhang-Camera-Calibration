function RT = Solve_Extrinsic(A,H)
lamda = (norm(inv(A)*H(:,1))+norm(inv(A)*H(:,2)))/2;
r1 = 1/lamda*inv(A)*H(:,1);
r2 = 1/lamda*inv(A)*H(:,2);
r3 = cross(r1, r2);
R = [r1 r2 r3];
[u,s,v] = svd(R);
R = u*v';
t = 1/lamda*inv(A)*H(:,3);
RT = [R(:,1) R(:,2) t];
end
clear,clc
syms deltax deltay deltaz dx dy dz real;
syms nx ox ax px ny oy ay py nz oz az pz pi real;
%ys010104
theta=sym('theta',[1,6]);
alpha=sym('alpha',[1,6]);
a=sym('a',[1,6]);
d=sym('d',[1,6]);
beta=sym('beta',[1,6]);
Dtheta=sym('Dtheta',[1,6]);
Dalpha=sym('Dalpha',[1,6]);
Da=sym('Da',[1,6]);
Dd=sym('Dd',[1,6]);
Dbeta=sym('Dbeta',[1,6]);

T=cell(1,6);
DT=cell(1,6);
Delta=cell(1,6);
J=cell(1,6);
JJ=cell(1,6);
JG=cell(1,6);
EI=cell(1,6);
E=cell(1,6);
G=cell(1,6);
RT=cell(1,6);
[ a1, a2, a3, a4, a5, a6]=deal(0,365 ,235 ,300 ,0,0 );
[ alpha1, alpha2, alpha3, alpha4, alpha5, alpha6]=deal(0, 0 ,0 , pi/2 ,pi/2,0);
[ d1, d2, d3, d4, d5, d6]=deal(500,d(2), -143.97, 0,0,100);
[theta1, theta2, theta3, theta4, theta5, theta6]=deal(theta(1), 0, theta(3), theta(4), theta(5), theta(6));%24 DH parameter
DD=[Da.',Dalpha.',Dtheta.',Dd.']; %24 DH parameter error setted for simulation
DD=[DD(1,1:4).';DD(2,1:4).';DD(3,1:4).';DD(4,1:4).';DD(5,1:4).';DD(6,1:4).'];
for i=1:6
     T{i}=sym(zeros(4,4));
     DT{i}=sym(zeros(4,4));
     Delta{i}=sym(zeros(4,4));
     J{i}=sym(zeros(6,6));
     JJ{i}=sym(zeros(6,6));
     EI{i}=sym(zeros(6,1));
     E{i}=sym(zeros(1,6));
     G{i}=sym(zeros(6,5));
     JG{i}=sym(zeros(6,5));
end
for i=1:6
    T{i}=subs(trotz(theta(i))*transl(0,0,d(i))*transl(a(i),0,0)*trotx(alpha(i)));%
end
for i=1:6
    GG1=diff(T{i},alpha(i));
    GG2=diff(T{i},theta(i));
    GG3=diff(T{i},a(i));
    GG4=diff(T{i},d(i));
 
    DT{i}=GG1*Dalpha(i)+GG2*Dtheta(i)+GG3*Da(i)+GG4*Dd(i);
    Delta{i}=T{i}\DT{i};
end
for i=1:6
    EI{i}(1:3,1)=Delta{i}(1:3,4);
    EI{i}(4,1)=Delta{i}(3,2);
    EI{i}(5,1)=Delta{i}(1,3);
    EI{i}(6,1)=Delta{i}(2,1);
     
end
M{1}= ((T{1}*T{2}*T{3}*T{4}*T{5}*T{6}));
M{2}= T{2}*T{3}*T{4}*T{5}*T{6};
M{3}= T{3}*T{4}*T{5}*T{6};
M{4}= T{4}*T{5}*T{6};
M{5}= T{5}*T{6};
M{6}= T{6};
J{7}=eye(6);
%simplify(M);
%T=[nx ox ax px;ny oy ay py;nz oz az pz;0 0 0 1];
%delta=[deltax;deltay;deltaz];
for i=1:6%J
    p=M{i}(1:3,4);
    nn=M{i}(1:3,1);
    o=M{i}(1:3,2);
    aa=M{i}(1:3,3);
    %B{i}=simplify(dot(delta,cross(p,n))+dot(n,D));
    J{i}(1:3,1:3)=(M{i}(1:3,1:3)).';
    J{i}(4:6,4:6)=(M{i}(1:3,1:3)).';
    J{i}(1,4:6)=cross(p,nn).';
    J{i}(2,4:6)=cross(p,o).';
    J{i}(3,4:6)=cross(p,aa).';
end
for i=1:6%J*G
    G{i}=[1 ,0 ,0 ,0;
         0,0 ,a(i)*cos(alpha(i)) ,sin(alpha(i));
          0 ,0, -a(i)*sin(alpha(i)) ,cos(alpha(i));
          0 ,1, 0 ,0;
          0, 0, sin(alpha(i)), 0;
          0 ,0,cos(alpha(i)),0];
   
end
EE=subs([J{2}*G{1},J{3}*G{2},J{4}*G{3},J{5}*G{4},J{6}*G{5},J{7}*G{6}]);
Rtheta=theta+Dtheta;
Ralpha=alpha+Dalpha;
Ra=a+Da;
Rd=d+Dd;
%Rbeta=beta+Dbeta;
[ Dbeta1, Dbeta2, Dbeta3, Dbeta4, Dbeta5, Dbeta6]=deal(0,0.00002,0,0 ,0,0);
[ Da1, Da2, Da3, Da4, Da5, Da6]=deal(1.7,1.2,-1.24,-1.61,1.77,1.13);
[ Dalpha1, Dalpha2, Dalpha3, Dalpha4, Dalpha5, Dalpha6]=deal(0.0066,-0.0033,-0.0077,-0.0018,0.0029,-0.0061);
%[ Dalpha1, Dalpha2, Dalpha3, Dalpha4, Dalpha5, Dalpha6]=deal(0.6,-0.3,-0.7,-0.18,0.9,-0.1);
[ Dtheta1, Dtheta2, Dtheta3, Dtheta4, Dtheta5, Dtheta6]=deal(0.0097,0.0066,-0.0014,-0.0053,0.0045,0.00188);
%[ Dtheta1, Dtheta2, Dtheta3, Dtheta4, Dtheta5, Dtheta6]=deal(0.9,0.6,-0.14,-0.5,0.4,0.1);
[ Dd1, Dd2, Dd3, Dd4, Dd5, Dd6]=deal(1.13,1.33,1.56,-2.4,2.13,-2.11);
%[ Dbeta1, Dbeta2, Dbeta3, Dbeta4, Dbeta5, Dbeta6]=deal(0,0.00002,0,0 ,0,0);
for i=1:6
    RT{i}=subs(trotz(Rtheta(i))*transl(0,0,Rd(i))*transl(Ra(i),0,0)*trotx(Ralpha(i)));
end
extraE=zeros(4,4);
extraE(1:3,4)=ones(3,1);
R= ((RT{1}*RT{2}*RT{3}*RT{4}*RT{5}*RT{6}))+extraE; 
n=50;% 50 random position for simulaiton
q=cell(1,n);
jgg=cell(1,n);
EER=cell(1,n);
ER=cell(1,n);
RR=cell(1,n);
g=-180 + (360).*rand(n,1);
h=200+ (1000-200).*rand(n,1);
c=-140 + (240).*rand(n,1);
ddd=-180 + (360).*rand(n,1);
l=-120 + (240).*rand(n,1);
f=-120 + (240).*rand(n,1);
yy=[g,h,c,ddd,l,f];
for i=1:n
    q{i}(1:6)=yy(i,1:6);
end
MMM=cell(1,n);
MM=M{1}(1:3,4);
for i=1:n 
    theta1=deg2rad(q{i}(1));
    d2=q{i}(2);
    %theta2=deg2rad(q{i}(2)-90);
    theta3=deg2rad(q{i}(3));
    theta4=deg2rad(q{i}(4));
    theta5=deg2rad(q{i}(5)+90);
    theta6=deg2rad(q{i}(6));
    K2=subs(R-M{1});
    ER{i}=K2(1:3,4);
    K3=subs(M{1});
    MMM{i}=K3(1:3,4);
     K1=K3\K2;
    %vpa(subs(R{1}));
    EER{i}(1:3,1)=K1(1:3,4);
    EER{i}(4,1)=K1(3,2);
    EER{i}(5,1)=K1(1,3);
    EER{i}(6,1)=K1(2,1);
    jgg{i}=subs(EE);
end
del=cell(1,n);
for i=1: length(q)
    del{i}=jgg{i}*subs(DD);
end
aaaaaa=cell(1,n);
for i=1: length(q)
    aaaaaa{i}=(del{i})-(EER{i});
end
jj=cat(1,jgg{:}); 
e=cat(1,EER{:}); 

figure(1)
eeeee=cat(2,ER{:});
x=eeeee(1,1:n);
y=eeeee(2,1:n);
z=eeeee(3,1:n);
plot3(x,y,z,'.');
title('補償前的誤差');
legend('補償前的誤差');
xlabel('x');
ylabel('y');
zlabel('z');
hold on


figure(2)%補償前的誤差（距離）
aaaa=zeros(1,n);
for i=1: length(q)
    aaaa(i)=sqrt(ER{i}(2).^2+ER{i}(1).^2+ER{i}(3).^2);
end

t=1:n;
plot(t,aaaa,'b+:');
title('補償前的誤差方差')
hold on
[aa2,b2]=max(aaaa);
plot(b2,aa2,'r*');
text(b2+0.1,aa2,num2str(aa2));
hold on
[aa2,b2]=min(aaaa);
plot(b2,aa2,'r*');
text(b2+0.1,aa2,num2str(aa2));
hold on
aa2=mean(aaaa);
line([1,n],[aa2,aa2])
text(0,aa2+0.1,num2str(aa2))
fprintf('Max=%f\n',max(aaaa))
fprintf('Min=%f\n',min(aaaa))
fprintf('Mean=%f\n',mean(aaaa))

figure(4)
aaaaa=zeros(1,n);
for i=1: n
    aaaaa(i)=sqrt(aaaaaa{i}(2).^2+aaaaaa{i}(1).^2+aaaaaa{i}(3).^2);
end
t=1:n;
plot(t,aaaaa,'.:');
title('JG*Delta-EER 方差');
hold on

figure(3)
eeeen=cat(2,MMM{:});
rx=eeeen(1,1:n);
ry=eeeen(2,1:n);
rz=eeeen(3,1:n);
figure(3)
plot3(rx,ry,rz,'g.');
title('取點理想位置');
hold on

figure(5)
eeeee=cat(2,ER{:});
x=eeeee(1,1:n);
y=eeeee(2,1:n);
z=eeeee(3,1:n);
plot3(x,y,z,'.');
title('補償前的誤差');
hold on
figure(5)
eeeee=cat(2,EER{:});
x=eeeee(1,1:n);
y=eeeee(2,1:n);
z=eeeee(3,1:n);
plot3(x,y,z,'g.');
title('補償前的誤差');
hold on

figure(5)
eeeee=cat(2,del{:});
x=eeeee(1,1:n);
y=eeeee(2,1:n);
z=eeeee(3,1:n);
plot3(x,y,z,'r.');
title('補償前的誤差');
hold on

figure(6)
aaaa=zeros(1,n);
for i=1: length(q)
    aaaa(i)=sqrt(EER{i}(2).^2+EER{i}(1).^2+EER{i}(3).^2);
end
t=1:n;
plot(t,aaaa,'b+:');
hold on
figure(6)
attt=zeros(1,n);
for i=1: length(q)
    attt(i)=sqrt(del{i}(2).^2+del{i}(1).^2+del{i}(3).^2);
end
t=1:n;
plot(t,attt,'g*');
legend('转换后的誤差方差','JG*Delta方差');
hold on

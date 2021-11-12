N = 10;
L = 1;
X = linspace(0,L,N+1);
Y = zeros(6*(N+1),1);
for i = 1:N+1
  Y(6*(i-1)+2) = X(i); 
end
for qwerty = 1:1
  E = 1;
r = 0.1;

I = (pi*r^4)/4;
nu = 0.3;


K = 1;0
G = E/(2*(1+nu));
delta = 0.1;
a = pi*r*r;

tol = 0.001;
F = zeros(6*N+6,1);
h = L/N;
P = qwerty*(3*E*I/L^3)*delta;
A = zeros(6*(N+1),6*(N+1));
A(1:3,1:3) = eye(3,3);
M = (1/2)*[0 1;-1 0];

while(true)
for i = 1:N 
  C = [(K*G*a) 0;0 (E*a)];
  ni = [Y(6*(i-1)+4) 0 Y(6*(i-1)+5)]; 
  ni1 = [Y(6*(i-1)+10) 0 Y(6*(i-1)+11)];
  ri = [Y(6*(i-1)+1) 0 Y(6*(i-1)+2)]; 
  ri1 = [Y(6*(i-1)+7) 0 Y(6*(i-1)+8)];
  thetai = Y(6*(i-1)+3);
  thetai1 = Y(6*(i-1)+9);
  mi = Y(6*(i-1)+6);
  mi1 = Y(6*(i-1)+12);
  R = [cos((thetai+thetai1)/2) sin((thetai+thetai1)/2);-sin((thetai+thetai1)/2) cos((thetai+thetai1)/2)];
  
  A(3+6*(i-1)+1:3+6*(i-1)+2,6*(i-1)+4:6*(i-1)+5)=-eye(2,2)/h;
  A(3+6*(i-1)+1:3+6*(i-1)+2,6*(i-1)+10:6*(i-1)+11)=eye(2,2)/h;
  
  A(3+6*(i-1)+3,6*(i-1)+1:6*(i-1)+2)=-cross(ni+ni1,[0 1 0])(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+4:6*(i-1)+5)=cross([0 1 0],ri1-ri)(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+7:6*(i-1)+8)=cross(ni+ni1,[0 1 0])(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+10:6*(i-1)+11)=cross([0 1 0],ri1-ri)(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+6) = -1/h;
  A(3+6*(i-1)+3,6*(i-1)+12) = 1/h;
  
  A(3+6*(i-1)+4:3+6*(i-1)+5,6*(i-1)+1:6*(i-1)+2) = (1/h)*R*C*R';
  A(3+6*(i-1)+4:3+6*(i-1)+5,6*(i-1)+7:6*(i-1)+8) = -(1/h)*R*C*R';
  Q = M*R*C*R';
  A(3+6*(i-1)+4:3+6*(i-1)+5,6*(i-1)+3)=-(1/(2*h))*(Q+Q')*((ri1-ri)(1,[1,3])')+(1/2)*M*R*C*[0 1]';
  A(3+6*(i-1)+4:3+6*(i-1)+5,6*(i-1)+9)=-(1/(2*h))*(Q+Q')*((ri1-ri)(1,[1,3])')+(1/2)*M*R*C*[0 1]';
  A(3+6*(i-1)+4:3+6*(i-1)+5,6*(i-1)+4:6*(i-1)+5) = (1/2)*eye(2,2);
  A(3+6*(i-1)+4:3+6*(i-1)+5,6*(i-1)+10:6*(i-1)+11) = (1/2)*eye(2,2);
  
  A(3+6*(i-1)+6,6*(i-1)+3) = E*I/h;
  A(3+6*(i-1)+6,6*(i-1)+6) = 1/2;
  A(3+6*(i-1)+6,6*(i-1)+9) = -E*I/h;
  A(3+6*(i-1)+6,6*(i-1)+12) = 1/2;
  
  F(3+6*(i-1)+1:3+6*(i-1)+2) = (1/h)*((ni1-ni)(1,[1 3])');
  F(3+6*(i-1)+3) = (1/h)*(mi1-mi)+(1/(2*h))*(dot([0 1 0],cross((ri1-ri),(ni+ni1))));
  F(3+6*(i-1)+6) = (mi+mi1)/2 - (E*I/h)*(thetai1-thetai);
  F(3+6*(i-1)+4:3+6*(i-1)+5) = (1/2)*((ni+ni1)(1,[1 3])')-R*C*(((R')/h)*((ri1-ri)(1,[1 3])')-[0;1]);
  n_n1 = [Y(6*(N)+4) Y(6*(N)+5)]';
theta_n1 = Y(6*N+3);
R_n1 = [cos(theta_n1) sin(theta_n1);-sin(theta_n1) cos(theta_n1)];
  F(3+6*(N)+1:6*N+6) = [n_n1-P*R_n1*[1 0]';mi1]';
  disp('ayush')
end
F(1:3) = [Y(1); Y(2); Y(3)];


A(3+6*(N-1)+7,6*(N-1)+10) = 1;
A(3+6*(N-1)+8,6*(N-1)+11) = 1;
A(3+6*(N-1)+7:3+6*(N-1)+8,6*(N)+3) = -P*M*R_n1*[1 0]';
A(3+6*(N-1)+9,6*(N-1)+12) = 1;
Yfinal = Y-A\F;
if(norm(F)<tol&&norm(Yfinal-Y)<tol)
break;
end
Y = Yfinal;
end
r1 = zeros(N+1,1);
r2 = zeros(N+1,1);
for i = 1:N+1
    r1(i) = Y(6*(i-1)+1); 
    r2(i) = Y(6*(i-1)+2);
endfor
plot(r2,r1,'k');
hold on
#axis([0 1 0 6]);
end
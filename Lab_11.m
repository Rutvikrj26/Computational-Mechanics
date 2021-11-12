N = 30;
E = 1;
r = 0.1;
L = 1;
I_hat = pi*(L)*(r^4)/4+(1/3)*pi*r^2*L^3;                              
rho = 1;
I = (pi*r^4)/4;
nu = 0.3;
Y = zeros(6*(N+1),1);
X = linspace(0,L,N+1);
alpha = 0.5;
beta = 0.5;
K = 1;
G = E/(2*(1+nu));
a = pi*r*r;
T = 200;
delta_t = 0.1;
r1 = zeros(N+1,(T/delta_t) + 1);
r1_dot = zeros(N+1,(T/delta_t) + 1);
r1_dotdot = zeros(N+1,(T/delta_t) + 1);
r3 = zeros(N+1,(T/delta_t) + 1);
r3_dot = zeros(N+1,(T/delta_t) + 1);
r3_dotdot = zeros(N+1,(T/delta_t) + 1);
theta = zeros(N+1,(T/delta_t) + 1);
theta_dot = zeros(N+1,(T/delta_t) + 1);
theta_dotdot = zeros(N+1,(T/delta_t) + 1);
for i = 1:N+1                                         
  Y(6*(i-1)+1) = (L/pi)*(1-cos((pi/N)*(i-1)));  
  r1(i,1) = Y(6*(i-1)+1);
  Y(6*(i-1)+3) = (pi/N)*(i-1); 
   theta(i,1) = Y(6*(i-1)+3);
  Y(6*(i-1)+2) = (L/pi)*(sin((pi/N)*(i-1)));
   r3(i,1) = Y(6*(i-1)+2);
   Y(6*(i-1)+6) = I*pi/L;
end
tol = 0.001;
F = zeros(6*N+6,1);
h = L/N;
A = zeros(6*(N+1),6*(N+1));
A(1:3,1:3) = eye(3,3);
M = [0 1;-1 0];
A(3+6*(N-1)+7,6*(N-1)+10) = 1;
A(3+6*(N-1)+8,6*(N-1)+11) = 1;
A(3+6*(N-1)+9,6*(N-1)+12) = 1;
for t = 1:T/delta_t+1
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
  A(3+6*(i-1)+1,6*(i-1)+1)=-a/(2*(beta*delta_t^2));
  A(3+6*(i-1)+1,6*(i-1)+7)=-a/(2*(beta*delta_t^2));
  A(3+6*(i-1)+2,6*(i-1)+2)=-a/(2*(beta*delta_t^2));
  A(3+6*(i-1)+2,6*(i-1)+8)=-a/(2*(beta*delta_t^2));
  
  A(3+6*(i-1)+3,6*(i-1)+1:6*(i-1)+2)=-cross(ni+ni1,[0 1 0])(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+4:6*(i-1)+5)=cross([0 1 0],ri1-ri)(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+7:6*(i-1)+8)=cross(ni+ni1,[0 1 0])(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+10:6*(i-1)+11)=cross([0 1 0],ri1-ri)(1,[1 3])/(2*h);
  A(3+6*(i-1)+3,6*(i-1)+6) = -1/h;
  A(3+6*(i-1)+3,6*(i-1)+12) = 1/h;
  A(3+6*(i-1)+3,6*(i-1)+3) = -I_hat/(2*(beta*delta_t^2));
  A(3+6*(i-1)+3,6*(i-1)+9) = -I_hat/(2*(beta*delta_t^2));
  
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
  
  F(3+6*(i-1)+1:3+6*(i-1)+2) = (1/h)*((ni1-ni)(1,[1 3])')-[(a/(2*beta))*((((-r1_dotdot(i,t)-r1_dotdot(i+1,t))/2)*(1-(2*beta)))+(2/delta_t^2)*(((ri(1)+ri1(1))/2)-((r1(i,t)+r1(i+1,t))/2)-((r1_dot(i,t)+r1_dot(i+1,t))/2)*delta_t));(a/(2*beta))*((((-r3_dotdot(i,t)-r3_dotdot(i+1,t))/2)*(1-(2*beta)))+(2/delta_t^2)*(((ri(3)+ri1(3))/2)-((r3(i,t)+r3(i+1,t))/2)-((r3_dot(i,t)+r3_dot(i+1,t))/2)*delta_t))];
  F(3+6*(i-1)+3) = (1/h)*(mi1-mi)+(1/(2*h))*(dot([0 1 0],cross((ri1-ri),(ni+ni1))))-(I_hat/(2*beta))*((((-theta_dotdot(i,t)-theta_dotdot(i+1,t))/2)*(1-(2*beta)))+(2/delta_t^2)*(((thetai+thetai1)/2)-((theta(i,t)+theta(i+1,t))/2)-((theta_dot(i,t)+theta_dot(i+1,t))/2)*delta_t));                                                                           
  F(3+6*(i-1)+6) = (mi+mi1)/2 - (E*I/h)*(thetai1-thetai);
  F(3+6*(i-1)+4:3+6*(i-1)+5) = (1/2)*((ni+ni1)(1,[1 3])')-R*C*(((R')/h)*((ri1-ri)(1,[1 3])')-[0;1]);
end
F(1:3) = [Y(1); Y(2); Y(3)];
F(3+6*(N)+1:6*N+6) = [ni1(1,[1,3]) mi1]';
Yfinal = Y-A\F;
if(norm(F)<tol && norm(Yfinal-Y)<tol)
break;
end
Y = Yfinal;
end
for j = 1:N+1
r1(j,t+1) = Y(6*(j-1)+1);
r1_dotdot(j,t+1) = (1/(2*beta))*((-r1_dotdot(j,t)*(1-(2*beta)))+(2/delta_t^2)*(r1(j,t+1)-r1(j,t)-r1_dot(j,t)*delta_t));
r1_dot(j,t+1) =  r1_dot(j,t)+delta_t*(alpha*r1_dotdot(j,t+1)+(1-alpha)*r1_dotdot(j,t));
theta(j,t+1) = Y(6*(j-1)+3);
theta_dotdot(j,t+1) = (1/(2*beta))*((-theta_dotdot(j,t)*(1-(2*beta)))+(2/delta_t^2)*(theta(j,t+1)-theta(j,t)-theta_dot(j,t)*delta_t)); 
theta_dot(j,t+1) = theta_dot(j,t)+delta_t*(alpha*theta_dotdot(j,t+1)+(1-alpha)*theta_dotdot(j,t));
r3(j,t+1) = Y(6*(j-1)+2);
r3_dotdot(j,t+1) = (1/(2*beta))*((-r3_dotdot(j,t)*(1-(2*beta)))+(2/delta_t^2)*(r3(j,t+1)-r3(j,t)-r3_dot(j,t)*delta_t));
r3_dot(j,t+1) =  r3_dot(j,t)+delta_t*(alpha*r3_dotdot(j,t+1)+(1-alpha)*r3_dotdot(j,t));
end
end
for t = 1:(T/delta_t)+1
  plot(r3(:,t),r1(:,t),'*');  
  axis([-.5 2 -1 1]);
  if(t==1)
  pause(2);
  end
  pause(0.0000001);
end

r = 0.01;
N = 100;
L = 1;
E = 1;
rho = 1;
w = 10;
I = (pi*r^4)/4;
beta = 0.1;
alpha = 0.7;
delta_x = L/N;
X = linspace(0,L,N+1);
T = input("Total time: ");
delta_t =0.01;
U_not = 1;
t = 0;
U = zeros(N+1,(T/delta_t)+1);
U_dot = zeros(N+1,(T/delta_t)+1);
U_dotdot = zeros(N+1,(T/delta_t)+1);
U(:,1) = 0;
U_dot(:,1) = 0;
A = zeros(N+1,N+1);
B = zeros(N+1,1);
for i = 2:N
   U_dotdot(i,1) = E*(U(i+1,1)+U(i-1,1)-2*U(i,1))/((delta_x^2)*rho);
end
for t = 1:T/delta_t
  B(1) = 0;
  B(N+1) = -U_not*(1-cos(w*t*delta_t));
  A(1,1) = 1;
  A(N+1,N+1) = 1;
  for j = 2:N
    B(j) = (rho/(2*beta))*((2/delta_t^2)*(U(j,t)+U_dot(j,t)*delta_t)+(1-2*beta)*U_dotdot(j,t));
    A(j,[j-1:j+1]) = [E/delta_x^2 -((2*E/delta_x^2)+(rho/(beta*delta_t^2))) E/delta_x^2];
  endfor
  U(:,t+1) = -A\(B);
  for i = 2:N
    U_dotdot(i,t+1) = (1/(2*beta))*((2*(U(i,t+1)-U(i,t)-U_dot(i,t)*delta_t)/delta_t^2)-(1-2*beta)*U_dotdot(i,t));
    U_dot(i,t+1) = U_dot(i,t) + delta_t*((1-alpha)*U_dotdot(i,t)+alpha*U_dotdot(i,t+1));  
  endfor
endfor
for i = 1:(T/delta_t)+1
plot(X,U(:,i));
axis([0 L -2*U_not 2*U_not]);
pause(0.001);
end
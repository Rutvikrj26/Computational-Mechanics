r = 0.04;
N = 10
L = 1
E = 1000
Khat = 0
K = 1
q = 0
I = (pi*r^4)/4;
a = pi*(r^2);
G = 1
Mhat = zeros(N,1);
Vhat = zeros(N,1);
A = zeros(2*N+4,2*N+4);
B = zeros(2*N+4,1);
h = L/N;
X = linspace(0,L,N+1);
Xb = zeros(N,1);
for i = 1:N
  Xb(i) = (X(i) + X(i+1))/2;
end
for i = 1:N
  Mhat(i) = (1/2)*q*(L-Xb(i))^2;
end
for i = 1:N
  Vhat(i) = q*(L - Xb(i));
endfor
A(1,1) = 1;
A(2,2) = 1;
for i = 1:N
  B(2*i+1:2*i+2) = [-Mhat(i); -Vhat(i)/(K*G*a)];
  A([2*i+1:2*i+2],[2*i-1:2*i+2]) = [0, -(E*I/h), 0, E*I/h; -1/h, -1/2, 1/h -1/2];
  A([2*i+1:2*i+2],[2*N+3:2*N+4]) = [-(L-Xb(i)), -1; -1/(K*G*a), 0];
endfor
A(2*N+3,2*N+4) = 1;
#A(2*N+4,2*N+1) = Khat;
A(2*N+4,2*N+3) = 1;
B(2*N+4) = -0.003;
Yhat = -A\B;
Yhat1 = zeros(N+1,1);
Yhat2 = zeros(N+1,1);
for i = 1:N+1
  Yhat1(i) = Yhat(2*i-1);
end
for i = 1:N+1
  Yhat2(i) = Yhat(2*i);
end
plot(X,Yhat1,'r');
hold on
#plot(X,Yhat2,'k');
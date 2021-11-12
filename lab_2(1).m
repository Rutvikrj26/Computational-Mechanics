r = input("Radius: ");
N = input("No. of elements in the beam: ");
L = input("Length of the beam: ");
E = input("Elastic modulus of the beam: ");
K = input("Spring constant: ");
q = input("What is the value of the distributed load: ");
I = (pi*r^4)/4;
Mhat = zeros(N,1);
A = zeros(N+3,N+3);
B = zeros(N+3,1);
h = L/N;
X = linspace(0,L,N+1);
Ydash = zeros(N+1,1);
for i = 1:N+1
  Ydash(i) = -(q/(E*I))*(((X(i)-L)^4)/24 + (L^3)*X(i)/6 - (L^4)/24);
end
for i = 2:N
  Mhat(i) = -(1/2)*q*(L-X(i))^2;
end
A(1,1) = 1; 
A(2,1) = -3/(2*h);
A(2,2) = (2/h);
A(2,3) = -1/(2*h);
for j = 3:N+1
  A(j,(j-2:j)) = [((E*I)/h^2) (-(2*E*I)/h^2) ((E*I)/h^2)];
  A(j,(N+2:N+3)) = [(-(L-X(j-1))) (-1)];
  B(j) = -Mhat(j-1);
endfor
A(N+2,(N+1:N+3)) = [K 1 0];
A(N+3,N+3) = 1;
Yhat = -A\B;
plot(X,Yhat(1:N+1));
hold on
plot(X,Ydash(1:N+1));

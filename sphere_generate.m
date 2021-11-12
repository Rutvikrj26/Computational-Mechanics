function ans = sphere_generate(rad,n)
a = linspace(0,pi,n);
b = linspace(0,2*pi,n);
for i = 1:length(a);
  x = rad*sin(a(i))*cos(b);
  y = rad*sin(a(i))*sin(b);
  z = rad*cos(a(i))*ones(length(b),1);
  plot3(x,y,z,'*-')
  hold on
end
end
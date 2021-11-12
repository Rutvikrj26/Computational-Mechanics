close all %closing all previous windows(plot)
x = linspace(1,10,100);
y = x;
z = x.^2;
plot(x,y,'ro-') %r for red, o for points, - for line
hold on
plot(x,z,'c*-')
figure()  %for new fig window
plot3(x,y,z)

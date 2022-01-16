x = -3 + 6*rand(50,1);
y = -3 + 6*rand(50,1);
v = sin(x).^4 .* cos(y);
F = scatteredInterpolant(x,y,v);
[xq,yq] = meshgrid(-3:0.1:3);
F.Method = 'linear';
vq1 = F(xq,yq);
figure(1), close(1), figure(1)
plot3(x,y,v,'mo')
hold on
mesh(xq,yq,vq1)
title('Nearest Neighbor')
legend('Sample Points','Interpolated Surface','Location','NorthWest')
view(2)
colormap de
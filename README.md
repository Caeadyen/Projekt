# Projekt
Modellierungspraktikum: Viskoelastizität
working setup:

para= Parameter(1,1,10,10,[0; 0], 1e20 , 0.3 , 3e10, 1e-3,-1e-3,1,10);
data = solvestoke(para);
plotgrid(data,para);
quiver(data.x_node[:,1],data.x_node[:,2],data.u[1:441],data.u[442:end]);


my = [1e21; 1e20; 1e23; 1e17; 1e10; 1e15; 1e17; 1e12; 1e10; 1e20]
para= Parameter(1,1,10,10,[0; 0], 10, my , 0.3 , 3e10, 1e-3,-1e-3,1,10);
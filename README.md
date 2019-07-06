# Projekt
Modellierungspraktikum: Viskoelastizität
working setup:

para= Parameter(1,1,10,10,[0; 0], 1e20 , 0.3 , 3e10, 1e-3,-1e-3,1,10);
data = solvestoke(para);
plotgrid(data,para);
quiver(data.x_node[:,1],data.x_node[:,2],data.u[1:441],data.u[442:end]);

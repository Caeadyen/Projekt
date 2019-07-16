# Projekt
Modellierungspraktikum: Viskoelastizität
working setup:

Einheitswürfel zusammengedrückt
20x20 Gitter, 1 Jahr zeitschritt, 10 jahre Gesammt

my = [1e20];
para= Parameter(1,1,20,20,[0; 0], 1, my , 0.25 , 3e10, 1e10, 0.01,-0.01,1,10);
data = builddata(para);
data = solvestoke(data,para);



beliebiges Gebiet mit verschiedenen Viskositäten:
100x50 gebiet, 50x25 Gitter, 5 verschiedene Viskositäten, 10 Jahre schritte, 2000 Jahre gesamt
my = [1e23; 1e20; 1e15; 1e17; 1e20] ;
para= Parameter(100,50,50,25,[0; 0], 5, my , 0.25 , 3e10, 1e10, 0.001,-0.001,10,2000);
data = builddata(para);
data = solvestoke(data,para);

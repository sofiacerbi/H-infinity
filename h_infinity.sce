clear
s=poly(0,'s')

////////////////////////////////////Planta//////////////////////////////////////
//Definición de parámetros
Mc = 1
Mp = 0.25
Ip = 7.88*10e-3
lp = 0.33
ng = 1
Kg = 3.7
nm = 1
Kt = 0.00767
Rm = 2.6
rmp = 6.35*10e-3
Beq = 5.4
Km = 0.00767
g = 9.81
Bp = 0.0024

//Constantes definidas
betaa = (Mc+Mp)*Ip + Mc*Mp*lp*lp;
gammaa = ng*Kg*nm*Kt/(Rm*rmp*rmp);

//Matriz A
A22 = -(Ip + Mp*lp*lp)*Beq/betaa - (Ip+Mp*lp*lp)*Kg*Km*gammaa/betaa;
A23 = Mp*Mp*lp*lp*g/betaa;
A24 = Mp*lp*Bp/betaa;
A42 = (Mp*lp*Beq + Mp*lp*Kg*Km*gammaa)/betaa;
A43 = -(Mc+Mp)*Mp*g*lp/betaa;
A44 = -(Mc+Mp)*Bp/betaa;

A = [0,1,0,0; 0,A22,A23,A24; 0,0,0,1; 0,A42,A43,A44];

//Matriz B
B2 = (Ip + Mp*lp*lp)*rmp*gammaa/betaa;
B4 = -Mp*lp*rmp*gammaa/betaa;

B = [0; B2; 0; B4];

//Matriz C
//Si se toma como salida la posición del carrito
C1 = [1, 0, 0, 0];
//Si se toma como salida el ángulo del péndulo
C2 = [0, 0, 1, 0];

//Matriz D: se asume cero
D = [0.001];

//Modelo de la planta
planta=syslin('c',A,B,C1,D);

////////////////////////Controlador H-infinito//////////////////////////////////
//Se usa el mixed-sensitivity approach

//A = error de seguimiento mínimo en estado estable
//ωB = ancho de banda mínimo
//M = magnitud máxima de S
a=0.0002; m=2; wb=100;

//Selección de pesos w
w1_n = ((1/m)*s+wb);
w1_d = (s+wb*a);
w1 = syslin('c',w1_n,w1_d);

w2_d = 2*(s/1000 + 1);
w2_n = s/300+1;
w2 = syslin('c',w2_n,w2_d);

//Se crea una planta generalizada
[Ap0,Bp0,Cp0,Dp0]=abcd(planta);
[Aw1,Bw1,Cw1,Dw1]=abcd(w1);
[Aw2,Bw2,Cw2,Dw2]=abcd(w2);

Ap = [Ap0 zeros(size(Ap0,1),size(Aw1,2)) zeros(size(Ap0,1),size(Aw2,2)); -Bw1*Cp0 Aw1 zeros(size(Aw1,1),size(Aw2,2)); Bw2*Cp0 zeros(size(Aw2,1),size(Aw1,2)) Aw2];
Bp = [zeros(size(Bp0,1),size(Bw1,2)) Bp0; Bw1 zeros(size(Bw1,1),size(Bp0,2)); zeros(size(Bw2,1),size(Bw1,2)) zeros(size(Bw2,1),size(Bp0,2))];
Cp = [-Dw1*Cp0 Cw1 zeros(size(Cw1,1),size(Cw2,2)); Dw2*Cp0 zeros(size(Cw2,1),size(Cw1,2)) Cw2; -Cp0 zeros(size(Cp0,1),size(Cw1,2)) zeros(size(Cp0,1),size(Cw2,2))];
Dp = [Dw1 0.001; 0 0; 1 0];  //D12 para que la matriz sea de rango completo

p_gen = syslin('c',Ap,Bp,Cp,Dp);

//Se sintetiza el controlador
Sk = ccontrg(p_gen,[1,1],5);
[Ak,Bk,Ck,Dk] = abcd(Sk)

//Se cierra el lazo
L = Sk*planta;
S = 1/(1+L); 
T = 1-S;
[Acl,Bcl,Ccl,Dcl] = abcd(T);


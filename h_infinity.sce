

function F=f(z,x,l,phi,theta)
    
endfunction

//Definición de parámetros
mt = 10;
mb = 10;
mc = 10;
g = 9.81;

//Condiciones iniciales
//se asumen cero
y0=0;
t0=0;
t=0:0.1:10;

//Se crea la EDO
y = ode(y0,t0,t,f);
plot(t,y)
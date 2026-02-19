# Seminario-de-modelado
function dx = Pendulo(~,y,z)

xc = y(1);
xcd = y(2);
a = y(3);
ad = y(4);

Ip  = z(1);
Mc  = z(2);
lp  = z(3);
Mp  = z(4);
Fc  = z(5);
Beq = z(6);
Bp  = z(7);
g   = z(8);

den = (Mc + Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(a)^2;

xdd = ( ...
    (Ip + Mp*lp^2)*Fc ...
  + Mp^2*lp^2*g*cos(a)*sin(a) ...
  - (Ip + Mp*lp^2)*Beq*xcd ...
  - (Ip*Mp*lp + Mp^2*lp^3)*ad^2*sin(a) ...
  - Mp*lp*a*cos(a)*Bp ...
  ) / den;

add = ( ...
    (Mc + Mp)*Mp*g*lp*sin(a) ...
  - (Mc + Mp)*Bp*ad ...
  + Fc*Mp*lp*cos(a) ...
  - Mp^2*lp^2*ad^2*sin(a)*cos(a) ...
  - Beq*Mp*lp*xcd*cos(a) ...
  ) / den;

dx = zeros(4,1);
dx(1) = xcd;
dx(2) = xdd;
dx(3) = ad;
dx(4) = add;
end
clear; clc; close all

Ip  = 0.0079;        
Mc  = 0.7031;        
lp  = 0.3302;        
Mp  = 0.23;          
Fc  = 0;             
Beq = 4.3;           
Bp  = 0.0024;        
g   = 9.81;          

params = [Ip Mc lp Mp Fc Beq Bp g];

a0 = deg2rad(1);   

y0 = [ ...
    0;        
    0;        
    a0;   
    0];       

tspan = [0 10];

[t, y0] = ode45(@(t,y) Pendulo(t,y,params), tspan, y0);

figure
subplot(2,1,1)
plot(t, y0(:,3),'LineWidth',1.5)
grid on
ylabel('\alpha (rad)')
title('Ángulo del péndulo')

subplot(2,1,2)
plot(t, y0(:,1),'LineWidth',1.5)
grid on
ylabel('xc (m)')
xlabel('Tiempo (s)')
title('Posición')

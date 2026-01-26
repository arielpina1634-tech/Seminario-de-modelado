# Seminario-de-modelado
/*Bloque 1*/
function xdd = xddot(xc_dot, alpha, alpha_dot)

Ip  = 0.0079;
Mc  = 0.7031;
lp  = 0.3302;
Mp  = 0.23;
Fc  = 0;
Beq = 4.3;
g   = 9.81;
Bp  = 0.0024;

D = (Mc + Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(alpha)^2;

N = (Ip + Mp*lp^2)*Fc ...
    + Mp^2*lp^2*g*cos(alpha)*sin(alpha) ...
    - (Ip + Mp*lp^2)*Beq*xc_dot ...
    - (Ip*Mp*lp - Mp^2*lp^3)*alpha_dot^2*sin(alpha) ...
    - Mp*lp*alpha_dot*cos(alpha)*Bp;

xdd = N / D;
end

/*Bloque 2*/
function alphadd = alphaddot(xc_dot, alpha, alpha_dot)

Ip  = 0.0079;
Mc  = 0.7031;
lp  = 0.3302;
Mp  = 0.23;
Fc  = 0;
Beq = 4.3;
g   = 9.81;
Bp  = 0.0024;

D = (Mc + Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(alpha)^2;

N = (Mc + Mp)*Mp*g*lp*sin(alpha) ...
    - (Mc + Mp)*Bp*alpha_dot ...
    + Fc*Mp*lp*cos(alpha) ...
    - Mp^2*lp^2*alpha_dot^2*sin(alpha)*cos(alpha) ...
    - Beq*Mp*lp*xc_dot*cos(alpha);

alphadd = N / D;
end

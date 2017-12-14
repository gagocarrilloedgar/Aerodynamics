
clear
clear workspace;
clc

for alpha= 0:10
    
M=100;    
N=M+1;
R=[cos(0) sin(0); -sin(0) cos(0)];
f=0.02;
p=0.4;
xh=10;
etha=0;

    
for i = 1:N
        x = (1/2)*(1-cos(((i-1)*pi)/(N-1)));
        X(i) = x;
        X_des(i)=x;
        
        if (x < xh)
            if (x <= p)
                z = (f/(p^2))*((2*p*x)-(x^2));
                z_des=z;
            else
                z = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
                z_des=z;
            end
        else
             zh=(f/((1-p)^2))*(1-(2*p)+(2*p*xh)-(xh^2));
             z = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
             x_rot = xh + R(1,1)*(x-xh) + R(1,2)*(z-zh);
             z_rot = zh + (R(2,1)*(x-xh) + R(2,2)*(z-zh));
             
             z_des= (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2))-(tan(etha)*(x-xh));
             
             x=x_rot;
             z=z_rot;
        end
        Z(i) = z;
        Z_des(i)=z_des;
        
end

% C�lculo de los parametros que definen cada panel 

% Vector normal, tangencial i m�dulo.

for i = 1:(N-1)
    C(i) = sqrt((X(i+1)-X(i))^2+(Z(i+1)-Z(i))^2);
    n(i,1) = -(Z(i+1)-Z(i))/C(i);
    n(i,2) = (X(i+1)-X(i))/C(i);
    t(i,1) = (X(i+1)-X(i))/C(i);
    t(i,2) = (Z(i+1)-Z(i))/C(i);
end 

% C�lculo d control points, punto de aplicaci�n del torbellino y sistema
% A*CIRCULACI�N = RHS
for i = 1:(N-1)
    xcp(i) = X(i) + (X(i+1)-X(i))*3/4;
    zcp(i) = Z(i) + (Z(i+1)-Z(i))*3/4;
    x0(i) = X(i) + (X(i+1)-X(i))*1/4;
    z0(i) = Z(i) + (Z(i+1)-Z(i))*1/4;
    
end
for i=1:(N-1)
    for j = 1:(N-1)
        r(i,j)= (xcp(i)-x0(j)).^2 + (zcp(i)-z0(j)).^2;
        u(i,j)= (zcp(i)-z0(j))/(2*pi*r(i,j));
        w(i,j)= -(xcp(i)-x0(j))/(2*pi*r(i,j));
        A(i,j) = u(i,j)*n(i,1) + w(i,j)*n(i,2);
    end
    RHS(i) = -(cos(deg2rad(alpha))*n(i,1) + sin(deg2rad(alpha)*n(i,2)));
end 


% C�lculo de la Circulaci�n, el Cl y el Cm
Circulacion =linsolve(A,RHS');

sumaCl = 0;
Cm = 0;
sumaCm = 0;

Cl=2*sum(Circulacion);
for j = 1:(N-1)
    Cm = -2*Circulacion(j)*X(j) + sumaCm;
    sumaCm = Cm;
end
Cm0=Cm+(1/4)*Cl;

Clvector(alpha+1)=Cl;
Cmvector(alpha+1)=Cm;
Cm0vector(alpha+1)=Cm0;


%CL,CM,CI ES el calculado a trav�s de la TAT.
end
for alpha= -10:0
    
M=100;    
N=M+1;
R=[cos(0) sin(0); -sin(0) cos(0)];
f=0.02;
p=0.4;
xh=10;
etha=0;

    
for i = 1:N
        x = (1/2)*(1-cos(((i-1)*pi)/(N-1)));
        X(i) = x;
        X_des(i)=x;
        
        if (x < xh)
            if (x <= p)
                z = (f/(p^2))*((2*p*x)-(x^2));
                z_des=z;
            else
                z = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
                z_des=z;
            end
        else
             zh=(f/((1-p)^2))*(1-(2*p)+(2*p*xh)-(xh^2));
             z = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
             x_rot = xh + R(1,1)*(x-xh) + R(1,2)*(z-zh);
             z_rot = zh + (R(2,1)*(x-xh) + R(2,2)*(z-zh));
             
             z_des= (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2))-(tan(etha)*(x-xh));
             
             x=x_rot;
             z=z_rot;
        end
        Z(i) = z;
        Z_des(i)=z_des;
        
end

% C�lculo de los parametros que definen cada panel 

% Vector normal, tangencial i m�dulo.

for i = 1:(N-1)
    C(i) = sqrt((X(i+1)-X(i))^2+(Z(i+1)-Z(i))^2);
    n(i,1) = -(Z(i+1)-Z(i))/C(i);
    n(i,2) = (X(i+1)-X(i))/C(i);
    t(i,1) = (X(i+1)-X(i))/C(i);
    t(i,2) = (Z(i+1)-Z(i))/C(i);
end 

% C�lculo d control points, punto de aplicaci�n del torbellino y sistema
% A*CIRCULACI�N = RHS
for i = 1:(N-1)
    xcp(i) = X(i) + (X(i+1)-X(i))*3/4;
    zcp(i) = Z(i) + (Z(i+1)-Z(i))*3/4;
    x0(i) = X(i) + (X(i+1)-X(i))*1/4;
    z0(i) = Z(i) + (Z(i+1)-Z(i))*1/4;
    
end
for i=1:(N-1)
    for j = 1:(N-1)
        r(i,j)= (xcp(i)-x0(j)).^2 + (zcp(i)-z0(j)).^2;
        u(i,j)= (zcp(i)-z0(j))/(2*pi*r(i,j));
        w(i,j)= -(xcp(i)-x0(j))/(2*pi*r(i,j));
        A(i,j) = u(i,j)*n(i,1) + w(i,j)*n(i,2);
    end
    RHS(i) = -(cos(deg2rad(alpha))*n(i,1) + sin(deg2rad(alpha)*n(i,2)));
end 


% C�lculo de la Circulaci�n, el Cl y el Cm
Circulacion =linsolve(A,RHS');

sumaCl = 0;
Cm = 0;
sumaCm = 0;

Cl=2*sum(Circulacion);
for j = 1:(N-1)
    Cm = -2*Circulacion(j)*X(j) + sumaCm;
    sumaCm = Cm;
end
Cm0=Cm+(1/4)*Cl;

Clvector1(alpha+11)=Cl;
Cmvector1(alpha+11)=Cm;
Cm0vector1(alpha+11)=Cm0;


%CL,CM,CI ES el calculado a trav�s de la TAT.
end

CL=zeros(22,1);
CM=zeros(22,1);
CM0=zeros(22,1);

for k=1:11
    CL(k)=Clvector1(1,k);
    CM(k)=Cmvector1(1,k);
    CM0(k)=Cm0vector1(1,k);
end

for k=1:11
    CL(k+10)=Clvector(1,k);
    CM(k+10)=Cmvector(1,k);
    CM0(k+10)=Cm0vector(1,k);
end

CL1=zeros(21,1);
CM1=zeros(21,1);
CM01=zeros(21,1);

for k=1:21
    CL1(k)=CL(k);
    CM1(k)=CM(k);
    CM01(k)=CM0(k);
end

Clalpha=((CL1(21)-CL1(20)))*(180/pi);
alpha0=(degtorad(10)-(CL1(21)/Clalpha))*(180/pi);
cm0= sum(CM01)/20;

disp('Cl alpha is: ');
disp(Clalpha);
disp('Cm0 is: ');
disp(cm0);
disp('zero lift alpha is: ');
disp(alpha0);

x=(-10:10);
plot(x,CL1,'b');
hold on;
plot(x,CM1,'r');
plot(x,CM01,'g--');

xlabel('Angle of attack \alpha (º)');
yyaxis right;
ylim([-0.3,0.1]);
ylabel('Pitching moment coefficient','interpreter','latex');

yyaxis left;
ylabel ('Lift coefficient','interpreter','latex');
ylim([-0.9,1.5]);

grid on;

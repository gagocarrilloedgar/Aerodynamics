
clear
clear workspace;
clc

%% Matrices y variables del sistema

% Geometry discretization


contador=1;
M=250;
Clvector=zeros(250,1);
Cmvector=zeros(250,1);
errorCl=zeros(250,1);
errorCm=zeros(250,1);
CL=0.6664;
CM=-0.219731;

for G=1:M
N=G+1;
X = zeros(N,1);
Z = zeros(N,1);
f=0.02;
p=0.4;
alpha=deg2rad(4);
C = zeros((N-1),1);
t = zeros((N-1),2);
n = zeros((N-1),2);
xcp = zeros(N-1,1);
zcp = zeros(N-1,1);
x0 = zeros(N-1,1);
z0 = zeros(N-1,1);
u = zeros(N-1,N-1);
w = zeros(N-1,N-1);
r = zeros(N-1,N-1);
A = zeros(N-1,N-1);
RHS = zeros(N-1,1);

  
for i = 1:N
        x = (1/2)*(1-cos(((i-1)*pi)/(N-1)));
        X(i)=x;
        
%Perfil no sim�trico                
            if (x <= p)
                Z(i) = (f/(p^2))*((2*p*x)-(x^2));
                
            else
                Z(i) = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
            end
            
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
    RHS(i) = -(cos(alpha)*n(i,1) + sin(alpha)*n(i,2));
end 


% C�lculo de la Circulaci�n, el Cl y el Cm
Circulacion = linsolve(A,RHS);


sumaCl=0;
sumaCm=0;


for j = 1:(N-1)
    Cl = 2*Circulacion(j) + sumaCl;
    sumaCl = Cl;
    Cm = -2*Circulacion(j)*X(j)*cos(alpha) + sumaCm;
    sumaCm = Cm;
end

Clvector(contador,1) = Cl;
Cmvector(contador,1)= Cm;

errorCl(contador,1)=abs((Cl-CL)/CL)*100;
errorCm(contador,1)=abs((Cm-CM)/CM)*100;

contador=contador+1;
%CL,CM,CI ES el calculado a trav�s de la TAT.
end

%% plots

%% Plot CL
figure
x=(1:250);
plot(x,Clvector,'r');
title('Cl - error','interpreter','latex');
xlabel('Number of panels M','interpreter','latex');
yyaxis right;
ylabel ('error','interpreter','latex');
hold on;
t=0:1:250;
plot(x,errorCl,'b');
yyaxis left;
ylabel ('Cl ','interpreter','latex');
plot(t,ones(size(t))*0.6664,'g')


hold off;

%% PLOT CM
figure

x=(1:250);
plot(x,Cmvector,'r');
xlabel('Number of panels M','interpreter','latex');
title('$ Cm_{LE} $ - error','interpreter','latex');
hold on;
yyaxis right; 
ylabel ('error','interpreter','latex');
plot(x,errorCm,'b');
yyaxis left;
ylabel ('$ Cm_{LE} $','interpreter','latex');
plot(t,ones(size(t))*-0.210786,'g')
hold off;

figure
plot(X,Z);






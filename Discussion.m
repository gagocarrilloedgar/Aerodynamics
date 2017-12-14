
function[]= Discussion()

clear
clear workspace;
clc
for f=1:6
    
    for p=1:6
    
M=100;  
N=M+1;
etha=0;
R=[cos(deg2rad(0)) sin(deg2rad(0)); -sin(deg2rad(0)) cos(deg2rad(0))];
alpha=deg2rad(4);
xh=10;  

f=f/100;
p=p/10;

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
             
             z_des= (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2))-(tan(deg2rad(etha))*(x-xh));
             
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
    RHS(i) = -(cos(alpha)*n(i,1) + sin(alpha)*n(i,2));
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

f=f*100;
p=p*10;

Clvector(f,p)=Cl;
Cmvector(f,p)=Cm;
Cm0vector(f,p)=Cm0;

CL=0.6762;
CM=-0.049;
CM0=0.12005;

    end   
end

%% plot CL
CL1=zeros(1,6);
CL2=zeros(1,6);
CL3=zeros(1,6);
CL4=zeros(1,6);
CL5=zeros(1,6);
CL6=zeros(1,6);

for k=1:6
    
    CL1(k)=Clvector(1,k);
    CL2(k)=Clvector(2,k);
    CL3(k)=Clvector(3,k);
    CL4(k)=Clvector(4,k);
    CL5(k)=Clvector(5,k);
    CL6(k)=Clvector(6,k);
    
end

for k=1:6
    CL01(k)= (deg2rad(4)-(CL1(k)/(2*pi)))*(180/pi);
    CL02(k)= (deg2rad(4)-(CL2(k)/(2*pi)))*(180/pi);
    CL03(k)= (deg2rad(4)-(CL3(k)/(2*pi)))*(180/pi);
    CL04(k)= (deg2rad(4)-(CL4(k)/(2*pi)))*(180/pi);
    CL05(k)= (deg2rad(4)-(CL5(k)/(2*pi)))*(180/pi);
    CL06(k)= (deg2rad(4)-(CL6(k)/(2*pi)))*(180/pi);
end




figure ;
x=(1:6);
plot(x,CL01,'DisplayName','f=0.01');
hold on;
plot(x,CL02,'DisplayName','f=0.02');
plot(x,CL03,'DisplayName','f=0.03');
plot(x,CL04,'DisplayName','f=0.04');
plot(x,CL05,'DisplayName','f=0.05');
plot(x,CL06,'DisplayName','f=0.06');

title('$\alpha_{l0}$ function of(f,p)','interpreter','latex');
ylabel('$\alpha_{l0}^{(o)} $ ','interpreter','latex');
xlabel('P','interpreter','latex');
grid on;
legend('Location','southwest');
hold off;


%% plot Cm
CM1=zeros(1,6);
CM2=zeros(1,6);
CM3=zeros(1,6);
CM4=zeros(1,6);
CM5=zeros(1,6);
CM6=zeros(1,6);

for k=1:6
    CM1(k)=Cmvector(1,k)+(1/4)*CL1(k);
    CM2(k)=Cmvector(2,k)+(1/4)*CL2(k);
    CM3(k)=Cmvector(3,k)+(1/4)*CL3(k);
    CM4(k)=Cmvector(4,k)+(1/4)*CL4(k);
    CM5(k)=Cmvector(5,k)+(1/4)*CL5(k);
    CM6(k)=Cmvector(6,k)+(1/4)*CL6(k);
end

figure ;
x=(1:6);
plot(x,CM1,'DisplayName','f=0.01');
hold on;
plot(x,CM2,'DisplayName','f=0.02');
plot(x,CM3,'DisplayName','f=0.03');
plot(x,CM4,'DisplayName','f=0.04');
plot(x,CM5,'DisplayName','f=0.05');
plot(x,CM6,'DisplayName','f=0.06');

title('$CM_0$ function of(f,p)','interpreter','latex');
ylabel('$CM_0$','interpreter','latex');
xlabel('P','interpreter','latex');
grid on;
legend('Location','southwest');
hold off;

end

clear workspace;
clear 
clc
Cl=zeros(4,20);

for xh=4:9
    
for etha=1:20
    
M=100;  
N=M+1;
R=[cos(deg2rad(etha)) sin(deg2rad(etha)); -sin(deg2rad(etha)) cos(deg2rad(etha))];
f=0.02;
p=0.4;
alpha=deg2rad(4);
xh=xh/10; 
CL=0.6664;
CM=-0.219731;



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

xh=xh*10;
Cl(xh,etha)=2*sum(Circulacion);
end

xhvector(xh)=xh;   

end

    CL1=zeros(1,20);
    CL2=zeros(1,20);
    CL3=zeros(1,20);
    CL4=zeros(1,20);
  
    
for k=1:20
    CL1(k)=Cl(6,k);
    CL2(k)=Cl(7,k);
    CL3(k)=Cl(8,k);
    CL4(k)=Cl(9,k);
end

alphaTAT=deg2rad(4)-(CL/(2*pi));

for k=1:20
    alpha01(k)=alphaTAT-(deg2rad(4)-(CL1(k)/(6.1609)));
    alpha02(k)=alphaTAT-(deg2rad(4)-(CL2(k)/(6.1609)));
    alpha03(k)=alphaTAT-(deg2rad(4)-(CL3(k)/(6.1609)));
    alpha04(k)=alphaTAT-(deg2rad(4)-(CL4(k)/(6.1609)));

end

    
EF=zeros(1,4);
    EF(4)=(alpha01(20)-alpha03(1))/deg2rad(20);
    EF(3)=(alpha02(20)-alpha04(1))/deg2rad(20);
    EF(2)=(alpha03(20)-alpha04(1))/deg2rad(20);
    EF(1)=(alpha04(20)-alpha04(1))/deg2rad(20);

%% zero lift alpha corregido 

EFC=zeros(1,4);


for i=1:4
    
    EFC(i)=EF(i)*0.7;
end




%% plot
figure;
x=(1:20);
plot(x,alpha01,'DisplayName','xh=0.6');
hold on;
plot(x,alpha02,'DisplayName','xh=0.7');
plot(x,alpha03,'DisplayName','xh=0.8');
plot(x,alpha04,'DisplayName','xh=0.9');

title(' $ \Delta \alpha_{l_0} $ zero lift alpha as function of eta-hinge position','interpreter','latex');
xlabel('\eta (º)');
ylabel('$ \Delta \alpha_{l0} $ ','interpreter','latex');

grid on;
legend('Location','northwest');
legend('show');
hold off;

figure
y=(0.1:0.1:0.4);
plot(y,EF,'DisplayName','Theorical');
hold on;
plot(y,EFC,'DisplayName','Corrected');

title('Flap efficency factor ','interpreter','latex');
xlabel(' $ flap-chord-ratio $','interpreter','latex');
ylabel('$\frac{\Delta \alpha_{l0}}{\Delta \eta} $ ','interpreter','latex');

legend('Location','southeast');
legend('show');
hold off;
grid on;


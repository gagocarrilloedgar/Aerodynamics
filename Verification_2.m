
clear workspace;
clear 
clc


inc=0.1;
a=(6/inc)+1;
Cl=zeros(a,20);
alpha0=zeros(a,20);
EF=zeros(a,1);
EFC=zeros(a,1);
CL=0.6664;
CM=-0.219731;
contador=1;
flapcord=zeros(a,1);




for g=(4:inc:10)
    
for etha=1:20
    
M=100;  
N=M+1;
R=[cos(deg2rad(etha)) sin(deg2rad(etha)); -sin(deg2rad(etha)) cos(deg2rad(etha))];
alpha=deg2rad(4);
xh=g/10; 
X = zeros(N,1);
Z = zeros(N,1);
f=0.02;
p=0.4;
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
        X(i) = x;
       
         if(x<xh)
            if (x <= p)
                Z(i) = (f/(p^2))*((2*p*x)-(x^2));
                
            else
                Z(i) = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
            end
            
          else
             z = (f/((1-p)^2))*(1-(2*p)+(2*p*x)-(x^2));
             zh=(f/((1-p)^2))*(1-(2*p)+(2*p*xh)-(xh^2));
             X(i) = xh + R(1,1)*(x-xh) + R(1,2)*(z-zh);
             Z(i) = zh + R(2,1)*(x-xh) + R(2,2)*(z-zh);
                    
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

Circulacion =linsolve(A,RHS);

Cl(contador,etha)=2*sum(Circulacion);
end
contador=contador+1;

end


for h=1:a
for k=1:20
    alpha0(h,k)=-0.0359-(deg2rad(4)-(Cl(h,k)/(6.2830)));
end
end
%Incremento del angulo de sustentacion nula (con falp - sinflap) es un
%incremento negativo pero lo contabilizamos como positivo


for i=1:a

    EF(i,1)=(alpha0(i,20)-alpha0(i,1))/deg2rad(19);
    
end

for i =1:a
    
    EFC(i,1)=EF(i,1)*0.7;
end 
for i=1:a
    flapcord(i,1)=0.6-((inc*(i-1))/10);
end



%% plot
figure;
x=(1:20);
plot(x,alpha0);
hold on

title(' $ \Delta \alpha_{l_0} $ zero lift alpha as function of eta-hinge position','interpreter','latex');
xlabel('\eta (�)');
ylabel('$ \Delta \alpha_{l0} $ ','interpreter','latex');

grid on;
legend('Location','northwest');
legend('show');
hold off;

figure
plot(flapcord,EF,'DisplayName','Theoretical');
hold on;

plot(flapcord,EFC,'DisplayName','Corrected');

title('Flap efficency factor ','interpreter','latex');
xlabel('flap-chord ratio','interpreter','latex');
ylabel('$\frac{\Delta \alpha_{l0}}{\Delta \eta} $ ','interpreter','latex');

legend('Location','southeast');
legend('show');
hold off;
grid on;



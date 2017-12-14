clear
clear workspace;
clc

%% Inputs del sistema 
% Se controlan los inputs del sistema para que el programa cumpla con todas las suposiciones 
% 4-digit NACA Airfoil
inputairfoil = 0; 
while (inputairfoil == 0)
    airfoil = input('4-digit NACA: ','s');
    if (length(airfoil) ~= 4)
        disp('Invalid NACA airfoil, must be a 4-digit NACA series');
        fprintf('\n');
        inputairfoil = 0;
    else
        if (str2double(airfoil(1)) == 0)
            if (str2double(airfoil(2)) ~= 0)
                disp('Invalid NACA airfoil');
                fprintf('\n');
                inputairfoil = 0;
            else
                % Camber line
                f = str2double(airfoil(1))/100; %Max camber
                p = str2double(airfoil(2))/10; %Max camber line position
                % Airfoil thickness
                t = str2double(airfoil(3:4))/100;
                if (t > 0.15)
                    disp('Thin Airfoil Theory may not apply for the chosen airfoil thickness');
                    fprintf('\n');
                    inputairfoil = 0;
                else
                    inputairfoil = 1;
                end
            end
        else
            % Camber line
            f = str2double(airfoil(1))/100; %Max camber
            p = str2double(airfoil(2))/10; %Max camber line position
            % Airfoil thickness
            t = str2double(airfoil(3:4))/100;
            if (t > 0.15)
                disp('Thin Airfoil Theory may not apply for the chosen airfoil thickness');
                fprintf('\n');
                inputairfoil = 0;
            else
                inputairfoil = 1;
            end
        end
    end
end 

% Angle of attack
inputalpha = 0;
while (inputalpha == 0)
    alpha = input('Angle of attack (�): ');
    if (alpha < -10 || alpha > 10 || alpha == j)
        disp('Thin Airfoil Theory may not apply for the chosen angle of attack');
        fprintf('\n');
        inputalpha = 0;
    else
        inputalpha = 1;
    end
end
alpha = deg2rad(alpha);   

% Number of panels
M = input('Number of panels: ');

% Geometry discretization
inputdist = 0;
while (inputdist == 0)
    dist = input( '(a) Uniform\n(b) Full cosine\nChoose between these different geometry discretizations: ','s');
    if (dist == 'A' || dist == 'a')
        inputdist = 1;
    elseif (dist == 'B' || dist == 'b')
        inputdist = 1;
    else
        disp('Invalid geometry discretization');
        fprintf('\n');
        inputdist = 0;
    end
end

% Flap geometry
inputflap = 0;
while (inputflap == 0)
    flap = input ('Include a flap ? [Y/N] ','s');
    if (flap == 'Y' || flap == 'y')
        % Hige location
        inputhige = 0;
        while (inputhige == 0)
            xh = input('Flap hinge location (in tenths of chord (0-10)): ');
            xh = xh/10;
            if (xh > 1)
                disp('Invalid hinge position');
                fprintf('\n');
                inputhige = 0;
            else
                % Flap deflection angle
                etha = input('Flap deflection angle(�): ');
                etha = deg2rad(etha);
                inputhige = 1;
            end
        end
        inputflap = 1;
    elseif (flap == 'N' || flap == 'n')
        xh = 1;
        etha = 0;
        inputflap = 1;
    else
        disp('Invalid input flap');
        fprintf('\n');
        inputflap = 0;
    end
end

%% Matrices y variables del sistema

N = M + 1;
X = zeros(N,1);
Z = zeros(N,1);
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
R=[cos(etha) sin(etha); -sin(etha) cos(etha)];

% Geometry discretization
for i = 1:N
    
    if (dist == 'A' || dist == 'a')
        x = (i-1)/M;
        X(i)=x;
%Perfil sim�trico  
        if (f == 0 && p == 0)
           if (x < xh)
              Z(i) = 0;
           else
              X(i)= xh + R(1,1)*(x-xh); 
              Z(i) = R(2,1)*(x-xh);
            end
%Perfil no sim�trico                
        else
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
       
       
        
    elseif (dist == 'B' || dist == 'b')
        x = (1/2)*(1-cos(((i-1)*pi)/(N-1)));
        X(i)=x;
        
%Perfil sim�trico  
        if (f == 0 && p == 0)
           if (x < xh)
              Z(i) = 0;
           else
              X(i)= xh + R(1,1)*(x-xh); 
              Z(i) = R(2,1)*(x-xh);
           end
%Perfil no sim�trico                
        else
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
    RHS(i) = -(cos(alpha)*n(i,1) + sin(alpha)*n(i,2));
end
for i = 1:(N-1)
    for j = 1:(N-1)
        r(i,j)= (xcp(i)-x0(j)).^2 + (zcp(i)-z0(j)).^2;
        u(i,j)= (zcp(i)-z0(j))/(2*pi*r(i,j));
        w(i,j)= -(xcp(i)-x0(j))/(2*pi*r(i,j));
        A(i,j) = u(i,j)*n(i,1) + w(i,j)*n(i,2);
    end
end

% C�lculo de la Circulaci�n, el Cl y el Cm
Circulacion = linsolve(A,RHS);

Cl = 0;
sumaCl = 0;
Cm = 0;
sumaCm = 0;


for j = 1:(N-1)
    Cl = 2*Circulacion(j) + sumaCl;
    sumaCl = Cl;
    Cm = -2*Circulacion(j)*X(j)*cos(alpha) + sumaCm;
    sumaCm = Cm;
end

%% Plots
figure
ax1=subplot(2,1,1);
plot(ax1,X,Z,'c*');
xlim([0,1]);
ylim([-0.2,0.1]);
title('Mean camber line discretization','interpreter','latex');
xlabel('x/c');
ylabel('z');
ax2=subplot(2,1,2);
X1=linspace(0,1,N-1);
plot(ax2,X1,Circulacion,'r*');
xlim([0,1]);
title('Circulation Distribution','interpreter','latex');
xlabel('x/c');
ylabel('$$\Gamma$$','Interpreter','latex');









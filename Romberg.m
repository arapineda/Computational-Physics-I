%% Romberg para funciones
% Este c�digo calcula la integral num�rica con el m�todo de Romberg
% implementado a la regla de trapecios para un intervalo y una funci�n
% dada.

clear

a = 0; % L�mite inferior
b = 10; % L�mite superior
n = 10; % N�mero final de intervalos
f = @(x) cos(x);
in = 1;
t = zeros(1,n);

% Trapecios m�ltiples intervalos

for j = 1:1:n
h = (b-a)/in;
xi = 0;
xk = 0;
res = 0;
    for i = 1:1:in % Calcula integral para diferentes cantidades de intervalos (in)
    
    xi = a + (i-1)*h;
    xk = a + i*h;
    preres = (h/2)*(f(xi)+f(xk));
    res = res + preres;
    
    end
    t(j) = res; % Asigna valores con diferentes intervalos a un vector
    in = in*2; % Duplica cantidad de intervalos
end

fprintf("La integral con m�todo de trapecios es: %4.4f\n", t(end));

% Romberg
vR = t;
for k = 1:1:n-1 % Nuevo vector
    for i = 1:1:n-1 % C�lculos con elementos del vector anterior
        vR(i)= (4^k * t(i+1)-(t(i))) /((4^k)-1); % Vector que guarda mejoras con Romberg
    end
    t = vR; % Nuevo vector t para c�lculo iterativo
end

fprintf("La integral con m�todo de Romberg es: %4.4f", t(end));

%% Romberg para puntos

clear

a = 0; % L�mite inferior
b = 10; % L�mite superior
n = 10; % N�mero final de intervalos
vIT = zeros(1,n);
in = 2;

% Trapecios m�ltiple intervalos
for i=1:1:n
    x = linspace(a,b,in);
    y = sin(x);
    v = (sum(y)-y(1)-y(in));
    vIT(i) = ((b-a)/(2*(in-1)))*(y(1)+2*v+y(in));
    in = in*2 - 1;
end
trapecio = vIT(end);
fprintf("Valor de la integral con m�todo de trapecios: %4.4f \n",trapecio);

% Metodo de Romberg
vR=vIT;
for k=1:1:n-1
    for i=1:1:n-1
        vR(i)= (4^k * vIT(i+1)-(vIT(i))) /((4^k)-1); % Vector que guarda mejoras con Romberg
    end
    vIT=vR; % Nuevo vector vIT para c�lculo iterativo
end
rom = vIT(1);
fprintf("Valor de la integral con m�todo de Romberg: %4.4f\n",rom);

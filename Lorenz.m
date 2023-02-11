% Este programa calcula numéricamente el resultado del sistema de
% ecuaciones de Lorenz después de cierto tiempo utilizando los métodos de
% Euler, Runge-Kutta, ODE45 y ODE113. Se asume que ODE45 obtiene la
% solución correcta por lo que se realiza la comparación con el método de
% Euler y el método de Runge-Kutta.
% En este programa, el usuario puede definir los parametros del sistema de
% ecuaciones, las condiciones iniciales, el intervalo de tiempo en el que se evalua y
% el diferencial de tiempo.

clear
close all
clc

% Se definen los parametros meteorológicos de Lorenz para el sistema de
% ecuaciones.
r=28;        % rho
s=10;        % sigma
b=8/3;       % beta

% Se definen las condiciones iniciales arbitrariamente
x0=10;
y0=11;
z0=12;

% Se define el intervalo de tiempo en el que se evaluará la ecuación y el
% diferencial de tiempo.
t0=0;              % Tiempo inicial
tf=10;             % Tiempo final
dt=0.0001;        % Diferencial de tiempo
t=t0:dt:tf;       % Intervalo de tiempo

fprintf('Valor final en t=%4.4f \n ',tf)

% Método de Euler: 
    % Se crea un vector para x, y, z los cuales definen la la solución al
    % sistema de ecuaciones de Lorenz en un tiempo (t) utilizando el método de
    % Euler. Finalmente se grafican las soluciones y se despliega la posición
    % final.
    
    % Se crean los vectores x,y,z comenzando en t0 utilizando las
    % condiciones iniciales
    x=x0;
    y=y0;
    z=z0;

    for it=1:length(t) % Se crean los vetores x, y, z en un tiempo (t>t0)
        x(it+1)=x(it)+(s * ( y(it) - x(it) ) )* dt; % Se actualiza el valor de x utilizando la ecuación diferencial de Lorenz de x en el tiempo.
        y(it+1)=y(it)+(x(it) * (r -z(it) ) - y(it)) *dt; % Se actualiza el valor de y utilizando la ecuación diferencial de Lorenz de y en el tiempo.
        z(it+1)=z(it)+( x(it) * y(it) - b * z(it) ) *dt; % Se actualiza el valor de z utilizando la ecuación diferencial de Lorenz de z en el tiempo.
    end
    
    % Se grafica la evolucion de las soluciones de x, y, z en el intervalo
    % de tiempo definido
    
    figure(1)
    subplot(2,2,1)
    set(gcf,'position',[1,1,900,1000])
    plot3(x,y,z,'w', 'LineWidth', 1.5)
    title('Método de Euler: Soluciones en el tiempo')
    xlabel('x'); ylabel('y'); zlabel('z')
    set(gca,'Color','k')

    % Se despliegan los valores finales de x, y, z despues de un tiempo (t)
    % para el método de Euler.
    fprintf('- Método de Euler: \n   ( %4.4f , %4.4f , %4.4f ) \n',x(end), y(end), z(end))
    
    %Se guardan los valores de x, y, z para futuros calculos
    xEuler=x(1:length(t));
    yEuler=y(1:length(t));
    zEuler=z(1:length(t));

%Runge-Kutta
    % Se crea un vector para x, y, z los cuales definen la la solución al
    % sistema de ecuaciones de Lorenz en un tiempo (t) utilizando el método de
    % Runge-Kutta de orden 4. Finalmente se grafican las soluciones y se despliega la posición
    % final.
    
    % Se crean los vectores x,y,z comenzando en t0 utilizando las
    % condiciones iniciales
    
    x=x0;
    y=y0;
    z=z0;

    fx = @(x,y) (s * ( y - x ) ); % Se define la ecuación diferencial de Lorenz para x.
    fy = @(x,y,z) (x * (r -z) - y); % Se define la ecuación diferencial de Lorenz para y.
    fz = @(x,y,z) ( x * y - b * z ); % Se define la ecuación diferencial de Lorenz para z.

    for i=1:(length(t)-1) % Se crean los vetores x, y, z en un tiempo (t)
       k1 = dt * fx (x(i),y(i));        % k1, k2, k3, k4 son las constantes que permiten definir x en t+dt
       m1 = dt * fy (x(i),y(i),z(i));   % m1, m2, m3, m4 son las constantes que permiten definir x en t+dt
       n1 = dt * fz (x(i),y(i),z(i));   % n1, n2, n3, n4 son las constantes que permiten definir x en t+dt
       k2 = dt * fx (x(i)+k1/2, y(i)+m1/2);  
       m2 = dt * fy (x(i)+k1/2, y(i)+m1/2, z(i)+n1/2);  
       n2 = dt * fz (x(i)+k1/2, y(i)+m1/2, z(i)+n1/2);  
       k3 = dt * fx (x(i)+k2/2, y(i)+m2/2);
       m3 = dt * fy (x(i)+k2/2, y(i)+m2/2, z(i)+n2/2);  
       n3 = dt * fz (x(i)+k2/2, y(i)+m2/2, z(i)+n2/2);
       k4 = dt * fx (x(i)+k3, y(i)+m3);
       m4 = dt * fy (x(i)+k3, y(i)+m3, z(i)+n3);
       n4 = dt * fz (x(i)+k3, y(i)+m3, z(i)+n3);
       
       x(i+1) = x(i) + (k1+(k2*2)+(k3*2)+k4)/6; % Se actualiza el valor de x utilizando las constantes previamente calculadas.
       y(i+1) = y(i) + (m1+(m2*2)+(m3*2)+m4)/6; % Se actualiza el valor de y utilizando las constantes previamente calculadas.
       z(i+1) = z(i) + (n1+(n2*2)+(n3*2)+n4)/6; % Se actualiza el valor de z utilizando las constantes previamente calculadas.
    end

    % Se grafica la evolucion de las soluciones de x, y, z en el intervalo
    % de tiempo definido
    
    subplot(2,2,2)
    plot3(x,y,z,'w', 'LineWidth', 1.5)
    title('Método Runge-Kutta: Soluciones en el tiempo')
    xlabel('x'); ylabel('y'); zlabel('z')
    set(gca,'Color','k')

    % Se despliegan los valores finales de x, y, z despues de un tiempo (t)
    % para el método de Runge-Kutta.
    fprintf('- Método de Runge-Kutta: \n   ( %4.4f , %4.4f , %4.4f ) \n',x(end), y(end), z(end))
    
    %Se guardan los valores de x, y, z para futuros calculos
    xRungeKutta=x(1:length(t));
    yRungeKutta=y(1:length(t));
    zRungeKutta=z(1:length(t));
%ODE45
    % Se crea una matriz (rlorenz) para x, y, z los cuales definen la la solución al
    % sistema de ecuaciones de Lorenz en un tiempo (t) utilizando el  comando de Matlab ODE45. Finalmente se grafican las soluciones y se despliega la posición
    % final.
    Parametros = [s; r; b]; % Se crea un vector con los parametros previamente definidos
    r0 = [x0; y0; z0];      % Se crea un vector con las condiciones iniciales previamente definidos
    options =odeset('RelTol',1e-12,'AbsTol', 1e-12*ones(1,3));  % Se establece una tolerancia relativa y absoluta para conseguir resultados más precisos.
    [tflorenz,rlorenz] = ode45( @(tflorenz,rflorenz)florenzz(tflorenz,rflorenz,Parametros),t,r0,options); % Se resuelve el sistema de ecuaciones diferenciales y como resultado una matriz con los valores de x,y,z en el tiempo

    % Se grafica la evolucion de las soluciones de x, y, z en el intervalo
    % de tiempo definido
    subplot(2,2,3)
    plot3( rlorenz(:,1),rlorenz(:,2),rlorenz(:,3),'w','LineWidth',1.5);
    title('Método ODE45: Soluciones en el tiempo')
    xlabel('x'); ylabel('y'); zlabel('z')
    set(gca,'Color','k')

    % Se despliegan los valores finales de x, y, z despues de un tiempo (t)
    % para el método de ODE45.
    fprintf('- Método ODE45: \n   ( %4.4f , %4.4f , %4.4f ) \n',rlorenz(end,1), rlorenz(end,2), rlorenz(end,3))
    
    %Se guardan los valores de x, y, z para futuros calculos
    xODE45=rlorenz((1:length(t)),1)';
    yODE45=rlorenz((1:length(t)),2)';
    zODE45=rlorenz((1:length(t)),3)';
%ODE113
    Parametros = [s; r; b]; % Se crea un vector con los parametros previamente definidos
    r0 = [x0; y0; z0];      % Se crea un vector con las condiciones iniciales previamente definidos
    options =odeset('RelTol',1e-12,'AbsTol', 1e-12*ones(1,3));  % Se establece una tolerancia relativa y absoluta para conseguir resultados más precisos.
    [tflorenz,rlorenz] = ode113( @(tflorenz,rflorenz)florenzz(tflorenz,rflorenz,Parametros),t,r0,options); % Se resuelve el sistema de ecuaciones diferenciales y como resultado una matriz con los valores de x,y,z en el tiempo

        % Se grafica la evolucion de las soluciones de x, y, z en el intervalo
        % de tiempo definido
        subplot(2,2,4)
        plot3( rlorenz(:,1),rlorenz(:,2),rlorenz(:,3),'w','LineWidth',1.5);
        title('Método ODE113: Soluciones en el tiempo')
        xlabel('x'); ylabel('y'); zlabel('z')
        set(gca,'Color','k')

        % Se despliegan los valores finales de x, y, z despues de un tiempo (t)
        % para el método de ODE45.
        fprintf('- Método ODE113: \n   ( %4.4f , %4.4f , %4.4f ) \n',rlorenz(end,1), rlorenz(end,2), rlorenz(end,3))
        
        %Se guardan los valores de x, y, z para futuros calculos
        xODE113=rlorenz((1:length(t)),1)';
        yODE113=rlorenz((1:length(t)),2)';
        zODE113=rlorenz((1:length(t)),3)';


 
% Errores
    % Error Euler vs ODE45
        % Se calcula el error del método de Euler en cada dimensión en
        % comparación con el ODE45
        xEE=abs((xEuler-xODE45)./xODE45);
        yEE=abs((yEuler-yODE45)./yODE45);
        zEE=abs((zEuler-zODE45)./zODE45);
        %Se grafican los errores
         figure (2)
         set(gcf,'position',[500,10,1000,1000])
         subplot(3,3,1)
            plot(xEE)
            title('Método Euler vs ODE 45: Error en X')
            xlabel('x')
            ylabel('% de Error')
        subplot(3,3,2)
            plot(yEE)
            title('Método Euler vs ODE 45: Error en Y')
            xlabel('y')
            ylabel('% de Error')
        subplot(3,3,3)
            plot(zEE)
            title('Método Euler vs ODE 45: Error en Z')
            xlabel('z')
            ylabel('% de Error')
        
    % Error Runge-Kutta vs ODE45
        % Se calcula el error del método de Runge-Kutta en cada dimensión en
        % comparación con el ODE45
        xERK=abs((xRungeKutta-xODE45)./xODE45);
        yERK=abs((yRungeKutta-yODE45)./yODE45);
        zERK=abs((zRungeKutta-zODE45)./zODE45);
        % Se grafican los errores
         subplot(3,3,4)
            plot(xERK)
            title('Método Runge-Kutta vs ODE 45: Error en X')
            xlabel('x')
            ylabel('% de Error')
        subplot(3,3,5)
            plot(yERK)
            title('Método Runge-Kutta vs ODE 45: Error en Y')
            xlabel('y')
            ylabel('% de Error')
        subplot(3,3,6)
            plot(zERK)
            title('Método Runge-Kutta vs ODE 45: Error en Z')
            xlabel('z')
            ylabel('% de Error')
    % Error ODE113 vs ODE45
        % Se calcula el error del método de OD113 en cada dimensión en
        % comparación con el ODE45
        xEODE113=abs((xODE113-xODE45)./xODE45);
        yEODE113=abs((yODE113-yODE45)./yODE45);
        zEODE113=abs((zODE113-zODE45)./zODE45);
        % Se grafican los errores
         subplot(3,3,7)
            plot(xEODE113)
            title('Método ODE 113 vs ODE 45: Error en X')
            xlabel('x')
            ylabel('% de Error')
        subplot(3,3,8)
            plot(yEODE113)
            title('Método ODE 113 vs ODE 45: Error en Y')
            xlabel('y')
            ylabel('% de Error')
        subplot(3,3,9)
            plot(zEODE113)
            title('Método ODE 113 vs ODE 45: Error en Z')
            xlabel('z')
            ylabel('% de Error')
        
% Comparación Gráfica de Euler, Runge-Kutta y 113 con ODE45
    % Euler vs ODE45
        figure(3)
        plot3(xEuler,yEuler,zEuler,'b--', xODE45, yODE45, zODE45,'r')
        title('Gráfica del método de Euler vs ODE45')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        legend('Euler','ODE45')
        set(gcf,'position',[1000,1,1000,1000])
    % Runge Kutta vs ODE45
        figure(4)
        plot3(xRungeKutta,yRungeKutta,zRungeKutta, 'b--', xODE45, yODE45, zODE45,'r')
        title('Gráfica del método de Runge-Kutta vs ODE45')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        legend('Runge-Kutta','ODE45')
        set(gcf,'position',[1,1,1000,1000])
    % ODE113 vs ODE45
        figure(5)
        plot3(xODE113,yODE113,zODE113, 'b--', xODE45, yODE45, zODE45,'r')
        title('Gráfica del método de ODE 113a vs ODE45')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        legend('ODE 113','ODE45')
        set(gcf,'position',[1,1,1000,1000])
    
        
% Funcion florenzz que se utiliza para calcular ODE45 y ODE113
  function dr= florenzz(tflorenz,rflorenz,Parametros)
dr= [
        Parametros(1) * ( rflorenz(2) - rflorenz(1) );
        rflorenz(1)*(Parametros(2)-rflorenz(3))-rflorenz(2);
        rflorenz(1)*rflorenz(2) - Parametros(3)*rflorenz(3);
        ];
  end

clear
clc

forg = @(var) (var^2); %Función a evaluar
a = 0; % Límite inferior de la integral
b = 6; % Límite superior de la integral

% Se define el orden de la serie del polinomio de Legendre
iol= 1; % Orden desde el cual se comienza a visualizar
iou= 10; % Orden hasta el cual se termina de visualizar
        
syms x;
xn = ((b-a)/2*x+(b+a)/2); % Se define nuevo dominio para cambiar límites de integración a [-1,1]
f = forg(xn); % Función con límites redefinidos de -1 a 1 a los nuevos (a,b)
fpl = forg(x); % Función usada para graficar
l=10000;
error=1;
    
 for in=iol:iou
     fg=f; % Se guarda la función a evaluar
     syms x
     ax = gca;
    % Polinomio de Legendre
        PL=legendreP(in,x);
        % Gráfica Legendre orden in
            subplot(3,1,1);
            fplot(PL)
            set(gcf,'position',[100,10,1000,1000])
            title(['Gráfica polinomio de Legendre Orden ' num2str(in)])
            xlim([-1,1])
            xlabel('x')
            ylim([-1,1])
            ylabel('P_l(x)')
            legend('P(x)')
            ax = gca;
            ax.XAxisLocation = 'origin';
            ay.YAxisLocation = 'origin';
        % Raices de Legendre
            rPL = vpasolve(legendreP(in,x) == 0);
            rPL = double(rPL');
%       % Se obtiene el valor de la derivada del polinomio de Legendre en
%       las raices del polinomio de Legendre
            dPL=1;
            for id=1:in
                difPL=diff(PL,x,1);
                x=rPL(1,id);
                dPL(1,id)=double(subs(difPL));
                syms x
            end
        % Pesos de la nueva función
            w=rPL(1,:);
            for iw=1:length(w)
                w(1,iw)= 2 / ( (1-(rPL (1,iw) )^2)* (( dPL(1,iw) )^2) );
            end
    % Serie de Taylor de la función a evaluar
        syms x;
        tpl=taylor(f, x, 'Order', in*2-1);
        t = tpl*((b-a)/2);
         % Gráfica de la función y la serie de Taylor
            subplot(3,1,2);
            fplot([f,tpl])
            title(['Aproximación de la serie de Taylor Grado ' num2str(in*2-1)])
            legend({'F(x)','Serie de Taylor de F(x)'})
            xlim([a-1,b+1])
            xlabel('x')
           ylim([-1,3])
            ylabel('F(x)')
            ax = gca;
            ax.XAxisLocation = 'origin';
            ay.YAxisLocation = 'origin';
            sum=0;
    % Se calcula el valor de la integral utilizando las raices del polinomio de Legendre 
        for is=1:in
            x=rPL(1,is);
            fl=subs(t);
            sum=sum+fl*w(1,is);    
        end
           
        %Se calcula el error por cada iteración y se grafica.
            subplot(3,1,3);
            errorx(1,error)=in;
            errory(1,error)= abs ((int(fpl, a,b) - sum) / (int(fpl, a,b)));
            hold on
            plot(errorx,errory, '-o')
            title('Error')
            axis([iol,iou,0, inf])
            set(gca, 'XTick', iol:iou)
            xlabel('Orden del polinomio de Legendre')
            ylabel('Error')
            hold on
            error=error+1;
            pause(0.05)

        f=fg; % Se recupera la función a evaluar
        
 end

sum=double(sum);
fprintf('El valor de la integral desde %4.3f hasta %4.3f es %4.3f',a,b,sum)
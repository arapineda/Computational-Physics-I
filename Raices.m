% Raices: se encuentra la x donde f(x)=0.
% Se despliega una matriz X donde el primer renglón es el valor de x y el
% segundo renglón es la función evaluada donde 

x1=linspace(0,2*pi,1000); % el vector x donde se define la función a evaluar
y1=cos(x1); % La función donde se busca el cero. Si se busca un minimo o un máximo, esta función debe de ser la derivada, en este caso, se pueden quitar los comentarios de F para evaluar dicha función donde se encuentra un minimo y un máximo.
%F=F % Función a maximizar o minimizar.
[h,l]=size(y1);
i=1;
x=1;
ix=1;
while i<l;
    if y1(i)>0;
        while y1(i)>0 & i<l;
            i=i+1;
        end
        if i~=l
            x(1,ix)=x1(i);
            %x(2,ix)=F(i);
            ix=ix+1;
        end
    else y1(i)<0;
        while y1(i)<0 & i<l;
            i=i+1;
        end
        if i~=l
            x(1,ix)=x1(i);
            %x(2,ix)=F(i);
            ix=ix+1;
        end
    end
end
x % Se dexpliega x